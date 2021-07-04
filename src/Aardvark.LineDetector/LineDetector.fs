namespace Aardvark.LineDetection


open Aardvark.Base
open Microsoft.FSharp.NativeInterop

#nowarn "9"

[<Struct>]
type LineDetectionConfig =
    {
        Threshold           : float //= 0.4
        GrowThreshold       : float //= 0.2
        Tolerance           : float //= 3.5
        MinLength           : float //= 4.0
        MinStability        : float //= 0.8
        MaxGuessCount       : int   //= 20
        MergeRadius         : float //= 20.0
        MinQuality          : float //= 0.8
        Threads             : int
    }

    static member Default =
        {
            Threshold           = 0.4
            GrowThreshold       = 0.1
            Tolerance           = 3.5
            MinLength           = 10.0
            MinStability        = 0.8
            MaxGuessCount       = 20
            MergeRadius         = 20.0
            MinQuality          = 0.8
            Threads             = 0
        }

module LineDetector =
    open System.Threading

    module private MatrixInfo =
        let iter2 (threads : int) (ct : CancellationToken) (action : int64 -> int64 -> unit) (a : MatrixInfo) (b : MatrixInfo) =
            let action = OptimizedClosures.FSharpFunc<_,_,_>.Adapt action
            let threads = if threads <= 0 then System.Environment.ProcessorCount else threads

            if threads > 1 then
                let blockSize = V2i(128, 128)
                let blockCount = V2d a.Size / V2d blockSize |> ceil |> V2i

                let blocks =
                    seq {
                        for by in 0 .. blockCount.Y - 1 do
                            for bx in 0 .. blockCount.X - 1 do
                                let o = V2i(bx,by) * blockSize
                                let s = 
                                    let e = min (o + blockSize) (V2i a.Size)
                                    e - o

                                yield a.SubMatrix(o, s), b.SubMatrix(o, s)

                    }

                let action =
                    System.Action<_>(fun (ab : MatrixInfo, bb : MatrixInfo) ->
                        ab.ForeachIndex(bb, fun ai bi ->
                            action.Invoke(ai, bi)
                        )
                    )

                let o = System.Threading.Tasks.ParallelOptions(MaxDegreeOfParallelism = threads, CancellationToken = ct)
                System.Threading.Tasks.Parallel.ForEach(blocks, o, action) |> ignore
            else
                a.ForeachIndex(b, fun ai bi ->
                    action.Invoke(ai, bi)
                    if ai &&& 1023L = 0L then ct.ThrowIfCancellationRequested()
                )

    let private edgeDetectLaplace (threads : int) (ct : CancellationToken) (image : PixImage<byte>) =
        let edges = PixImage<byte>(Col.Format.Gray, image.Size)

        NativeMatrix.using (edges.GetChannel 0L) (fun pEdges ->
            NativeVolume.using (image.Volume) (fun pSrc ->

                let offsets =
                    [|
                        struct( 8, 0L)
                        struct(-1, pSrc.DX)
                        struct(-1,-pSrc.DX)
                        struct(-1, pSrc.DY)
                        struct(-1,-pSrc.DY)
                        struct(-1, pSrc.DX + pSrc.DY)
                        struct(-1, pSrc.DX - pSrc.DY)
                        struct(-1,-pSrc.DX + pSrc.DY)
                        struct(-1,-pSrc.DX - pSrc.DY)
                    |]

                let edgeInfo = pEdges.Info.SubMatrix(V2l.II, pEdges.Size - 2L*V2l.II)
                let srcInfo = pSrc.Info.SubXYMatrix(0L).SubMatrix(V2l.II, pEdges.Size - 2L*V2l.II)

                let readGrayValue =
                    match image.Format with
                    | Col.Format.Gray -> 
                        fun (ptr : nativeptr<byte>) -> NativePtr.read ptr
                    | Col.Format.RGB | Col.Format.RGBA ->
                        let go = int pSrc.DZ
                        let bo = int pSrc.DZ <<< 1
                        fun (ptr : nativeptr<byte>) -> C3b(NativePtr.read ptr, NativePtr.get ptr go, NativePtr.get ptr bo).ToGrayByte()
                            
                    | Col.Format.BGR | Col.Format.BGRA ->
                        let go = int pSrc.DZ
                        let ro = int pSrc.DZ <<< 1
                        fun (ptr : nativeptr<byte>) -> C3b(NativePtr.get ptr ro, NativePtr.get ptr go, NativePtr.read ptr).ToGrayByte()
                    | fmt ->
                        failwithf "bad format: %A" fmt

                (edgeInfo, srcInfo) ||> MatrixInfo.iter2 threads ct (fun (ei : int64) (si : int64) ->
                //edgeInfo.ForeachIndex(srcInfo, fun (ei : int64) (si : int64) ->
                        
                    let mutable maxValue = 0uy
                    let mutable oo = 0L
                    for i in 0L .. pSrc.SZ - 1L do
                        let src = NativePtr.add pSrc.Pointer (int (si + oo))
                        let mutable sum = 0
                        for struct(w, o) in offsets do
                            let v = int (NativePtr.get src (int o))
                            sum <- sum + w*v
                        let value = sum |> clamp 0 255 |> byte

                        if value > maxValue then maxValue <- value
                        oo <- oo + pSrc.DZ

                    NativePtr.set pEdges.Pointer (int ei) maxValue

                )
            )
        )

        edges
       
    let private neighbours =
        [|
            for x in -1 .. 1 do
                for y in -1 .. 1 do
                    if x <> 0 || y <> 0 then yield V2i(x,y)
        |]
        
    let private cmp = System.Func<_,_,_>(fun struct(a,_) struct(b,_) -> System.Collections.Generic.Comparer<float>.Default.Compare(b,a))

    let private detectLinesPar (config : LineDetectionConfig) (ct : CancellationToken) (image : PixImage<byte>) =  
        let threads = if config.Threads <= 0 then System.Environment.ProcessorCount else config.Threads
        
        let inline iter (size : V2i) (action : V2i -> unit) =
            if threads = 1 then
                for y in 0 .. size.Y - 1 do
                    for x in 0 .. size.X - 1 do
                        action (V2i(x,y))
                    ct.ThrowIfCancellationRequested()
            else
                let coords = Dither.getCoords size
                let o = System.Threading.Tasks.ParallelOptions(MaxDegreeOfParallelism = threads, CancellationToken = ct)
                System.Threading.Tasks.Parallel.ForEach(coords, o, action) |> ignore

        //Log.startTimed "detecting edges"
        let edges = edgeDetectLaplace threads ct image
        //Log.stop()

        let lines = System.Collections.Concurrent.ConcurrentBag()
        //Log.startTimed "finding lines"
        NativeMatrix.using (edges.GetChannel 0L) (fun pEdges ->
            let size = V2i pEdges.Size
            let coords = Dither.getCoords size
            iter size (fun (c : V2i) ->
                let v = float pEdges.[c] / 255.0
                if v >= config.Threshold then
                    let mutable r = WeightedRegression2d.Empty

                    let inline regressionStable (r : WeightedRegression2d) =
                        if r.Count >= 15 then
                            r.GetQuality() > config.MinStability
                        else
                            false

                    let queue = System.Collections.Generic.List(100)
                    let queueHash = IntSet(1024)

                    let inline addHash (v : V2i) =
                        queueHash.Add (v.Y * size.X + v.X)

                    do
                        queue.HeapEnqueue(cmp, struct(v, c))
                        addHash c |> ignore

                        while r.Count < config.MaxGuessCount && queue.Count > 0 && not (regressionStable r) do
                            let struct(v, e) = queue.HeapDequeue(cmp)
                            r <- r.Add(V2d e + V2d.Half, v)

                            for no in neighbours do
                                let n = e + no
                                if n.AllGreaterOrEqual 0 && n.AllSmaller size then
                                    if addHash n then
                                        let v = float pEdges.[n] / 255.0
                                        if v >= config.GrowThreshold then
                                            queue.HeapEnqueue(cmp, struct(v, n))

                    if regressionStable r then
                        match r.TryGetPlaneInfo() with
                        | Some info when abs (info.Plane.Height(V2d c + V2d.Half)) <= config.Tolerance ->
                            let mutable info = info
                            r <- WeightedRegression2d.Empty

                            queue.Clear()
                            queueHash.Clear()
                            queue.HeapEnqueue(cmp, struct(v, c))
                            addHash c |> ignore

                            let all = System.Collections.Generic.HashSet()

                            while queue.Count > 0 do
                                let struct(vi,pi) = queue.HeapDequeue(cmp)
                                        
                                let err = abs (info.Plane.Height (V2d pi + V2d.Half))
                                if vi >= config.GrowThreshold && err <= config.Tolerance then
                                    all.Add pi |> ignore
                                    r <- r.Add(V2d pi + V2d.Half, vi)
                                    if regressionStable r then
                                        match r.TryGetPlaneInfo() with
                                        | Some i -> info <- i
                                        | None -> ()
                                    for no in neighbours do
                                        let n = pi + no
                                        if n.AllGreaterOrEqual 0 && n.AllSmaller size then
                                            if addHash n then 
                                                let vn = float pEdges.[n] / 255.0
                                                queue.HeapEnqueue(cmp, struct(vn, n))

                            if all.Count > 1 then
                                let ray = Ray2d(info.Center, info.XAxis)
                                let mutable tRange = Range1d.Invalid

                                for a in all do
                                    let t = ray.GetClosestPointTOn (V2d a + V2d.Half)
                                    tRange.ExtendBy t
                                    pEdges.[a] <- 0uy

                                if tRange.Size >= config.MinLength then
                                    let line = Line2d(ray.GetPointOnRay tRange.Min, ray.GetPointOnRay tRange.Max)

                                    lines.Add {
                                        regression = r
                                        info = info
                                        line = line
                                    }
                        | _ ->
                            ()
            )

        )
        //Log.stop()

        lines.ToArray()


    let private mergeLines (config : LineDetectionConfig) (ct : CancellationToken) (lines : DetectedLine[])  =
        let radius = config.MergeRadius

        let graph = LineGraph.Build(lines, config.MergeRadius)


        let tryMerge (config : LineDetectionConfig) (l : WeightedRegression2d) (ll : Line2d) (r : WeightedRegression2d) (rl : Line2d) =
            
            let h = l + r

            //let lq = l.GetQuality()
            //let rq = r.GetQuality()
            //let hq = h.GetQuality()

            //if hq > lq && hq > rq then
            //    let lInfo = l.TryGetPlaneInfo() |> Option.get
            //    let rInfo = r.TryGetPlaneInfo() |> Option.get
            //    let hInfo = h.TryGetPlaneInfo() |> Option.get
                
            //    let angle = Constant.DegreesPerRadian * acos (abs (Vec.dot lInfo.YAxis rInfo.YAxis) |> min 1.0)
            //    if angle < 5.0 then
                   
            //        let ray = Ray2d(hInfo.Center, hInfo.XAxis)
            //        let mutable tRange = Range1d.Invalid
            //        tRange.ExtendBy(ray.GetClosestPointTOn ll.P0)
            //        tRange.ExtendBy(ray.GetClosestPointTOn ll.P1)
            //        tRange.ExtendBy(ray.GetClosestPointTOn rl.P0)
            //        tRange.ExtendBy(ray.GetClosestPointTOn rl.P1)

            //        Some (h, hInfo, Line2d(ray.GetPointOnRay tRange.Min, ray.GetPointOnRay tRange.Max), hq)
            //    else
            //        None
            //else
            //    None

            ct.ThrowIfCancellationRequested()

            match l.TryGetPlaneInfo() with
            | Some lInfo ->
                match r.TryGetPlaneInfo() with
                | Some rInfo ->
                    match h.TryGetPlaneInfo() with
                    | Some hInfo ->
                        let lAngle = acos (abs (Vec.dot lInfo.XAxis hInfo.XAxis)) 
                        let rAngle = acos (abs (Vec.dot rInfo.XAxis hInfo.XAxis)) 
                        let lDist = max (abs (hInfo.Plane.Height ll.P0)) (abs (hInfo.Plane.Height ll.P1))
                        let rDist = max (abs (hInfo.Plane.Height rl.P0)) (abs (hInfo.Plane.Height rl.P1))

                        if lDist < config.Tolerance && rDist < config.Tolerance && lAngle < lInfo.AngularError && rAngle < rInfo.AngularError then
                            let ray = Ray2d(hInfo.Center, hInfo.XAxis)
                            let mutable tRange = Range1d.Invalid
                            tRange.ExtendBy(ray.GetClosestPointTOn ll.P0)
                            tRange.ExtendBy(ray.GetClosestPointTOn ll.P1)
                            tRange.ExtendBy(ray.GetClosestPointTOn rl.P0)
                            tRange.ExtendBy(ray.GetClosestPointTOn rl.P1)

                            Some (h, hInfo, Line2d(ray.GetPointOnRay tRange.Min, ray.GetPointOnRay tRange.Max))
                        else
                            None
                            

                    | None ->
                        None
                | None ->
                    None
            | _ ->
                None

        let output = System.Collections.Generic.List()

        let top = System.Collections.Generic.Queue()

        let rec traverse (config : LineDetectionConfig) (visited : System.Collections.Generic.HashSet<int>) (li : int) (lInfo : DetectedLine) =
            if visited.Add li then
                let conn = graph.GetNeighbours li
                if conn.Length = 0 then
                    output.Add lInfo
                else
                    let merged =
                        let d = lInfo.info.XAxis
                        conn 
                        |> Array.sortBy (fun c -> acos (abs (Vec.dot d lines.[c.Key].info.XAxis) |> min 1.0))
                        |> Array.choose (fun (KeyValue(ri, _k)) ->
                            if not (visited.Contains ri) then
                                let rInfo = lines.[ri]
                                match tryMerge config lInfo.regression lInfo.line rInfo.regression rInfo.line with
                                | Some(hReg, hInfo, hLine) -> 
                                    Some(ri, hReg, hInfo, hLine)
                                | None -> None
                            else
                                None
                        )

                    if merged.Length > 0 then   
                        let (ri, hReg, hInfo, hLine) = merged |> Array.maxBy (fun (_,_,i,_) -> i.Quality)
               
                        let (lint, rint) =
                            graph.TryGetConnection(li, ri) |> Option.get

                        let ropp = LinePart.opposite rint

                        for KeyValue(ii, ip) in graph.GetNeighbours(li, lint) do
                            if ii <> li && ii <> ri then
                                graph.Remove(li, ii) |> ignore
                                graph.Add(li, LinePart.Line, ii, ip)

                        for KeyValue(ii, ip) in graph.GetNeighbours(ri, rint) do
                            if ii <> li && ii <> ri then
                                graph.Add(li, LinePart.Line, ii, ip)

                        for KeyValue(ii, ip) in graph.GetNeighbours(ri, ropp) do
                            if ii <> li && ii <> ri then
                                graph.Add(li, lint, ii, ip)
                            

                        graph.Remove ri |> ignore
                        visited.Remove li |> ignore
                        visited.Add ri |> ignore
                        traverse config visited li { regression = hReg; info = hInfo; line = hLine }
                    else
                        output.Add lInfo
                        for KeyValue(ri, _) in conn do
                            top.Enqueue ri
                            //let rInfo = lines.[ri]
                            //traverse visited ri rInfo
                        
        let visited = System.Collections.Generic.HashSet()
        for i in 0 .. lines.Length - 1 do
            while top.Count > 0 do
                let ri = top.Dequeue()
                let rInfo = lines.[ri]
                traverse config visited ri rInfo
            let lInfo = lines.[i]
            traverse config visited i lInfo


        CSharpList.toArray output

    let detect (config : LineDetectionConfig) (ct : CancellationToken) (image : PixImage<byte>) =
       
        let rawLines = detectLinesPar config ct image
        
        //Log.startTimed "merge lines"
        let merged = mergeLines config ct rawLines
        //Log.stop()
     
        merged
        |> Array.filter (fun l -> l.info.Quality >= config.MinQuality)

    
    let detectAsync (config : LineDetectionConfig) (image : PixImage<byte>) =
        async {
            do! Async.SwitchToThreadPool()
            let! ct = Async.CancellationToken
            return detect config ct image
        }

