namespace LineDetector

open Aardvark.Base
open Aardvark.Geometry
open SixLabors.ImageSharp
open SixLabors.ImageSharp.Processing
open SixLabors.ImageSharp.Processing.Processors.Convolution

[<Struct>]
type PlaneFitInfo =
    {
        Trafo       : Trafo2d
        Plane       : Plane2d
        StdDev      : V2d
        Center      : V2d
        Count       : int
        Mass        : float
    }

    member x.Quality =
        let r = abs x.StdDev.X / (abs x.StdDev.X + abs x.StdDev.Y) 
        2.0 * r - 1.0

    member x.XAxis = x.Trafo.Forward.C0.XY
    member x.YAxis = x.Trafo.Forward.C1.XY

    member x.ToRegression2d() =
        let s = M22d.FromDiagonal(sqr x.StdDev * x.Mass)
        let u = x.Trafo.Forward.UpperLeftM22()
        let i = u * s * u.Transposed

        let sumSq = V3d(i.M00, i.M11, i.M01)
        Regression2d(sumSq, V2d.Zero, x.Center, x.Count, x.Mass)


    
and [<Struct>] Regression2d(sumSq : V3d, sum : V2d, ref : V2d, count : int, mass : float) =
        
    static member Empty = Regression2d(V3d.Zero, V2d.Zero, V2d.Zero, 0, 0.0)


            

    member x.Add(point : V2d, m : float) =
        if m <= 0.0 then
            x
        elif count = 0 then
            Regression2d(V3d.Zero, V2d.Zero, point, 1, m)
        else
            let p = point - ref
            let sumSq = sumSq + m * V3d(sqr p.X, sqr p.Y, p.X * p.Y)
            let sum = sum + m * p
            Regression2d(sumSq, sum, ref, count + 1, mass + m)

    member x.Remove(point : V2d, m : float) =
        if m <= 0.0 then
            x
        elif count = 1 then
            Regression2d.Empty
        else
            let p = point - ref
            let sumSq = sumSq - m * V3d(sqr p.X, sqr p.Y, p.X * p.Y)
            let sum = sum - m * p
            Regression2d(sumSq, sum, ref, count - 1, mass - m)

    member x.IsEmpty =
        count = 0

    member x.Center = 
        if count = 1 then ref
        else ref + sum / mass

    member x.Mass = mass

    member x.Count = count

    member x.InertiaTensor =
        if mass <= 0.0 || count < 2 then
            M22d.Zero
        else
            // https://en.wikipedia.org/wiki/Parallel_axis_theorem
            let a = sum / mass
            let ixx = sumSq.Y - mass*sqr a.Y
            let iyy = sumSq.X - mass*sqr a.X
            let ixy = mass*a.X*a.Y - sumSq.Z
            M22d(ixx, ixy, ixy, iyy)
                
    member x.GetQuality() =
        if count >= 2 then
            match SVD.Decompose x.InertiaTensor with
            | Some (_u, s, _vt) -> 
                let stddev = abs s.Diagonal / mass |> sqrt
                let r = abs stddev.X / (abs stddev.X + abs stddev.Y) 
                2.0 * r - 1.0
            | None -> 0.0
        else
            0.0
            

    member x.GetStdDevRatio() =
        if count >= 2 then
            match SVD.Decompose x.InertiaTensor with
            | Some (_u, s, _vt) -> 
                let stddev = abs s.Diagonal / mass |> sqrt
                abs stddev.X / abs stddev.Y
            | None -> 0.0
        else
            0.0
    member x.TryGetPlaneInfo() =
        if count >= 2 then
            match SVD.Decompose x.InertiaTensor with
            | Some (u, s, _vt) ->
                let avg = x.Center

                let x = u.C1
                let y = u.C0
                Some {
                    Plane = Plane2d(y, avg)
                    Trafo = Trafo2d.FromBasis(x, y, avg)
                    StdDev = abs s.Diagonal / mass |> sqrt
                    Center = avg
                    Count = count
                    Mass = mass
                }
            | None ->
                None
        else
            None
                
    member private x.Ref = ref
    member private x.Sum = sum
    member private x.SumSq = sumSq

    static member (+) (l : Regression2d, r : Regression2d) =
        if l.IsEmpty then r
        elif r.IsEmpty then l
        elif l.Count = 1 then r.Add(l.Center, l.Mass)
        elif r.Count = 1 then l.Add(r.Center, r.Mass)
        elif l.Count > r.Count then 
            let ref = l.Ref
            let o = r.Ref - ref
            let cnt = l.Count + r.Count
            let mass = l.Mass + r.Mass
            let sum = l.Sum + r.Sum + r.Mass * o

            let sumSq =
                let m = r.Sum
                l.SumSq + V3d(
                    r.SumSq.X   + 2.0*o.X*m.X       + r.Mass*sqr o.X,
                    r.SumSq.Y   + 2.0*o.Y*m.Y       + r.Mass*sqr o.Y,
                    r.SumSq.Z   + o.X*m.Y + o.Y*m.X + r.Mass*o.X*o.Y
                )

            Regression2d(sumSq, sum, ref, cnt, mass)
        else
            let ref = r.Ref
            let o = l.Ref - ref
            let cnt = l.Count + r.Count
            let mass = l.Mass + r.Mass
            let sum = r.Sum + l.Sum + l.Mass * o

            let sumSq =
                let m = l.Sum
                r.SumSq + V3d(
                    l.SumSq.X   + 2.0*o.X*m.X       + l.Mass*sqr o.X,
                    l.SumSq.Y   + 2.0*o.Y*m.Y       + l.Mass*sqr o.Y,
                    l.SumSq.Z   + o.X*m.Y + o.Y*m.X + l.Mass*o.X*o.Y
                )

            Regression2d(sumSq, sum, ref, cnt, mass)

[<AutoOpen>]
module private LineDetectionHelpers =
    let neighbours =
        [|
            for x in -2 .. 2 do
                for y in -2 .. 2 do
                    if x <> 0 || y <> 0 then yield V2i(x,y)
        |]

    let cmp =
        System.Func<_,_,_>(fun struct(a,_) struct(b,_) -> compare b a)
        
    let cmp3 =
        System.Func<_,_,_>(fun struct(a,_,_) struct(b,_,_) -> compare b a)

type DetectedLine =
    {
        regression  : Regression2d
        line        : Line2d
        info        : PlaneFitInfo
    }
    
type LinePart =
    | P0        = 0x0001
    | P1        = 0x0002
    | Line      = 0x0004
    | All       = 0x0007

module LinePart =
    let opposite (p : LinePart) =
        match p with
        | LinePart.P0 -> LinePart.P1
        | LinePart.P1 -> LinePart.P0
        | _ -> p

type private LineConnectionInfo() =
    let p0 = IntDict<LinePart>()
    let p1 = IntDict<LinePart>()
    let line = IntDict<LinePart>()

    member x.Count =
        p0.Count + p1.Count + line.Count

    member x.GetCount(flags : LinePart) =
        let mutable c = 0
        if flags.HasFlag LinePart.P0 then c <- c + p0.Count
        if flags.HasFlag LinePart.P1 then c <- c + p1.Count
        if flags.HasFlag LinePart.Line then c <- c + line.Count
        c

    member x.IsEmpty =
        p0.Count = 0 && p1.Count = 0 && line.Count = 0

    member x.P0 = p0
    member x.P1 = p1
    member x.Line = line

    member x.Add(lp : LinePart, ri : int, rp : LinePart) =
        match lp with
        | LinePart.P0 -> p0.[ri] <- rp
        | LinePart.P1 -> p1.[ri] <- rp
        | _ -> line.[ri] <- rp

    member x.Remove(ri : int) =
        if not (p0.Remove ri) then
            if not (p1.Remove ri) then
                line.Remove ri
            else
                true
        else
            true

    member x.TryRemove(ri : int) =
        match p0.TryRemove ri with
        | (true, rp) -> Some (LinePart.P0, rp)
        | _ ->
            match p1.TryRemove ri with
            | (true, rp) -> Some (LinePart.P1, rp)
            | _ ->
                match line.TryRemove ri with
                | (true, rp) -> Some (LinePart.Line, rp)
                | _ -> None

type LineGraph() =
    let neighbours = IntDict<LineConnectionInfo>()
    let nodes = IntSet()

    member x.Add(li : int, lp : LinePart, ri : int, rp : LinePart) =
        if li <> ri then
            let ln = neighbours.GetOrCreate(li, fun _ -> LineConnectionInfo())
            ln.Add(lp, ri, rp)
            nodes.Add li |> ignore

            let rn = neighbours.GetOrCreate(ri, fun _ -> LineConnectionInfo())
            rn.Add(rp, li, lp)
            nodes.Add ri |> ignore

    member x.Remove(li : int, ri : int) =
        match neighbours.TryGetValue li with
        | (true, ln) -> 
            if ln.Remove ri && ln.IsEmpty then neighbours.Remove li |> ignore
        | _ ->
            ()
                
        match neighbours.TryGetValue ri with
        | (true, rn) -> 
            if rn.Remove li then
                if rn.IsEmpty then neighbours.Remove ri |> ignore
                true
            else
                false
        | _ ->
            false

    member x.TryGetConnection(li : int, ri : int) =
        match neighbours.TryGetValue li with
        | (true, ns) ->
            match ns.P0.TryGetValue ri with
            | (true, rs) -> Some(LinePart.P0, rs)
            | _ ->
                match ns.P1.TryGetValue ri with
                | (true, rs) -> Some(LinePart.P1, rs)
                | _ ->
                    match ns.Line.TryGetValue ri with
                    | (true, rs) -> Some(LinePart.Line, rs)
                    | _ -> None
        | _ ->
            None

    member x.Remove(li : int, lp : LinePart, ri : int) =
        match neighbours.TryGetValue li with
        | (true, ln) -> 
            let result = 
                match lp with
                | LinePart.P0 -> ln.P0.Remove ri
                | LinePart.P1 -> ln.P1.Remove ri
                | _ -> ln.Line.Remove ri
            if result && ln.IsEmpty then neighbours.Remove li |> ignore
            result
        | _ ->
            false
            
    member x.GetNeighbours(li : int, part : LinePart) =
        match neighbours.TryGetValue li with
        | (true, ln) ->
            let result = Array.zeroCreate (ln.GetCount part)
            let mutable offset = 0
                
            if part.HasFlag LinePart.P0 then
                ln.P0.CopyTo(result, offset)
                offset <- offset + ln.P0.Count

            if part.HasFlag LinePart.P1 then
                ln.P1.CopyTo(result, offset)
                offset <- offset + ln.P1.Count
                    
            if part.HasFlag LinePart.Line then
                ln.Line.CopyTo(result, offset)
                offset <- offset + ln.Line.Count

            result
        | _ ->
            [||]
                
    member x.GetNeighbours(li : int) =
        x.GetNeighbours(li, LinePart.All)

    member x.Remove(node : int) =
        if nodes.Remove node then
            match neighbours.TryRemove node with
            | (true, ns) ->
                for KeyValue(ri, _rp) in ns.P0 do
                    match neighbours.TryGetValue ri with
                    | (true, rn) -> 
                        if rn.Remove node && rn.IsEmpty then neighbours.Remove ri |> ignore
                    | _ -> ()
                        
                for KeyValue(ri, _rp) in ns.P1 do
                    match neighbours.TryGetValue ri with
                    | (true, rn) -> 
                        if rn.Remove node && rn.IsEmpty then neighbours.Remove ri |> ignore
                    | _ -> ()

                for KeyValue(ri, _rp) in ns.Line do
                    match neighbours.TryGetValue ri with
                    | (true, rn) -> 
                        if rn.Remove node && rn.IsEmpty then neighbours.Remove ri |> ignore
                    | _ -> ()
            | _ ->
                ()
            true
        else
            false
        
    member x.Nodes =
        nodes

    static member Build(lines : DetectedLine[], radius : float) =
            
        let linePoints = Array.zeroCreate (2 * lines.Length)
        let mutable oi = 0
        for i in 0 .. lines.Length - 1 do
            let { line = l } = lines.[i]
            linePoints.[oi] <- l.P0
            linePoints.[oi+1] <- l.P1
            oi <- oi + 2

        let planes = lines |> Array.map (fun { line = l } -> l.Plane2d.Normalized)
        let tree = linePoints.CreateRkdTree(Metric.Euclidean, 1E-2)

        let graph = LineGraph()

        for li in 0 .. lines.Length - 1 do
            let { line = l } = lines.[li]
            let lp = planes.[li]

            let q = tree.CreateClosestToLineQuery(radius, 0)
            let closest = tree.GetClosest(q, l.P0, l.P1)

            for id in closest do
                let ri = int id.Index / 2
                if ri > li then
                    let { line = r } = lines.[ri]
                    let rp = planes.[ri]

                    let d00 = Vec.distance l.P0 r.P0
                    let d01 = Vec.distance l.P0 r.P1
                    let d10 = Vec.distance l.P1 r.P0
                    let d11 = Vec.distance l.P1 r.P1

                    let struct(distance, ls, rs) =
                        Array.minBy (fun struct(a,_,_) -> a) [|
                            struct (d00, LinePart.P0, LinePart.P0) 
                            struct (d01, LinePart.P0, LinePart.P1) 
                            struct (d10, LinePart.P1, LinePart.P0) 
                            struct (d11, LinePart.P1, LinePart.P1) 
                        |]

                    if distance <= radius then
                        graph.Add(li, ls, ri, rs)
                    else
                        let dl0 = abs (lp.Height r.P0)
                        let dl1 = abs (lp.Height r.P1)
                        let d0l = abs (rp.Height l.P0)
                        let d1l = abs (rp.Height l.P1)
                        let struct(distance, ls, rs) =
                            Array.minBy (fun struct(a,_,_) -> a) [|
                                struct (dl0, LinePart.Line, LinePart.P0) 
                                struct (dl1, LinePart.Line, LinePart.P1) 
                                struct (d0l, LinePart.P0, LinePart.Line) 
                                struct (d1l, LinePart.P1, LinePart.Line) 
                            |]
                        if distance <= radius then
                            graph.Add(li, ls, ri, rs)
                        elif l.Intersects r then
                            graph.Add(li, LinePart.Line, ri, LinePart.Line)
                
        graph


module LineDetector =
    let detectLines (image : PixImage<byte>) =  
        
        let threshold = 0.4
        let growThreshold = 0.1
        let tolerance = 3.5
        let minLength = 4.0
        let minStability = 0.8
        let maxGuessCount = 40

        Log.startTimed "detecting edges"
        use img = image.ToImage()
        img.Mutate (fun ctx ->
            
            ctx.Grayscale()
                //.HistogramEqualization()
                //.GaussianBlur(2.0f)
                .DetectEdges(EdgeDetectorKernel.Laplacian3x3, true)
            |> ignore
        )

        img.SaveAsJpeg @"C:\Users\Schorsch\Desktop\bla.png"
        let edges = img.ToPixImage().ToPixImage<byte>()
        Log.stop()

        let size = edges.Size

        let lines = System.Collections.Generic.List()


        Log.startTimed "finding lines"
        NativeMatrix.using (edges.GetChannel 0L) (fun pEdges ->
            
            pEdges |> NativeMatrix.iter (fun c v ->
                let c = V2i c
                let v = float v / 255.0
                if v >= threshold then
                    let mutable r = Regression2d.Empty

                    let regressionStable (r : Regression2d) =
                        if r.Count >= 4 then
                            match r.TryGetPlaneInfo() with
                            | Some i ->
                                i.Quality > minStability
                                //if Fun.IsTiny i.StdDev.Y then
                                //    i.StdDev.X >= minStability
                                //else
                                //    let r = i.StdDev.X / i.StdDev.Y
                                //    r >= minStability
                            | None ->
                                false
                        else
                            false

                    do
                        let queue = System.Collections.Generic.List()
                        let queueHash = System.Collections.Generic.HashSet()
                        queue.HeapEnqueue(cmp, struct(v, c))
                        queueHash.Add c |> ignore

                        while r.Count < maxGuessCount && queue.Count > 0 && not (regressionStable r) do
                            let struct(v, e) = queue.HeapDequeue(cmp)
                            r <- r.Add(V2d e + V2d.Half, v)

                            for no in neighbours do
                                let n = e + no
                                if n.AllGreaterOrEqual 0 && n.AllSmaller size then
                                    if queueHash.Add n then
                                        let v = float pEdges.[n] / 255.0
                                        if v >= growThreshold then
                                            queue.HeapEnqueue(cmp, struct(v, n))

                    if regressionStable r then
                        match r.TryGetPlaneInfo() with
                        | Some info when abs (info.Plane.Height(V2d c + V2d.Half)) <= tolerance ->
                            let mutable info = info
                            r <- Regression2d.Empty
                            let queue = System.Collections.Generic.Queue()
                            let queueHash = System.Collections.Generic.HashSet()
                            queue.Enqueue(c)
                            queueHash.Add c |> ignore

                            let all = System.Collections.Generic.HashSet()

                            while queue.Count > 0 do
                                let pi = queue.Dequeue()
                                let vi = float pEdges.[pi] / 255.0
                                let err = abs (info.Plane.Height (V2d pi))
                                if vi >= growThreshold && err <= tolerance then
                                    all.Add pi |> ignore
                                    r <- r.Add(V2d pi + V2d.Half, vi)
                                    if regressionStable r then
                                        match r.TryGetPlaneInfo() with
                                        | Some i -> info <- i
                                        | None -> ()
                                    for no in neighbours do
                                        let n = pi + no
                                        if n.AllGreaterOrEqual 0 && n.AllSmaller size then
                                            if queueHash.Add n then queue.Enqueue n

                            if all.Count > 1 then
                                let n = 0.5 * info.Plane.Normal + V2d.Half


                                let ray = Ray2d(info.Center, info.XAxis)
                                let mutable tRange = Range1d.Invalid

                                for a in all do 
                                    let t = ray.GetClosestPointTOn (V2d a)
                                    tRange.ExtendBy t
                                    pEdges.[a] <- 0uy
                                if tRange.Size >= minLength then
                                    lines.Add {
                                        regression = r
                                        info = info
                                        line = Line2d(ray.GetPointOnRay tRange.Min, ray.GetPointOnRay tRange.Max)
                                    }
                        | _ ->
                            ()

            )

        )
        Log.stop()

        CSharpList.toArray lines

    let mergeLines (lines : DetectedLine[])  =
    
        let radius = 20.0

        let linePoints = Array.zeroCreate (2 * lines.Length)
        let mutable oi = 0
        for i in 0 .. lines.Length - 1 do
            let { line = l } = lines.[i]
            linePoints.[oi] <- l.P0
            linePoints.[oi+1] <- l.P1
            oi <- oi + 2

        let planes = lines |> Array.map (fun { line = l } -> l.Plane2d.Normalized)
        let tree = linePoints.CreateRkdTree(Metric.Euclidean, 1E-2)

        let graph = LineGraph()

        for li in 0 .. lines.Length - 1 do
            let { line = l } = lines.[li]
            let lp = planes.[li]

            let q = tree.CreateClosestToLineQuery(radius, 0)
            let closest = tree.GetClosest(q, l.P0, l.P1)

            for id in closest do
                let ri = int id.Index / 2
                if ri > li then
                    let { line = r } = lines.[ri]
                    let rp = planes.[ri]

                    let d00 = Vec.distance l.P0 r.P0
                    let d01 = Vec.distance l.P0 r.P1
                    let d10 = Vec.distance l.P1 r.P0
                    let d11 = Vec.distance l.P1 r.P1

                    let struct(distance, ls, rs) =
                        Array.minBy (fun struct(a,_,_) -> a) [|
                            struct (d00, LinePart.P0, LinePart.P0) 
                            struct (d01, LinePart.P0, LinePart.P1) 
                            struct (d10, LinePart.P1, LinePart.P0) 
                            struct (d11, LinePart.P1, LinePart.P1) 
                        |]

                    if distance <= radius then
                        graph.Add(li, ls, ri, rs)
                    else
                        let dl0 = abs (lp.Height r.P0)
                        let dl1 = abs (lp.Height r.P1)
                        let d0l = abs (rp.Height l.P0)
                        let d1l = abs (rp.Height l.P1)
                        let struct(distance, ls, rs) =
                            Array.minBy (fun struct(a,_,_) -> a) [|
                                struct (dl0, LinePart.Line, LinePart.P0) 
                                struct (dl1, LinePart.Line, LinePart.P1) 
                                struct (d0l, LinePart.P0, LinePart.Line) 
                                struct (d1l, LinePart.P1, LinePart.Line) 
                            |]
                        if distance <= radius then
                            graph.Add(li, ls, ri, rs)
                        elif l.Intersects r then
                            graph.Add(li, LinePart.Line, ri, LinePart.Line)
                        

        let tryMerge (l : Regression2d) (ll : Line2d) (r : Regression2d) (rl : Line2d) =
            let h = l + r
            match l.TryGetPlaneInfo() with
            | Some lInfo ->
                let lScore = lInfo.Quality //lInfo.StdDev.X / lInfo.StdDev.Y
                match r.TryGetPlaneInfo() with
                | Some rInfo ->
                    let rScore = rInfo.Quality // rInfo.StdDev.X / rInfo.StdDev.Y
                    match h.TryGetPlaneInfo() with
                    | Some hInfo ->
                        let hScore = hInfo.Quality // hInfo.StdDev.X / hInfo.StdDev.Y
                        if hScore >= lScore && hScore >= rScore then
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
            | None ->
                None

        let output = System.Collections.Generic.List()

        let rec traverse (visited : System.Collections.Generic.HashSet<int>) (li : int) (lInfo : DetectedLine) =
            if visited.Add li then
                let conn = graph.GetNeighbours li
                if conn.Length = 0 then
                    output.Add lInfo
                else
                    let merged =
                        conn |> Array.tryPick (fun (KeyValue(ri, _k)) ->
                            if not (visited.Contains ri) then
                                let rInfo = lines.[ri]
                                match tryMerge lInfo.regression lInfo.line rInfo.regression rInfo.line with
                                | Some(hReg, hInfo, hLine) -> 
                                    Some(ri, hReg, hInfo, hLine)
                                | None -> None
                            else
                                None
                        )

                    match merged with
                    | Some (ri, hReg, hInfo, hLine) ->
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
                        traverse visited li { regression = hReg; info = hInfo; line = hLine }
                    | None ->
                        output.Add lInfo
                        for KeyValue(ri, _) in conn do
                            let rInfo = lines.[ri]
                            traverse visited ri rInfo
                        
        let visited = System.Collections.Generic.HashSet()
        for i in 0 .. lines.Length - 1 do
            let lInfo = lines.[i]
            traverse visited i lInfo


        CSharpList.toArray output



