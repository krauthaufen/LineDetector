namespace LineDetector

open Aardvark.Base
open Aardvark.Geometry
open SixLabors.ImageSharp
open SixLabors.ImageSharp.Processing
open SixLabors.ImageSharp.Processing.Processors.Convolution
open Microsoft.FSharp.NativeInterop

#nowarn "9"

[<Struct>]
type PlaneFitInfo =
    {
        Trafo           : Trafo2d
        Plane           : Plane2d
        StdDev          : V2d
        Center          : V2d
        Count           : int
        Mass            : float
        AngularError    : float
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
        
    static let getEigenvalues (sumSq : V3d) (sum : V2d) (mass : float) = 
        let a = sum / mass
        let ixx = sumSq.Y - mass*sqr a.Y
        let iyy = sumSq.X - mass*sqr a.X
        let ixy = mass*a.X*a.Y - sumSq.Z
        let mutable struct(a, b) = Polynomial.RealRootsOfNormed(-(iyy+ixx), ixx*iyy - sqr ixy)
        if a < b then Fun.Swap(&a, &b)
        V2d(a,b)

    static let decomposeInertia (sumSq : V3d) (sum : V2d) (mass : float) =
        // https://en.wikipedia.org/wiki/Parallel_axis_theorem
        let a = sum / mass
        let ixx = sumSq.Y - mass*sqr a.Y
        let iyy = sumSq.X - mass*sqr a.X
        let ixy = mass*a.X*a.Y - sumSq.Z

        //// (m.M00 - a) * (m.M11 - a) - m.M01*m.M10 = 0
        //// m.M00*m.M11 - a*(m.M11 + m.M00) + a^2 - m.M01*m.M10
        //// a^2 + a*-(m.M11 + m.M00) + (m.M00*m.M11 - m.M01*m.M10)
        let mutable struct(a, b) = Polynomial.RealRootsOfNormed(-(iyy+ixx), ixx*iyy - sqr ixy)

        if a < b then Fun.Swap(&a, &b)

        if System.Double.IsNaN a || System.Double.IsNaN b then
            struct(M22d.Identity, V2d.Zero)
        else
            let a0 = V2d(ixx - a, ixy).Rot90
            let a1 = V2d(ixy, iyy - a).Rot90
            //let b0 = V2d(ixx - b, ixy).Rot90
            //let b1 = V2d(ixy, iyy - b).Rot90

            let la0 = Vec.lengthSquared a0
            let la1 = Vec.lengthSquared a1
            //let lb0 = Vec.lengthSquared b0
            //let lb1 = Vec.lengthSquared b1

            //let mutable x = V2d.Zero
            //let mutable y = V2d.Zero

            let x = 
                if la0 > la1 then a0 / sqrt la0
                else a1 / sqrt la1
            let y = x.Rot90


            //if la0 > la1 then
            //    if la0 > lb0 then
            //        if la0 > lb1 then  
            //            x <- a0 / sqrt la0
            //            y <- x.Rot90
            //        else 
            //            y <- b1 / sqrt lb1
            //            x <- y.Rot270
            //    else
            //        if lb0 > lb1 then 
            //            y <- b0 / sqrt lb0
            //            x <- y.Rot270
            //        else 
            //            y <- b1 / sqrt lb1
            //            x <- y.Rot270
            //else
            //    if la1 > lb0 then
            //        if la1 > lb1 then 
            //            x <- a1 / sqrt la1
            //            y <- x.Rot90
            //        else 
            //            y <- b1 / sqrt lb1
            //            x <- y.Rot270
            //    else
            //        if lb0 > lb1 then 
            //            y <- b0 / sqrt lb0
            //            x <- y.Rot270
            //        else 
            //            y <- b1 / sqrt lb1
            //            x <- y.Rot270

            struct(M22d.FromCols(x,y), V2d(a, b))



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
            let s = getEigenvalues sumSq sum mass
            let stddev = abs s / mass |> sqrt
            let r = abs stddev.X / (abs stddev.X + abs stddev.Y) 
            2.0 * r - 1.0
        else
            0.0
            

    member x.GetStdDevRatio() =
        if count >= 2 then
            let s = getEigenvalues sumSq sum mass
            let stddev = abs s / mass |> sqrt
            abs stddev.X / abs stddev.Y
        else
            0.0


    member x.TryGetPlaneInfo() =
        if count >= 2 then
            let struct(u, s) = decomposeInertia sumSq sum mass
            if s.X <= 0.0 then
                None
            else
                let avg = x.Center
                
                let ev = abs s
                let measurementNoise = ev.Y / (x.Mass * float x.Count)
                let noiseCov = 4.0 * ev * measurementNoise
                let inline fppf (x : float) (d1 : int) (d2 : int) =
                    FDistr.invCDF (float d1) (float d2) x

                let inline fisherStatistic (n : int) (confidence : float) (dof : int) : float =
                    fppf confidence dof (n - dof)
                
                let inline applyErrorScaling (nominal : V2d) (err : V2d) (n : int) : V2d =
                    nominal * V2d.PN - err |> abs


                let z = fisherStatistic x.Count 0.95 x.Count
                let err = z * sqrt noiseCov
                
                let hyp = applyErrorScaling ev err x.Count
                let n = hyp |> sqrt
                let angle = atan2 n.Y n.X

                let x = u.C1
                let y = u.C0


                Some {
                    Plane = Plane2d(y, avg)
                    Trafo = Trafo2d.FromBasis(x, y, avg)
                    StdDev = abs s / mass |> sqrt
                    Center = avg
                    Count = count
                    Mass = mass
                    AngularError = angle
                }
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
            for x in -1 .. 1 do
                for y in -1 .. 1 do
                    if x <> 0 || y <> 0 then yield V2i(x,y)
        |]

    let cmp =
        System.Func<_,_,_>(fun struct(a,_) struct(b,_) -> compare b a)
        
    let cmp3 =
        System.Func<_,_,_>(fun struct(a,_,_) struct(b,_,_) -> compare b a)

type DetectedLine =
    {
        regression  : Regression2d
        stddev      : float
        avgResponse : float
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

module Dither =
    let matrix =
        let ass = typeof<LineDetectionConfig>.Assembly

        let resName = 
            ass.GetManifestResourceNames()
            |> Array.find (fun n -> n.EndsWith "forced-1024.bin")

        use s = ass.GetManifestResourceStream resName
        use r = new System.IO.BinaryReader(s)
        let res = Array.zeroCreate (int (s.Length / 4L))

        for i in 0 .. res.Length - 1 do
            res.[i] <- r.ReadInt32()

        Matrix<int>(res, 0L, V2l(1024, 1024), V2l(1L, 1024L))

    let coords =
        let l = Array.zeroCreate (int matrix.SX * int matrix.SY)
        matrix.ForeachXYIndex(fun x y i ->
            let v = matrix.[i]
            l.[v] <- V2i(int x, int y)
        )
        l

    let getCoords (size : V2i) =
        let blocks = ceil (V2d size / V2d matrix.Size) |> V2i

        seq {
            for by in 0 .. blocks.Y-1 do
                for bx in 0 .. blocks.X-1 do
                    let o = V2i(bx, by) * V2i matrix.Size

                    if (o + V2i matrix.Size).AllSmaller size then
                        if o = V2i.Zero then yield! coords
                        else yield! coords |> Seq.map (fun c -> c + o)
                    else
                        for c in coords do
                            let cc = o + c
                            if cc.AllSmaller size then
                                yield cc
                        
        }


module LineDetector =

    let private par (blocks : V2i) (mat : NativeMatrix<'a>) (action : V2i -> NativeMatrix<'a> -> unit) =
        if blocks.AllSmallerOrEqual 1 then
            action V2i.Zero mat
        else
            let s = V2d mat.Size / V2d blocks |> ceil |> V2i
            let matrices =
                [|
                    for y in 0 .. blocks.Y - 1 do
                        for x in 0 .. blocks.X - 1 do
                            let offset = s * V2i(x,y)
                            let size = min s (V2i mat.Size - offset)
                            offset, mat.SubMatrix(offset, size)
                |]

            System.Threading.Tasks.Parallel.ForEach(matrices, fun (offset, mat) -> action offset mat) |> ignore

    module MatrixInfo =
        let iter2 (threads : int) (action : int64 -> int64 -> unit) (a : MatrixInfo) (b : MatrixInfo) =
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

                let o = System.Threading.Tasks.ParallelOptions(MaxDegreeOfParallelism = threads)
                System.Threading.Tasks.Parallel.ForEach(blocks, o, action) |> ignore
            else
                a.ForeachIndex(b, fun ai bi ->
                    action.Invoke(ai, bi)
                )

    
    let private edgeDetectLaplace (threads : int) (image : PixImage<byte>) =
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

                (edgeInfo, srcInfo) ||> MatrixInfo.iter2 threads (fun (ei : int64) (si : int64) ->
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
       

    open System.Threading

    let private detectLines (config : LineDetectionConfig) (image : PixImage<byte>) =  
        
        //let configThreshold = 0.4
        //let configGrowThreshold = 0.2
        //let configTolerance = 3.5
        //let configMinLength = 4.0
        //let configMinStability = 0.8
        //let configMaxGuessCount = 20

        Log.startTimed "detecting edges"
        let edges = edgeDetectLaplace 1 image
        Log.stop()
        
        let globalLines = System.Collections.Generic.List()

        //let regions = PixImage<byte>(Col.Format.RGBA, image.Size)
        //let mutable rm = regions.GetMatrix<C4b>()
        //rm.SetMap(edges.GetChannel 0L, fun v -> C4b(v,v,v,255uy)) |> ignore


        Log.startTimed "finding lines"
        NativeMatrix.using (edges.GetChannel 0L) (fun pEdges ->

            let blocks = V2i.II
                //let threads = 8 * System.Environment.ProcessorCount
                //let mutable x = sqrt (float threads) |> floor |> int
                //while threads % x <> 0 do
                //    x <- x - 1
                //let y = threads / x
                //V2i(y, x)

            par blocks pEdges (fun offset pEdges ->
                let lines = System.Collections.Generic.List()
                
                let size = V2i pEdges.Size

                pEdges |> NativeMatrix.iter (fun c v ->
                    let c = V2i c
                    let v = float v / 255.0
                    if v >= config.Threshold then
                        let mutable r = Regression2d.Empty

                        let regressionStable (r : Regression2d) =
                            if r.Count >= 15 then
                                r.GetQuality() > config.MinStability
                            else
                                false

                        do
                            let queue = System.Collections.Generic.List()
                            let queueHash = System.Collections.Generic.HashSet()


                            queue.HeapEnqueue(cmp, struct(v, c))
                            queueHash.Add c |> ignore

                            while r.Count < config.MaxGuessCount && queue.Count > 0 && not (regressionStable r) do
                                let struct(v, e) = queue.HeapDequeue(cmp)
                                r <- r.Add(V2d e + V2d offset + V2d.Half, v)

                                for no in neighbours do
                                    let n = e + no
                                    if n.AllGreaterOrEqual 0 && n.AllSmaller size then
                                        if queueHash.Add n then
                                            let v = float pEdges.[n] / 255.0
                                            if v >= config.GrowThreshold then
                                                queue.HeapEnqueue(cmp, struct(v, n))

                        if regressionStable r then
                            match r.TryGetPlaneInfo() with
                            | Some info when abs (info.Plane.Height(V2d c + V2d offset + V2d.Half)) <= config.Tolerance ->
                                let mutable info = info
                                r <- Regression2d.Empty
                                let queue = System.Collections.Generic.List()
                                let queueHash = System.Collections.Generic.HashSet()
                                queue.HeapEnqueue(cmp, struct(1.0, c))
                                queueHash.Add c |> ignore

                                let all = System.Collections.Generic.HashSet()

                                let mutable sum = 0.0
                                let mutable sumSq = 0.0

                                while queue.Count > 0 do
                                    let struct(vi,pi) = queue.HeapDequeue(cmp)
                                    
                                    let err = abs (info.Plane.Height (V2d pi + V2d offset + V2d.Half))
                                    if vi >= config.GrowThreshold && err <= config.Tolerance then
                                        sum <- sum + vi
                                        sumSq <- sumSq + sqr vi
                                        all.Add pi |> ignore
                                        r <- r.Add(V2d pi + V2d offset + V2d.Half, vi)
                                        if regressionStable r then
                                            match r.TryGetPlaneInfo() with
                                            | Some i -> info <- i
                                            | None -> ()
                                        for no in neighbours do
                                            let n = pi + no
                                            if n.AllGreaterOrEqual 0 && n.AllSmaller size then
                                                if queueHash.Add n then 
                                                    let vn = float pEdges.[n] / 255.0
                                                    queue.HeapEnqueue(cmp, struct(vn, n))

                                if all.Count > 1 then
                                    let ray = Ray2d(info.Center, info.XAxis)
                                    let mutable tRange = Range1d.Invalid

                                    for a in all do 
                                        let t = ray.GetClosestPointTOn (V2d a + V2d offset + V2d.Half)
                                        tRange.ExtendBy t
                                        let v = float pEdges.[a] / 255.0
                                        pEdges.[a] <- 0uy
                                        //rm.[a] <- (c.ToC3f() * float32 v).ToC4b()

                                    if tRange.Size >= config.MinLength then
                                        let avg = sum / float all.Count
                                        let avgSq = sumSq / float all.Count

                                        let var = (sumSq - float all.Count * sqr avg) / float (all.Count - 1)

                                        let line = Line2d(ray.GetPointOnRay tRange.Min, ray.GetPointOnRay tRange.Max)

                                        //rm.SetLine(line.P0, line.P1, C4b(255uy - c.R, 255uy - c.G, 255uy - c.B, 255uy))

                                        lines.Add {
                                            regression = r
                                            info = info
                                            avgResponse = avg
                                            stddev = sqrt var
                                            line = Line2d(ray.GetPointOnRay tRange.Min, ray.GetPointOnRay tRange.Max)
                                        }
                            | _ ->
                                ()

                )

                lock globalLines (fun () ->
                    globalLines.AddRange lines
                )

            )

        )
        Log.stop()

        //regions.SaveAsImage (sprintf @"C:\Users\Schorsch\Desktop\debug\%A.jpg" (System.Guid.NewGuid()))


        CSharpList.toArray globalLines
   
    let private detectLinesPar (config : LineDetectionConfig) (image : PixImage<byte>) =  
        let threads = if config.Threads <= 0 then System.Environment.ProcessorCount else config.Threads
        
        if threads = 1 then
            detectLines config image
        else
            Log.startTimed "detecting edges"
            let edges = edgeDetectLaplace threads image
            Log.stop()

            let lines = System.Collections.Concurrent.ConcurrentBag()
            Log.startTimed "finding lines"
            NativeMatrix.using (edges.GetChannel 0L) (fun pEdges ->
                let size = V2i pEdges.Size
                let coords = Dither.getCoords size
                let o = System.Threading.Tasks.ParallelOptions(MaxDegreeOfParallelism = threads)
                System.Threading.Tasks.Parallel.ForEach(coords, o, fun (c : V2i) ->
                    let v = float pEdges.[c] / 255.0
                    if v >= config.Threshold then
                        let mutable r = Regression2d.Empty

                        let inline regressionStable (r : Regression2d) =
                            if r.Count >= 15 then
                                r.GetQuality() > config.MinStability
                            else
                                false

                        let queue = System.Collections.Generic.List()
                        let queueHash = System.Collections.Generic.HashSet()


                        do
                            queue.HeapEnqueue(cmp, struct(v, c))
                            queueHash.Add c |> ignore

                            while r.Count < config.MaxGuessCount && queue.Count > 0 && not (regressionStable r) do
                                let struct(v, e) = queue.HeapDequeue(cmp)
                                r <- r.Add(V2d e + V2d.Half, v)

                                for no in neighbours do
                                    let n = e + no
                                    if n.AllGreaterOrEqual 0 && n.AllSmaller size then
                                        if queueHash.Add n then
                                            let v = float pEdges.[n] / 255.0
                                            if v >= config.GrowThreshold then
                                                queue.HeapEnqueue(cmp, struct(v, n))

                        if regressionStable r then
                            match r.TryGetPlaneInfo() with
                            | Some info when abs (info.Plane.Height(V2d c + V2d.Half)) <= config.Tolerance ->
                                let mutable info = info
                                r <- Regression2d.Empty

                                queue.Clear()
                                queueHash.Clear()
                                queue.HeapEnqueue(cmp, struct(v, c))
                                queueHash.Add c |> ignore

                                let all = System.Collections.Generic.HashSet()

                                let mutable sum = 0.0
                                let mutable sumSq = 0.0

                                while queue.Count > 0 do
                                    let struct(vi,pi) = queue.HeapDequeue(cmp)
                                        
                                    let err = abs (info.Plane.Height (V2d pi + V2d.Half))
                                    if vi >= config.GrowThreshold && err <= config.Tolerance then
                                        sum <- sum + vi
                                        sumSq <- sumSq + sqr vi
                                        all.Add pi |> ignore
                                        r <- r.Add(V2d pi + V2d.Half, vi)
                                        if regressionStable r then
                                            match r.TryGetPlaneInfo() with
                                            | Some i -> info <- i
                                            | None -> ()
                                        for no in neighbours do
                                            let n = pi + no
                                            if n.AllGreaterOrEqual 0 && n.AllSmaller size then
                                                if queueHash.Add n then 
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
                                        let avg = sum / float all.Count
                                        //let avgSq = sumSq / float all.Count

                                        let var = (sumSq - float all.Count * sqr avg) / float (all.Count - 1)

                                        let line = Line2d(ray.GetPointOnRay tRange.Min, ray.GetPointOnRay tRange.Max)

                                        lines.Add {
                                            regression = r
                                            info = info
                                            avgResponse = avg
                                            stddev = sqrt var
                                            line = line
                                        }
                            | _ ->
                                ()
                    else
                        ()
                ) |> ignore

            )
            Log.stop()


            lines.ToArray()


    let private mergeLines (config : LineDetectionConfig) (lines : DetectedLine[])  =
        let radius = config.MergeRadius

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
                        

        let tryMerge (config : LineDetectionConfig) (l : Regression2d) (ll : Line2d) (r : Regression2d) (rl : Line2d) =
            
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
                        traverse config visited li { regression = hReg; info = hInfo; line = hLine; stddev = lInfo.stddev; avgResponse = lInfo.avgResponse }
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


    let run (config : LineDetectionConfig) (image : PixImage<byte>) =
       
        let rawLines = detectLinesPar config image
        
        Log.startTimed "merge lines"
        let merged = mergeLines config rawLines
        Log.stop()

        merged
        |> Array.filter (fun l -> l.info.Quality >= config.MinQuality)

    


