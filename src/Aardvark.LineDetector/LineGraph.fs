namespace Aardvark.LineDetection

open Aardvark.Base
open Aardvark.Geometry
open FSharp.Data.Adaptive

type DetectedLine =
    {
        regression  : WeightedRegression2d
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
            
        //let linePoints = Array.zeroCreate (2 * lines.Length)
        //let mutable oi = 0
        //for i in 0 .. lines.Length - 1 do
        //    let { line = l } = lines.[i]
        //    linePoints.[oi] <- l.P0
        //    linePoints.[oi+1] <- l.P1
        //    oi <- oi + 2

        let graph = LineGraph()

        
        
        let mutable tree = Aardvark.Geometry.BvhTree2d.Empty 8
        for li in 0 .. lines.Length - 1 do
            let { line = l } = lines.[li]
            tree <- tree.Add(li, l.BoundingBox2d, l)
            

        let query (li : int) (l : Line2d) =
            let bb = l.BoundingBox2d
            let q = Box2d(bb.Min - V2d radius, bb.Max + V2d radius)
            tree.GetIntersecting(q, fun ri _bb rl ->
                if ri > li then
                    let d = l.GetMinimalDistanceTo rl
                    if d <= radius then Some ri
                    else None
                else
                    None
            )
            |> HashMap.keys

        for li in 0 .. lines.Length - 1 do
            let { line = l } = lines.[li]
            let lp = l.Plane2d

            for ri in query li l do
                
                let { line = r } = lines.[ri]

                let d = r.GetMinimalDistanceTo l
                if d <= radius then

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
                        let rp = r.Plane2d
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
