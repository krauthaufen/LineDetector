namespace Aardvark.LineDetection

open Aardvark.Base
open System

module internal Dither =
    type private Marker = class end

    let matrix =
        let ass = typeof<Marker>.Assembly

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

        let res = Array.zeroCreate (size.X * size.Y)
        let mutable index = 0

        let inline add v = 
            res.[index] <- v
            index <- index + 1

        for by in 0 .. blocks.Y-1 do
            for bx in 0 .. blocks.X-1 do
                let o = V2i(bx, by) * V2i matrix.Size
                let e = o + V2i matrix.Size
                if e.AllSmallerOrEqual size then
                    if o.X = 0 && o.Y = 0 then 
                        coords.CopyTo(res, index)
                        index <- index + coords.Length
                    else 
                        for c in coords do add (c+o)
                else
                    for c in coords do
                        let cc = o + c
                        if cc.AllSmaller size then add cc
                
        res
