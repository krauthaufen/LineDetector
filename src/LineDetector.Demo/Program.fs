open Aardvark.Base
open Aardvark.LineDetection
open System.IO

let trySeparate (m : Matrix<float>) =
    match SVD.Decompose m with
    | Some (u, s, vt) ->
        if Fun.IsTiny s.[1,1] then
            let f = sqrt s.[0,0]

            // (u00, u01)    *    (a, 0)     *     (v00, v01)
            // (u10, u11)         (0, 0)           (v10, v11)
            //
            // (u00, u01)    *    (a * v00, a * v01)  
            // (u10, u11)         (0,             0)  
            //
            // (u00*a*v00,    u00*a*v01)
            // (u10*a*v00,    u10*a*v01)
            //
            // a * (u00*v00,    u00*v01)
            //     (u10*v00,    u10*v01)
            //
            // a * outer(u.C0, v.R0)

            let f0 = vt.GetRow(0).Map (fun v -> f * v)
            let f1 = u.GetCol(0).Map (fun v -> f * v)
            Some(f0,f1)
        else
            None

    | None ->
        None


    // (a0, a1, a2, ...) * (b0, b1, b2, ..)^T = m

    //a0*b0 = m00
    //a1*b0 = m01

    //a[i]*b[j] = m[j][i]

    //a * b0 = m[0]


    //a0*b0 = -1
    //a1*b0 = -1
    //a2*b0 = -1

module HeatMap =
    let private heatMapColors =
        let fromInt (i : int) =
            C4b(
                byte ((i >>> 16) &&& 0xFF),
                byte ((i >>> 8) &&& 0xFF),
                byte (i &&& 0xFF),
                255uy
            )

        Array.map fromInt [|
            0x1639fa
            0x2050fa
            0x3275fb
            0x459afa
            0x55bdfb
            0x67e1fc
            0x72f9f4
            0x72f8d3
            0x72f7ad
            0x71f787
            0x71f55f
            0x70f538
            0x74f530
            0x86f631
            0x9ff633
            0xbbf735
            0xd9f938
            0xf7fa3b
            0xfae238
            0xf4be31
            0xf29c2d
            0xee7627
            0xec5223
            0xeb3b22
        |]

    let color (tc : float) =
        let tc = clamp 0.0 1.0 tc
        let fid = tc * float heatMapColors.Length - 0.5

        let id = int (floor fid)
        if id < 0 then 
            heatMapColors.[0]
        elif id >= heatMapColors.Length - 1 then
            heatMapColors.[heatMapColors.Length - 1]
        else
            let c0 = heatMapColors.[id].ToC3f()
            let c1 = heatMapColors.[id + 1].ToC3f()
            let t = fid - float id
            (c0 * (1.0f - float32 t) + c1 * float32 t).ToC4b()



[<EntryPoint>]
let main _args =
    
    //let filter =
    //    let s = V2l.II * 5L
    //    let m = Matrix<float>(s)
    //    let center = s / 2L
    //    m.SetByCoord (fun (c : V2l) ->
    //        let d = V2d (c - center) / (0.5 * V2d s)

    //        exp -d.LengthSquared

    //    )

    //let filter =
    //    Matrix<float>([|
    //        0.0; -1.0; 0.0
    //        -1.0;  4.0; -1.0
    //        0.0; -1.0;  0.0
    //    |], 3L, 3L)

    //match trySeparate filter with
    //| Some (a, b) ->

    //    let m = Matrix<float>(V2l(a.Size, b.Size))
    //    m.SetByCoord(fun (c : V2l) ->
    //        a.[c.Y]*b.[c.X]
    //    ) |> ignore

    //    printfn "original: "
    //    for y in 0L .. filter.SY - 1L do
    //        for x in 0L .. filter.SX - 1L do
    //            printf "% .3f " filter.[x,y]
    //        printfn ""

    //    printfn "separate: "

    //    printfn " a: %s" (Array.init (int a.S) (fun i -> sprintf "% .3f" a.[i]) |> String.concat " ")
    //    printfn " b: %s" (Array.init (int b.S) (fun i -> sprintf "% .3f" b.[i]) |> String.concat " ")

    //    printfn "outer(a,b): "

    //    for y in 0L .. m.SY - 1L do
    //        for x in 0L .. m.SX - 1L do
    //            printf "% .3f " m.[x,y]
    //        printfn ""

    //    //let a = V3d a
    //    //let b = V3d b
    //    //Log.warn "%A" a
    //    //Log.warn "%A" b

    //    //let test = Vec.Outer(a,b)
    //    //Log.warn "%A" test

    //    ()
    //| None ->
    //    ()

    //exit 0


    Aardvark.Init()

    let inDir = @"C:\Users\Schorsch\Desktop\LineDetector"
    let outDir = @"C:\Users\Schorsch\Desktop\LineDetector\results"

    let testImages =
        Directory.GetFiles inDir
        |> Array.choose (fun p ->
            try 
                let img = PixImage.Create p
                Some(Path.GetFileName p, img)
            with _ -> None
        )


     //0: detection 7952x5304 ............................................. 10.782 s
     //0: detection 4912x3264 .............................................. 1.400 s
     //0: detection 4912x3264 .............................................. 1.143 s
     //0: detection 4928x3264 .............................................. 1.184 s
     //0: detection 1600x1200 .............................................. 0.581 s
     //0: detection 4000x1800 .............................................. 2.295 s
     //0: detection 6000x4000 .............................................. 2.296 s
     //0: detection 6000x4000 .............................................. 2.085 s
     //0: detection 5472x3648 .............................................. 1.655 s
     //0: detection 4608x2074 .............................................. 0.739 s


    for (name, img) in testImages do
        let img = img.ToPixImage<byte>()
        Log.startTimed "detection %dx%d" img.Size.X img.Size.Y
        let cfg = { LineDetectionConfig.Default with MinQuality = 0.8 }
        let lines = LineDetector.detect cfg Unchecked.defaultof<_> img
        //Log.line "%d lines found" lines.Length
        let lines = 
            lines |> Array.filter (fun d ->
                d.info.AngularError * Constant.DegreesPerRadian < 10.0
            )
        //Log.line "%d lines" lines.Length
        Log.stop()

        let img = img.ToPixImage<byte>(Col.Format.RGBA)
        let mat = img.GetMatrix<C4b>()

        let rand = RandomSystem()
        for li, d in Seq.indexed lines do
            let info = d.info
            let color = rand.UniformC3f().ToC4b()
            //let q = abs d.info.StdDev.Y / cfg.Tolerance
            //let color = HeatMap.color q
            let n = d.line.Plane2d.Normal |> Vec.normalize
            //let color = HeatMap.color (d.info.AngularError * Constant.DegreesPerRadian / 15.0)
            for o in [V2d.Zero;-n * Constant.Sqrt2Half;n* Constant.Sqrt2Half] do 
                mat.SetLine(o+d.line.P0, o+d.line.P1, color)



        img.SaveAsImage (Path.Combine(outDir, name))

    0