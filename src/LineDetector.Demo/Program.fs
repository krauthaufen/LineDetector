open Aardvark.Base
open LineDetector
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



    for (name, img) in testImages do
        Log.startTimed "detection %dx%d" img.Size.X img.Size.Y
        let cfg = { LineDetectionConfig.Default with GrowThreshold = 0.05 }
        let lines = LineDetector.run LineDetectionConfig.Default (img.ToPixImage<byte>())
        Log.line "%d lines found" lines.Length
        Log.stop()
    
        let img = img.ToPixImage<byte>(Col.Format.RGBA)
        let mat = img.GetMatrix<C4b>()

        let rand = RandomSystem()
        for li, d in Seq.indexed lines do
            let info = d.info
            let color = rand.UniformC3f().ToC4b()
            let n = d.line.Plane2d.Normal |> Vec.normalize
            for o in [V2d.Zero;-n * Constant.Sqrt2Half;n* Constant.Sqrt2Half] do 
                mat.SetLine(o+d.line.P0, o+d.line.P1, color)



        img.SaveAsImage (Path.Combine(outDir, name))

    0