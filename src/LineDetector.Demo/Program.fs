open Aardvark.Base
open LineDetector


[<EntryPoint>]
let main _args =
    Aardvark.Init()

    let img = PixImage.Create @"C:\Users\Schorsch\Desktop\Brusher\Datasets\H1_29\IMG_7761.JPG"
    //let img = 
    //    let dst = PixImage<byte>(Col.Format.Gray, img.Size)
    //    let src = img.ToPixImage<byte>().GetMatrix<C4b>()
    //    dst.GetChannel(0L).SetMap(src, System.Func<_,_>(fun (c : C4b) -> c.ToGrayByte())) |> ignore
    //    dst

    Log.startTimed "detection %dx%d" img.Size.X img.Size.Y
    Log.startTimed "detecting lines"
    let lines = LineDetector.detectLines (img.ToPixImage<byte>())
    Log.line "%d lines found" lines.Length
    Log.stop()
    Log.startTimed "merging lines"
    let newLines = LineDetector.mergeLines 20.0 lines
    Log.line "%d lines" newLines.Length
    Log.stop()
    
    Log.startTimed "filtering lines"
    let newLines = newLines |> Array.filter (fun d ->  Vec.distance d.line.P0 d.line.P1 >= 40.0 && d.info.Quality > 0.7 && d.info.Quality < 0.999)
    Log.line "%d lines" newLines.Length
    Log.stop()
    Log.stop()
    
    let img = img.ToPixImage<byte>(Col.Format.RGBA)
    let mat = img.GetMatrix<C4b>()

    let rand = RandomSystem()
    for li, d in Seq.indexed newLines do
        let info = d.info
        let color = rand.UniformC3f().ToC4b()
        let n = d.line.Plane2d.Normal |> Vec.normalize
        for o in [V2d.Zero;-n;n] do 
            mat.SetLine(o+d.line.P0, o+d.line.P1, color)



    img.SaveAsImage @"C:\Users\Schorsch\Desktop\grad.jpg"

    0