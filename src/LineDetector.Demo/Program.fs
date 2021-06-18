open Aardvark.Base
open LineDetector


[<EntryPoint>]
let main _args =
    Aardvark.Init()

    let img = PixImage.Create @"C:\Users\Schorsch\Desktop\Brusher\Datasets\orange\DSC02289.JPG"

    Log.startTimed "detecting lines"
    let lines = LineDetector.detectLines (img.ToPixImage<byte>())
    Log.line "%d lines found" lines.Length
    Log.stop()
    Log.startTimed "merging lines"
    let newLines = LineDetector.mergeLines lines
    Log.line "%d lines" lines.Length
    Log.stop()
    
    Log.startTimed "filtering lines"
    let newLines = newLines |> Array.filter (fun d -> d.info.Quality > 0.9 && d.info.Quality < 0.999)
    Log.line "%d lines" newLines.Length
    Log.stop()
    
    let img = img.ToPixImage<byte>(Col.Format.RGBA)
    let mat = img.GetMatrix<C4b>()

    let rand = RandomSystem()
    for li, d in Seq.indexed newLines do
        let info = d.info
        let color = rand.UniformC3f().ToC4b()
        let n = d.line.Plane2d.Normal |> Vec.normalize
        for o in [V2d.Zero] do 
            mat.SetLine(o+d.line.P0, o+d.line.P1, color)



    img.SaveAsImage @"C:\Users\Schorsch\Desktop\grad.png"

    0