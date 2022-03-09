open System.IO
open Aardvark.Base
open Aardvark.LineDetection
open System.Threading
open System
open System.Text

type Args =
    {
        Inputs      : list<string>
        Output      : string
        SaveImages  : bool
    }

module Args =
    let rec private run (acc : Args) (l : list<string>) =    
        match l with
        | "-o" :: path :: rest ->
            run { acc with Output = path } rest

        | "-i" :: rest ->
            run { acc with SaveImages = true } rest

        | a :: rest when a.StartsWith "-" ->
            Log.warn "unknown option %A" a
            run acc rest
        | path :: rest ->
            let files = 
                if File.Exists path then 
                    [path]
                elif Directory.Exists path then 
                    let imageExts = Set.ofList [".jpg"; ".png"]
                    Directory.GetFiles(path)
                    |> Array.filter (fun p -> Set.contains (Path.GetExtension p) imageExts)
                    |> Array.toList
                else
                    Log.warn "bad input: %A" path
                    []

            run { acc with Inputs = acc.Inputs @ files } rest
        | [] ->
            acc
            
    let tryParse (a : string[]) =
        let args = run { Inputs = []; Output = null; SaveImages = false } (Array.toList a)

        if not (List.isEmpty args.Inputs) then
            let args =
                if isNull args.Output then { args with Output = Environment.CurrentDirectory }
                else args

            Some args
        else
            None

let help() = ()
    

[<EntryPoint>]
let main argv =
    match Args.tryParse argv with
    | Some a ->
        
        for f in a.Inputs do
            try
                let img = PixImageSharp.Create(f).ToPixImage<byte>()
                let lines = LineDetector.detect LineDetectionConfig.Default CancellationToken.None img
            
                let json = StringBuilder()
                json.AppendLine "[" |> ignore

                for l in lines do
                    let p0 = l.line.P0
                    let p1 = l.line.P1
                    json.AppendLine(
                        $"    {{ X0: %.3f{p0.X}, Y0: %.3f{p0.Y}, X1: %.3f{p1.X}, Y1: %.3f{p1.Y} }},"
                    ) |> ignore


                json.AppendLine "]" |> ignore

                let outfile = Path.Combine(a.Output, Path.ChangeExtension(Path.GetFileName(f), ".json"))
                File.WriteAllText(outfile, json.ToString())

            with _ ->
                Log.error "could not load %s" f
    | None ->
        help()
    0