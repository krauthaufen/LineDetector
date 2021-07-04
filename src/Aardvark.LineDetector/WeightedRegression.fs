namespace Aardvark.LineDetection

open Aardvark.Base
open System

#nowarn "9"

[<Struct>]
type PlaneFitInfo =
    {
        Trafo           : Trafo2d
        Plane           : Plane2d
        StdDev          : V2d
        Center          : V2d
        Count           : int
        AverageMass     : float
        MassStdDev      : float
        AngularError    : float
    }

    member x.Quality =
        let r = abs x.StdDev.X / (abs x.StdDev.X + abs x.StdDev.Y) 
        2.0 * r - 1.0

    member x.XAxis = x.Trafo.Forward.C0.XY
    member x.YAxis = x.Trafo.Forward.C1.XY

    member x.ToWeightedRegression2d() =

        let mass = x.AverageMass * float x.Count
        
        let massSq =
            sqr x.MassStdDev * float (x.Count - 1) + float x.Count * sqr x.AverageMass

        let s = M22d.FromDiagonal(sqr x.StdDev * mass)
        let u = x.Trafo.Forward.UpperLeftM22()
        let i = u * s * u.Transposed

        let sumSq = V3d(i.M00, i.M11, i.M01)
        WeightedRegression2d(sumSq, V2d.Zero, x.Center, x.Count, mass, massSq)


/// Represents an immutable, incremental regression-plane with per-point weights
and [<Struct>] WeightedRegression2d(sumSq : V3d, sum : V2d, ref : V2d, count : int, mass : float, massSq : float) =
        
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

            let la0 = Vec.lengthSquared a0
            let la1 = Vec.lengthSquared a1
            let x = 
                if la0 > la1 then a0 / sqrt la0
                else a1 / sqrt la1
            let y = x.Rot90

            struct(M22d.FromCols(x,y), V2d(a, b))


    /// The empty regression without any points
    static member Empty = WeightedRegression2d(V3d.Zero, V2d.Zero, V2d.Zero, 0, 0.0, 0.0)
     

            
    /// Adds a point to the regression and returns the result
    member x.Add(point : V2d, m : float) =
        if m <= 0.0 then
            x
        elif count = 0 then
            WeightedRegression2d(V3d.Zero, V2d.Zero, point, 1, m, sqr m)
        else
            let p = point - ref
            let sumSq = sumSq + m * V3d(sqr p.X, sqr p.Y, p.X * p.Y)
            let sum = sum + m * p
            WeightedRegression2d(sumSq, sum, ref, count + 1, mass + m, massSq + sqr m)

    /// Removes a point from the regression and returns the result.
    /// NOTE that the implementation assumes that this point has previously been added 
    ///      and may result in invalid state whennot.
    member x.Remove(point : V2d, m : float) =
        if m <= 0.0 then
            x
        elif count = 1 then
            WeightedRegression2d.Empty
        else
            let p = point - ref
            let sumSq = sumSq - m * V3d(sqr p.X, sqr p.Y, p.X * p.Y)
            let sum = sum - m * p
            WeightedRegression2d(sumSq, sum, ref, count - 1, mass - m, massSq - sqr m)

    /// Checks if the regression is empty
    member x.IsEmpty =
        count = 0

    /// The weighted average of all points (raises IndexOutOfRangeException it empty)
    member x.Center = 
        if count < 1 then raise <| IndexOutOfRangeException()
        if count = 1 then ref
        else ref + sum / mass

    /// The total "mass" of all points added so far
    member x.Mass = mass

    /// The number of points added so far
    member x.Count = count

    /// The inertia tensor for all points added
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
                
    /// Gets a quality measure in [0..1] based on comparing the "deviations" along the principal axes
    member x.GetQuality() =
        if count >= 2 then
            let s = getEigenvalues sumSq sum mass
            let stddev = abs s / mass |> sqrt
            let r = abs stddev.X / (abs stddev.X + abs stddev.Y) 
            2.0 * r - 1.0
        else
            0.0
            
    /// Gets the (optional) current regression plane including several statistics (AngularError, etc.)
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

                let avgMass = mass / float count
                let massStdDev = (massSq - float count * sqr avgMass) / float (count - 1) |> sqrt


                Some {
                    Plane = Plane2d(y, avg)
                    Trafo = Trafo2d.FromBasis(x, y, avg)
                    StdDev = abs s / mass |> sqrt
                    Center = avg
                    Count = count
                    AverageMass = avgMass
                    MassStdDev = massStdDev
                    AngularError = angle
                }
        else
            None
                
    member private x.Ref = ref
    member private x.Sum = sum
    member private x.SumSq = sumSq
    member private x.MassSq = massSq

    static member Zero = WeightedRegression2d.Empty

    static member Union (l : WeightedRegression2d, r : WeightedRegression2d) =
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

            WeightedRegression2d(sumSq, sum, ref, cnt, mass, l.MassSq + r.MassSq)
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

            WeightedRegression2d(sumSq, sum, ref, cnt, mass, l.MassSq + r.MassSq)
        
    static member Difference (l : WeightedRegression2d, r : WeightedRegression2d) =
        if l.IsEmpty then l
        elif r.IsEmpty then l
        elif r.Count = 1 then l.Remove(r.Center, r.Mass)
        elif l.Count >= r.Count then 
            let ref = l.Ref
            let o = r.Ref - ref
            let cnt = l.Count - r.Count
            let mass = l.Mass - r.Mass
            let sum = l.Sum - r.Sum - r.Mass * o

            let sumSq =
                let m = r.Sum
                l.SumSq - V3d(
                    r.SumSq.X   + 2.0*o.X*m.X       + r.Mass*sqr o.X,
                    r.SumSq.Y   + 2.0*o.Y*m.Y       + r.Mass*sqr o.Y,
                    r.SumSq.Z   + o.X*m.Y + o.Y*m.X + r.Mass*o.X*o.Y
                )

            WeightedRegression2d(sumSq, sum, ref, cnt, mass, l.MassSq - r.MassSq)
        else
            WeightedRegression2d.Empty

    static member (+) (l : WeightedRegression2d, r : WeightedRegression2d) =
        WeightedRegression2d.Union(l, r)
    
    static member (-) (l : WeightedRegression2d, r : WeightedRegression2d) =
        WeightedRegression2d.Difference(l, r)


module WeightedRegression2d =
    
    /// The empty regression.
    let empty = WeightedRegression2d.Empty

    /// Creates a regression from the given points and masses.
    let inline ofSeq (s : seq<V2d * float>) =
        let mutable res = empty
        for (pt, m) in s do res <- res.Add(pt, m)
        res
        
    /// Creates a regression from the given points and masses.
    let inline ofSeqV (s : seq<struct(V2d * float)>) =
        let mutable res = empty
        for struct(pt, m) in s do res <- res.Add(pt, m)
        res
        
    /// Creates a regression from the given points and masses.
    let inline ofList (s : list<V2d * float>) =
        let mutable res = empty
        for (pt, m) in s do res <- res.Add(pt, m)
        res
        
    /// Creates a regression from the given points and masses.
    let inline ofListV (s : list<struct(V2d * float)>) =
        let mutable res = empty
        for struct(pt, m) in s do res <- res.Add(pt, m)
        res
        
    /// Creates a regression from the given points and masses.
    let inline ofArray (s : array<V2d * float>) =
        let mutable res = empty
        for (pt, m) in s do res <- res.Add(pt, m)
        res
        
    /// Creates a regression from the given points and masses.
    let inline ofArrayV (s : array<struct(V2d * float)>) =
        let mutable res = empty
        for struct(pt, m) in s do res <- res.Add(pt, m)
        res

    /// Adds a point to the regression.
    let inline add (point : V2d) (mass : float) (regression : WeightedRegression2d) =
        regression.Add(point, mass)
        
    /// Removes a point from the regression and returns the result.
    /// NOTE that the implementation assumes that this point has previously been added 
    ///      and may result in invalid state whennot.
    let inline remove (point : V2d) (mass : float) (regression : WeightedRegression2d) =
        regression.Remove(point, mass)
        
    /// Gets the (optional) current regression plane including several statistics (AngularError, etc.)
    let inline tryGetPlaneFitInfo (r : WeightedRegression2d) =
        r.TryGetPlaneInfo()




