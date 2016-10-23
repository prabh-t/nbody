(*
 File: Allpairs.fs
 Description: nbody computation: barnes-hut algoritms
 Author: Prabhat Totoo
*)

module Allpairs

open System
open System.IO
open System.Diagnostics

open MapSkel

// (x,y,z,vx,vy,vz,m)
type Body = Body of double * double * double * double * double * double * double

// (ax,ay,az)
type Accel = Accel of double * double * double

let timeStep = 0.001
let eps = 0.01

let getBodies n generate =
    if generate then
        let newBody seed =
            let rnd = new Random(seed)
            let [x;y;z;vx;vy;vz;m] = List.init 7 (fun _ -> rnd.NextDouble())
            Body (x,y,z,vx,vy,vz,m)
        Seq.map newBody [0..(n-1)]
        |> Seq.toList
    else
        let bodyFromLine (line:string) =
            let words = line.Split [|' '|]
            let word (i:int) = double (Seq.nth i words)
            Body ((word 1),(word 2),(word 3),(word 4),(word 5),(word 6),(word 0))
        File.ReadLines("input.txt")
        |> Seq.map bodyFromLine
        |> Seq.toList

// checking function
let f acc (Body (x,y,z,vx,vy,vz,m)) = acc+x+y+z+vx+vy+vz+m

let mutable v = 0

let rec doSteps s bs =
    match s with
    | 0 -> bs
    | _ ->
        let updateVel (b:Body) =
            let accel (Body (ix,iy,iz,_,_,_,_) as bi) (Body (jx,jy,jz,_,_,_,jm) as bj) =
                if (bi = bj) then Accel (0.0,0.0,0.0)
                else
                    let dx = ix - jx
                    let dy = iy - jy
                    let dz = iz - jz
                    let dSquared = dx*dx + dy*dy + dz*dz + eps
                    let distance = sqrt (dSquared)
                    let mag = timeStep / (dSquared * distance)
                    let mmag = jm*mag
                    Accel ((dx*mmag),(dy*mmag),(dz*mmag))
            let deductChange (Body (x,y,z,vx,vy,vz,m)) (Accel (ax,ay,az)) = Body (x,y,z,(vx-ax),(vy-ay),(vz-az),m)
            
            //List.fold deductChange b (List.map (accel b) bs)  //mmxp230: 25.76s
            //fold/map merge
            List.fold (fun state b' -> deductChange state (accel b b')) b bs  //mmxp230: 14.45s

        let updatePos (Body (x,y,z,vx,vy,vz,m)) = Body ((x+timeStep*vx),(y+timeStep*vy),(z+timeStep* z),vx,vy,vz,m)
        let mapFunc =
            match v with
            | 0 -> List.map
            | 1 -> pmap_async_lst
            | 2 -> pmap_plinq_lst
            | _ -> List.map // default is sequential map 
        let new_bs = mapFunc (updatePos << updateVel) bs
        doSteps (s-1) new_bs

let run n steps pv =
    v <- pv
    printfn "ver: Allpairs \nparver: %d \nn: %d \nsteps: %d" v n steps
    let bs = getBodies n true
    let res1 = List.fold f 0.0 bs
    printfn "%f" res1
    let t1 = Stopwatch.StartNew()
    let bs' = doSteps steps bs
    let res2 = List.fold f 0.0 bs'
    printfn "%f" res2
    let t2 = t1.Elapsed.TotalSeconds
    printfn "time taken: %.2fs" t2