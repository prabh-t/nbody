(*
 File: BHList.fs
 Description: nbody computation: barnes-hut algoritms
 Author: Prabhat Totoo
*)

module BHList

open System
open System.IO
open System.Diagnostics

open MapSkel

// a body consists of pos, vel and mass
// (x,y,z,vx,vy,vz,m)
type Body = Body of double * double * double * double * double * double * double

// centroid
// (cx,cy,cz,cm)
type Centroid = Centroid of double * double * double * double

// acceleration
// (ax,ay,az)
type Accel = Accel of double * double * double

// bounding box - the min and max pts in 3D space, making a region containing points
// (minx,miny,minz,maxx,maxy,maxz)
type Bbox = Bbox of double * double * double * double * double * double

// the BH tree
// node consists of size, center X, Y, Z, total mass, and subtrees
// (size,centerx,centery,centerz,mass,subTrees)
type BHTree = BHT of double * double * double * double * double * BHTree list


let timeStep = 0.001

// If the distance between the points is smaller than this
let eps = 0.01     // then ignore the forces between them.

// threshold is the ratio of size of a region and the distance between
// a body and the region centre of mass.
// If s / d < threshold, then the internal node is sufficiently far away.
let threshold = 0.25 // threshold=0 degenerates to brute force

let getBodies n generate =
    if generate then
        let newBody seed =
            let rnd = new Random(seed)
            let [x;y;z;vx;vy;vz;m] = List.init 7 (fun _ -> rnd.NextDouble())
            Body (x,y,z,vx,vy,vz,m)
        List.map newBody [0..(n-1)]
    else
        let bodyFromLine (line:string) =
            let words = line.Split [|' '|]
            let word (i:int) = double (Seq.nth i words)
            Body ((word 1),(word 2),(word 3),(word 4),(word 5),(word 6),(word 0))
        File.ReadLines("input.txt")
        |> Seq.map bodyFromLine
        |> Seq.toList

// find the coordinates of the bounding box that contains the bodies
let inline findBounds bs =
    let f (Bbox (minx,miny,minz,maxx,maxy,maxz)) (Body (x,y,z,_,_,_,_)) =
        Bbox ((min minx x),(min miny y),(min minz z),(max maxx x),(max maxy y),(max maxz z))
    List.fold f (Bbox (0.0,0.0,0.0,0.0,0.0,0.0)) bs
    
// calculate the centroid of points
let inline calcCentroid bs =
    let f (Centroid (x,y,z,m)) (Body (x',y',z',_,_,_,m')) = Centroid ((m'*x'+x),(m'*y'+y),(m'*z'+z),(m'+m))
    let (Centroid (xs,ys,zs,mass)) = List.fold f (Centroid (0.0,0.0,0.0,0.0)) bs
    Centroid ((xs/mass),(ys/mass),(zs/mass),mass)


// check if point is in box
let inline inBox (Bbox (minx,miny,minz,maxx,maxy,maxz)) (Body (x,y,z,_,_,_,_)) =
    (x > minx) && (x <= maxx) && (y > miny) && (y <= maxy) && (z > minz) && (z <= maxz)
    
// split points according to their locations in the box
let inline splitPoints (Bbox (minx,miny,minz,maxx,maxy,maxz) as bb) bs =
    match bs with
    | [] -> [(bb,[])]
    | b1::[] -> [(bb,[b1])]
    | _ -> let (midx,midy,midz) = ((minx+maxx)/2.0,(miny+maxy)/2.0,(minz+maxz)/2.0)
           let box1 = Bbox (minx,miny,minz,midx,midy,midz)
           let box2 = Bbox (minx,midy,minz,midx,maxy,midz)
           let box3 = Bbox (midx,miny,minz,maxx,midy,midz)
           let box4 = Bbox (midx,midy,minz,maxx,maxy,midz)
           let box5 = Bbox (minx,miny,midz,midx,midy,maxz)
           let box6 = Bbox (minx,midy,midz,midx,maxy,maxz)
           let box7 = Bbox (midx,miny,midz,maxx,midy,maxz)
           let box8 = Bbox (midx,midy,midz,maxx,maxy,maxz)
           let p1 = List.filter (inBox box1) bs
           let p2 = List.filter (inBox box2) bs
           let p3 = List.filter (inBox box3) bs
           let p4 = List.filter (inBox box4) bs
           let p5 = List.filter (inBox box5) bs
           let p6 = List.filter (inBox box6) bs
           let p7 = List.filter (inBox box7) bs
           let p8 = List.filter (inBox box8) bs
           let boxes = box1::box2::box3::box4::box5::box6::box7::box8::[]
           let splitPts = p1::p2::p3::p4::p5::p6::p7::p8::[]
           let boxesAndPts = List.filter (fun (_,pts) -> not (List.isEmpty pts)) (List.zip boxes splitPts)
           boxesAndPts

// build the Barnes-Hut tree
let rec buildTree (Bbox (minx,miny,minz,maxx,maxy,maxz) as bb) bs =
    let size = List.min [abs (maxx - minx); abs (maxy - miny); abs (maxz - minz)]
    let (Centroid (centerx,centery,centerz,mass)) = calcCentroid bs
    match bs with
    | [] -> BHT (size,centerx,centery,centerz,mass,[])
    | b1::[] -> BHT (size,centerx,centery,centerz,mass,[])
    | _ ->
        let boxesAndPts = splitPoints bb bs
        let subTrees = List.map (fun (bb',bs') -> buildTree bb' bs') boxesAndPts
        BHT (size,centerx,centery,centerz,mass,subTrees)

// calculate the acceleration on a point due to some other point
let inline accel (BHT (_,centerx,centery,centerz,mass,_)) (Body (x,y,z,_,_,_,m)) =
    let dx = x - centerx
    let dy = y - centery
    let dz = z - centerz
    let dsqr = (dx * dx) + (dy * dy) + (dz * dz) + eps
    let d = sqrt dsqr 
    let mag = timeStep / (dsqr * d)
    let mmag = mass * mag
    Accel ((dx * mmag),(dy * mmag),(dz * mmag))

// use centroid as approximation if the point is far from a cell
let inline isFar (BHT (size,centerx,centery,centerz,_,_)) (Body (x,y,z,_,_,_,_)) =
    let dx = centerx - x
    let dy = centery - y
    let dz = centerz - z
    let dist = sqrt (dx * dx + dy * dy + dz * dz)
    (size / dist) < threshold

let inline addAccel (Accel (ax1,ay1,az1)) (Accel (ax2,ay2,az2)) = Accel ((ax1+ax2),(ay1+ay2),(az1+az2))

// calculate the accelleration of a point due to the points in the given tree
let rec calcAccel (BHT (_,_,_,_,_,subtrees) as t) (b:Body) =
    match subtrees with
    | [] -> accel t b
    | _ ->
        if (isFar t b) then accel t b
        else 
            List.fold addAccel (Accel (0.0,0.0,0.0)) (List.map (fun st -> 
                    calcAccel st b) subtrees)

// this version merge the fold and map as Tomas suggested (performs better than the other)
let rec calcAccel1 (BHT (_,_,_,_,_,subtrees) as t) (b:Body) =
    match subtrees with
    | [] -> accel t b
    | _ ->
        if (isFar t b) then accel t b
        else 
            List.fold (fun acc1 st -> addAccel acc1 (calcAccel1 st b)) (Accel (0.0,0.0,0.0)) subtrees

// this version uses an inner recursive function that is tail optimised.
// the initial calcAccel code is changed in order to make it possible/easy
// to make tail recursive.
let calcAccel2 (t:BHTree) (b:Body) =
    let rec calcAccelRec (subtrees:BHTree list) (Accel (ax1,ay1,az1) as acc) =
        match subtrees with
        | [] -> acc
        | (BHT (_,_,_,_,_,t1sub) as t1)::ts ->
            if (List.isEmpty t1sub || isFar t1 b) then 
                let (Accel (ax2,ay2,az2)) = accel t1 b
                calcAccelRec ts (Accel ((ax1+ax2),(ay1+ay2),(az1+az2)))
            else calcAccelRec (List.append t1sub ts) acc
    calcAccelRec [t] (Accel (0.0,0.0,0.0))

// continuation-based
let accelCont (t:BHTree) (b:Body) cont = cont (accel t b)

let addAccelCont (acc1:Accel) (acc2:Accel) cont = cont (addAccel acc1 acc2)

let rec foldCont f (state:Accel) (list:BHTree list) cont =
    match list with
    | [] -> cont state
    | x::xs -> foldCont f state xs (fun state -> f x state cont)

let rec calcAccelCont (BHT (_,_,_,_,_,subtrees) as t) (b:Body) cont =
    match subtrees with
    | [] -> accelCont t b cont
    | _ ->
        if (isFar t b) then accelCont t b cont
        else
            foldCont (fun st state cont -> 
                calcAccelCont st b (fun res ->
                addAccelCont state res cont)) (Accel (0.0,0.0,0.0)) subtrees cont

// checking function
let f acc (Body (x,y,z,vx,vy,vz,m)) = acc+x+y+z+vx+vy+vz+m

let mutable v = 0

let rec doSteps s bs =
    match s with
    | 0 -> bs
    | _ ->
        let box = findBounds bs
        let tree = buildTree box bs
        let updateVel (Body (x,y,z,vx,vy,vz,m)) =
            let (Accel (ax,ay,az)) = calcAccel1 tree (Body (x,y,z,vx,vy,vz,m)) //(fun a -> a)
            Body (x,y,z,(vx-ax),(vy-ay),(vz-az),m)
        let updatePos (Body (x,y,z,vx,vy,vz,m)) = Body ((x + timeStep * vx),(y + timeStep * vy),(z + timeStep * vz),vx,vy,vz,m)
        let mapFunc =
            match v with
            | 0 -> List.map
            // in order in which they are presented in the paper
            
            // Async Workflow
            | 1 -> pmap_async_lst
            
            // explicit tasks creation/with chunking (custom)/partitions (inbuilt)
            | 2 -> pmap_tpl_tasks_lst
            | 3 -> pmap_tpl_tasks_lstchunk
            | 4 -> pmap_tpl_tasks_lstparts
            | 44 -> pmap_tpl_tasks_lstparts2 //loadBalance=true
            | 444 -> pmap_tpl_tasks_lstparts3 //dynamicPartitions
            
            // PLINQ
            | 5 -> pmap_plinq_lst
            | 55 -> pmap_plinq_lst2
            | 6 -> pmap_plinq_lstdyn

            // Parallel.For (convert list to array, then back)
            | 7 -> pmap_tpl_parfor_lst

            | _ -> List.map // default is sequential map 
        let new_bs = mapFunc (updatePos << updateVel) bs
        doSteps (s-1) new_bs

let run n steps pv =
    v <- pv
    printfn "ver: BHList \nparver: %d \nn: %d \nsteps: %d" v n steps
    let bs = getBodies n true
    let res1 = List.fold f 0.0 bs
    printfn "%f" res1
    let t1 = Stopwatch.StartNew()
    let bs' = doSteps steps bs
    let res2 = List.fold f 0.0 bs'
    printfn "%f" res2
    let t2 = t1.Elapsed.TotalSeconds
    printfn "time taken: %.2fs" t2
