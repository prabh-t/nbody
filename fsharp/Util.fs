module Util

open System
open System.IO
open System.Diagnostics
open Microsoft.FSharp.Collections

// splitAt function similar as Haskell version
let splitAt n xs =  (Seq.truncate n xs, if Seq.length xs < n then Seq.empty else Seq.skip n xs)

// group the elements into chunks of size n
let rec chunk n xs =
    if Seq.isEmpty xs then Seq.empty
    else
        let (ys,zs) = splitAt n xs
        Seq.append (Seq.singleton ys) (chunk n zs)

let chunksOf n items =
    let rec loop i acc items = seq {
        match i, items, acc with
        //exit if chunk size is zero or input list is empty
        | _, [], [] | 0, _, [] -> ()
        //counter=0 so yield group and continue looping
        | 0, _, _::_ -> yield List.rev acc; yield! loop n [] items 
        //decrement counter, add head to group, and loop through tail
        | _, h::t, _ -> yield! loop (i-1) (h::acc) t
        //reached the end of input list, yield accumulated elements
        //handles items.Length % n <> 0
        | _, [], _ -> yield List.rev acc
    }
    loop n [] items

let flatten l =
    seq {
        yield Seq.head (Seq.head l) (* first item of first list *)
        for a in l do yield! (Seq.skip 1 a) (* other items *)
    }

let time msg f args =
    let stopwatch = System.Diagnostics.Stopwatch.StartNew()
    let temp = f args
    stopwatch.Stop()
    printfn "(%f s) %s: %A" stopwatch.Elapsed.TotalSeconds msg temp

let time2 msg f args =
    let stopwatch = System.Diagnostics.Stopwatch.StartNew()
    let temp = f args
    stopwatch.Stop()
    printfn "(%f s) %s" stopwatch.Elapsed.TotalSeconds msg
    temp

let inline myfold f a (xs: _ []) =
    let mutable a = a
    for i=0 to xs.Length-1 do
        a <- f a xs.[i]
    a

let rec lazyfoldl f z (l:LazyList<'a>) =
    match l with
    | LazyList.Nil -> z
    | LazyList.Cons (x,xs) ->
        let z' = f z x
        lazyfoldl f z' xs
