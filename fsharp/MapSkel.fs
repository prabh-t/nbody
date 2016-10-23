module MapSkel

open System.Linq
open System.Threading
open System.Threading.Tasks
open System.Collections.Concurrent

open Util

let numProc = System.Environment.ProcessorCount

// Async Workflows
let pmap_async_arr f xs =
    seq { for x in xs -> async { return f x } }
    |> Async.Parallel
    |> Async.RunSynchronously

let pmap_async_seq f xs = pmap_async_arr f xs |> Array.toSeq

let pmap_async_lst f xs = pmap_async_arr f xs |> Seq.toList

// TPL

//// Tasks
let taskResult (t:Task<_>) = t.Result

////// Task:  does not return a result, but performs an imperative action

let pmap_tpl_tasks_arr_noreturn f (xs:array<_>) =
    let len = xs.Length
    //let tasks = Array.zeroCreate len
    for i in 0 .. len-1 do
        xs.[i] <- Task<_>.Factory.StartNew(fun () -> f xs.[i]).Result

////// Task<'T>: returns a result of type 'T

let pmap_tpl_tasks_arr f (xs:array<_>) =
    let tasks = xs |> Array.map (fun x -> Task<_>.Factory.StartNew(fun () -> f x).Result)
    tasks

let pmap_tpl_tasks_lst f (xs:list<_>) =
    let createTask x = Task<_>.Factory.StartNew(fun () -> f x).Result
    let tasks = xs |> List.map createTask
    tasks

let pmap_tpl_tasks_lstchunk f (xs:list<_>) =
    let chunks = chunksOf (xs.Length / (numProc * 2)) xs
    let chunkTask chunk = Task<_>.Factory.StartNew(fun () -> List.map f chunk).Result
    let tasks = List.map chunkTask (chunks |> Seq.toList)
    tasks |> List.concat
   
let pmap_tpl_tasks_lstparts f (xs:list<_>) =
    let partitions = Partitioner.Create(xs).GetPartitions(numProc)
    let tasks = [
        for partition in partitions do
            yield Task<_>.Factory.StartNew(fun () -> 
               [ while (partition.MoveNext()) do yield (f partition.Current) ]
            ).Result
    ]
    tasks |> List.concat

//loadBalance=true
let pmap_tpl_tasks_lstparts2 f (xs:list<_>) =
    let partitions = Partitioner.Create(xs.ToList(),true).GetPartitions(numProc)
    let tasks = [
        for partition in partitions do
            yield Task<_>.Factory.StartNew(fun () -> 
               [ while (partition.MoveNext()) do yield (f partition.Current) ]
            ).Result
    ]
    tasks |> List.concat

// not really sure what this does
let pmap_tpl_tasks_lstparts3 f (xs:list<_>) =
    let partitions = Partitioner.Create(xs).GetDynamicPartitions()
    let tasks = [
        for partition in partitions do
            yield Task<_>.Factory.StartNew(fun () -> (f partition)).Result
    ]
    tasks


//// Parallel.For

//// a bit of a workaround for list
let pmap_tpl_parfor_lst f (xs:list<_>) =
    let xs_arr = xs.ToArray()
    Parallel.For(0, xs_arr.Length, (fun i -> xs_arr.[i] <- f (xs_arr.[i]) )) |> ignore
    xs_arr |> Array.toList

//// for arrays

let pmap_tpl_parfor_noreturn f (xs:array<_>) =
    Parallel.For(0, xs.Length, (fun i -> xs.[i] <- f (xs.[i]) )) |> ignore

// return new_xs
let pmap_tpl_parfor f (xs:array<_>) =
    let new_xs = Array.zeroCreate xs.Length
    Parallel.For(0, xs.Length, (fun i -> new_xs.[i] <- f (xs.[i]) )) |> ignore
    new_xs

// return same (updated) input array: xs
let pmap_tpl_parfor2 f (xs:array<_>) =
    Parallel.For(0, xs.Length, (fun i -> xs.[i] <- f (xs.[i]) )) |> ignore
    xs


let pmap_tpl_parforchunk f (xs:array<_>) =
  let chunks = chunk 20000 [0..(Array.length xs - 1)]
  let noChunks = (Seq.length chunks)
  Parallel.For(0,noChunks,fun c ->
    let chunki = (Seq.nth c chunks)
    for i in chunki do
      xs.[i] <- f (xs.[i])
  ) |> ignore
  xs

let pmap_tpl_parfor_parts f (xs:array<_>) =
    let len = xs.Length
    let new_xs = Array.zeroCreate xs.Length
    let partitioner = Partitioner.Create(0,len-1,len / (numProc))
    Parallel.ForEach(partitioner, 
        (fun range loopState -> 
            let interval1, interval2 = fst range, snd range
            for i in interval1 .. interval2 do
                xs.[i] <- f (xs.[i])
        )) |> ignore
    xs


//// PLINQ - uses TPL internally

let pmap_plinq_seq f (xs:seq<_>) = xs.AsParallel().Select(fun x -> f x)
let pmap_plinq_lst f (xs:list<_>) = xs.AsParallel().Select(fun x -> f x) |> Seq.toList

let pmap_plinq_lst2 f (xs:list<_>) = xs.AsParallel()
                                                .WithDegreeOfParallelism(System.Environment.ProcessorCount * 2)
                                                .WithExecutionMode(ParallelExecutionMode.ForceParallelism)
                                                .Select(fun x -> f x) |> Seq.toList

let pmap_plinq_lstdyn f (xs:list<_>) = Partitioner.Create(xs).GetDynamicPartitions().AsParallel().Select(fun x -> f x) |> Seq.toList

let pmap_plinq_arr f (xs:array<_>) = xs.AsParallel().Select(fun x -> f x).ToArray()
