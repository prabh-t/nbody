{-# LANGUAGE CPP, BangPatterns, NoMonomorphismRestriction #-}

-- File: allpairs.hs
-- Description: nbody computation: full all-pairs version
-- Prabhat Totoo

import System.IO
import System.Environment
import System.Time
import System.Random
import Text.Printf
import Data.List

#if defined (PAR)
import Control.DeepSeq
import GHC.Conc (numCapabilities)
import Control.Parallel
import Control.Parallel.Strategies
#endif

-- a body consists of pos, vel and mass
data Body
	= Body
	{ x	 :: {-# UNPACK #-} !Double	-- pos of x
	, y	 :: {-# UNPACK #-} !Double	-- pos of y
	, z	 :: {-# UNPACK #-} !Double	-- pos of z
	, vx :: {-# UNPACK #-} !Double	-- vel of x
	, vy :: {-# UNPACK #-} !Double	-- vel of y
	, vz :: {-# UNPACK #-} !Double	-- vel of z
	, m :: {-# UNPACK #-} !Double } -- mass
	deriving (Show,Eq)

-- acceleration
data Accel
	= Accel
	{ ax :: {-# UNPACK #-} !Double
	, ay :: {-# UNPACK #-} !Double
	, az :: {-# UNPACK #-} !Double }

#if defined (PAR)
instance NFData Body
instance NFData Accel
#endif

timeStep = 0.001

-- epsilon
eps = 0.01

randomList :: (Random a) => Int -> [a]
randomList seed = randoms (mkStdGen seed)
   
genBody s = Body x y z vx vy vz m
	where [x,y,z,vx,vy,vz,m] = take 7 $ randomList s

-- func to calc time taken between t0 and t1 in sec
secDiff::ClockTime -> ClockTime -> Float
secDiff (TOD secs1 psecs1) (TOD secs2 psecs2) = fromInteger (psecs2 - psecs1) / 1e12 + fromInteger (secs2 - secs1)

main =
	do
		args <- getArgs
		let 
			n = read (args!!0) :: Int -- num of bodies
			steps = read (args!!1) :: Int -- num of steps (iterations)
			bs = [ (genBody i) | i <- [0..n-1] ]
			res1 = foldl' f 0 bs -- do meaningless sum...
		print (res1)
		t1 <- getClockTime
		let bs' = doSteps steps bs
		print $ foldl' f 0 bs'
		t2 <- getClockTime
		putStrLn $ printf "time taken: %.2fs" $ secDiff t1 t2

f acc (Body x y z vx vy vz m) = acc+x+y+z+vx+vy+vz+m

doSteps 0 bs = bs
doSteps s bs = doSteps (s-1) new_bs
	where
#if defined (PAR)
#if defined (NOCHUNK)
		new_bs = map (updatePos . updateVel) bs `using` parList rdeepseq  -- no chunking
#else
		new_bs = map (updatePos . updateVel) bs `using` parListChunk chunksize rdeepseq
		chunksize = (length bs) `quot` (numCapabilities * 4)
#endif
#else
		new_bs = map f bs
#endif
		f = (updatePos . updateVel)
		updatePos (Body x y z vx vy vz m) = Body (x+timeStep*vx) (y+timeStep*vy) (z+timeStep*vz) vx vy vz m

		--runtimes on mmxp230: 16k,1iter,3runs
		
		--(1)as in journal paper: 15.55s
		updateVel b = foldl' deductChange b (map (accel b) bs)
		
		--foldl: 15.77s
		--updateVel b = foldl deductChange b (map (accel b) bs)
		--foldr: 16.14s
		--updateVel b = foldr deductChangeR b (map (accel b) bs)
		
		--(2)use list comprehension instead of map: 16.85s
		--updateVel b = foldl' deductChange b [(accel b b')| b' <- bs] 
		
		--foldr/build: 16.18s
		--updateVel b = foldr deductChangeR b [(accel b b')| b' <- bs]
		
		--(3)foldl/map merge (manual): 52.43s
		--updateVel b = foldl' (\state b' -> deductChange state (accel b b')) b bs
		
		--foldr/map merge (manual): 26.43s
		--updateVel b = foldr (\state b' -> deductChange state (accel b b')) b bs

		deductChange (Body x y z vx vy vz m) (Accel ax ay az) = Body x y z (vx-ax) (vy-ay) (vz-az) m
		deductChangeR (Accel ax ay az) (Body x y z vx vy vz m) = Body x y z (vx-ax) (vy-ay) (vz-az) m

		accel b_i b_j
			| b_i == b_j = Accel 0 0 0
			| otherwise = Accel (dx*mmag) (dy*mmag) (dz*mmag)
				where
					mmag = jm*mag
					mag = timeStep / (dSquared * distance)
					distance = sqrt (dSquared)
					dSquared = dx*dx + dy*dy + dz*dz + eps
					dx = ix - jx
					dy = iy - jy
					dz = iz - jz
					Body ix iy iz _ _ _ _  = b_i
					Body jx jy jz _ _ _ jm = b_j
