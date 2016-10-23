{-# LANGUAGE CPP, BangPatterns, NoMonomorphismRestriction #-}

-- File: bh.hs
-- Description: nbody computation: barnes-hut algoritms
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

instance NFData Body
instance NFData Accel
instance NFData Centroid
instance NFData Bbox
instance NFData BHTree
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
	deriving Show

-- centroid
data Centroid
	= Centroid
	{ cx :: {-# UNPACK #-} !Double
	, cy :: {-# UNPACK #-} !Double
	, cz :: {-# UNPACK #-} !Double
	, cm :: {-# UNPACK #-} !Double }

-- acceleration
data Accel
	= Accel
	{ ax :: {-# UNPACK #-} !Double
	, ay :: {-# UNPACK #-} !Double
	, az :: {-# UNPACK #-} !Double }

-- bounding box - the min and max pts in 3D space, making a region containing points
data Bbox
	= Bbox
	{ minx :: {-# UNPACK #-} !Double
	, miny :: {-# UNPACK #-} !Double
	, minz :: {-# UNPACK #-} !Double
	, maxx :: {-# UNPACK #-} !Double
	, maxy :: {-# UNPACK #-} !Double
	, maxz :: {-# UNPACK #-} !Double }
	deriving Show

-- the BH tree
-- node consists of size, center X, Y, Z, total mass, and subtrees
data BHTree
	= BHT
	{ size     :: {-# UNPACK #-} !Double
	, centerx  :: {-# UNPACK #-} !Double
	, centery  :: {-# UNPACK #-} !Double
	, centerz  :: {-# UNPACK #-} !Double
	, mass     :: {-# UNPACK #-} !Double
	, subTrees :: ![BHTree] }


timeStep = 0.001

-- If the distance between the points is smaller than this
eps = 0.01 	-- then ignore the forces between them.

-- threshold is the ratio of size of a region and the distance between
-- a body and the region centre of mass.
-- If s / d < threshold, then the internal node is sufficiently far away.
threshold = 0.25 -- threshold=0 degenerates to brute force

getBodies n generate =
	if generate then
		do
			let
				newBody seed = Body x y z vx vy vz m
					where [x,y,z,vx,vy,vz,m] = take 7 $ randomList seed
				
				randomList :: (Random a) => Int -> [a]
				randomList seed = randoms (mkStdGen seed)
				
			return $ map newBody [0..n-1]
	else
		do
			str <- readFile "../input.txt" -- line: m x y z vx vy vz
			let
				ls = lines str
				bs = [bodyFromLine l|l<-ls]
				bodyFromLine line = Body (val 1) (val 2) (val 3) (val 4) (val 5) (val 6) (val 0)
					where
						val i = read (ws!!i) :: Double
						ws = (words line)
			return bs

-- find the coordinates of the bounding box that contains the bodies
findBounds bs = foldl' f (Bbox 0 0 0 0 0 0) bs
	where
		f (Bbox minx miny minz maxx maxy maxz) (Body x y z _ _ _ _) =
			Bbox (min minx x) (min miny y) (min minz z) (max maxx x) (max maxy y) (max maxz z)

-- calculate the centroid of points
calcCentroid :: [Body] -> Centroid
calcCentroid bs = Centroid (xs/mass) (ys/mass) (zs/mass) mass
  where
    Centroid xs ys zs mass = foldl' f (Centroid 0 0 0 0) bs
    f (Centroid x y z m) (Body x' y' z' _ _ _ m') = Centroid (m'*x'+x) (m'*y'+y) (m'*z'+z) (m'+m)

-- check if point is in box
inBox :: Bbox -> Body -> Bool
inBox (Bbox minx miny minz maxx maxy maxz) (Body x y z _ _ _ _)
 	= (x > minx) && (x <= maxx) && (y > miny) && (y <= maxy) && (z > minz) && (z <= maxz)

-- split points according to their locations in the box
splitPoints :: Bbox -> [Body] -> [(Bbox, [Body])]
splitPoints bb []  	= [(bb,[])]
splitPoints bb [b] = [(bb,[b])]
splitPoints bb bs  = boxesAndPts
	where	
		Bbox minx miny minz maxx maxy maxz = bb
		p1		= [b|b<-bs,inBox box1 b]
		p2		= [b|b<-bs,inBox box2 b]
		p3		= [b|b<-bs,inBox box3 b]
		p4		= [b|b<-bs,inBox box4 b]
		p5		= [b|b<-bs,inBox box5 b]
		p6		= [b|b<-bs,inBox box6 b]
		p7		= [b|b<-bs,inBox box7 b]
		p8		= [b|b<-bs,inBox box8 b]
		box1		= Bbox minx miny minz midx midy midz
		box2		= Bbox minx midy minz midx maxy midz
		box3		= Bbox midx miny minz maxx midy midz
		box4		= Bbox midx midy minz maxx maxy midz
		box5		= Bbox minx miny midz midx midy maxz
		box6		= Bbox minx midy midz midx maxy maxz
		box7		= Bbox midx miny midz maxx midy maxz
		box8		= Bbox midx midy midz maxx maxy maxz
		boxes		= box1:box2:box3:box4:box5:box6:box7:box8:[]
		splitPts	= p1:p2:p3:p4:p5:p6:p7:p8:[]
		boxesAndPts = [ (box,pts) | (box,pts) <- zip boxes splitPts, not (null pts) ]
		(midx, midy, midz)	= ((minx + maxx) / 2.0 , (miny + maxy) / 2.0 , (minz + maxz) / 2.0) 

-- build the Barnes-Hut tree
buildTree :: Bbox -> [Body] -> BHTree
buildTree bb bs 
 | length bs <= 1	= BHT size centerx centery centerz mass []
 | otherwise		= BHT size centerx centery centerz mass subTrees
 where
 	Centroid centerx centery centerz mass = calcCentroid bs
	boxesAndPts		= splitPoints bb bs
	subTrees		= map (\(bb',bs') -> buildTree bb' bs') boxesAndPts
	size = minimum [abs (maxx - minx), abs (maxy - miny), abs (maxz - minz)]
	Bbox minx miny minz maxx maxy maxz = bb

-- calculate the acceleration on a point due to some other point
accel :: BHTree -> Body -> Accel
accel (BHT _ centerx centery centerz mass _) (Body x y z _ _ _ m) = Accel (dx * mmag ) (dy * mmag) (dz * mmag)
	where
		dsqr = (dx * dx) + (dy * dy) + (dz * dz) + eps
		d    = sqrt dsqr 
		dx   = x - centerx
		dy   = y - centery
		dz   = z - centerz
		mag = timeStep / (dsqr * d)
		mmag = mass * mag

-- use centroid as approximation if the point is far from a cell
isFar :: BHTree -> Body -> Bool
isFar (BHT size centerx centery centerz _ _) (Body x y z _ _ _ _)
 = let	dx	= centerx - x
	dy	= centery - y
	dz	= centerz - z
	dist	= sqrt (dx * dx + dy * dy + dz * dz)
   in	(size / dist) < threshold

-- calculate the accelleration of a point due to the points in the given tree

--29.05.2012
--optimisation used: foldr/build used in conjunction with list comprehension - because laziness
calcAccel :: BHTree -> Body -> Accel
calcAccel t@(BHT _ _ _ _ _ subtrees) b
	| null subtrees	= accel t b
	| isFar t b		= accel t b
	| otherwise		= foldr addAccel (Accel 0 0 0) [calcAccel st b | st <- subtrees] -- fold is merged with the list comprehension (results in increased no. of converted sparks)

--combine fold with list manually
calcAccel1 :: BHTree -> Body -> Accel
calcAccel1 t@(BHT _ _ _ _ _ subtrees) b
	| null subtrees	= accel t b
	| isFar t b		= accel t b
	| otherwise		= foldl (\acc1 st -> addAccel acc1 (calcAccel1 st b)) (Accel 0 0 0) subtrees

--tail optimised as in scala/f# but does not perform better than foldr/build above
calcAccel2 :: BHTree -> Body -> Accel
calcAccel2 t b = calcAccelRec [t] (Accel 0 0 0)
	where
		calcAccelRec :: [BHTree] -> Accel -> Accel
		calcAccelRec [] acc = acc
		calcAccelRec ((t1@(BHT _ _ _ _ _ t1sub)):ts) acc = 
		    if (null t1sub || isFar t1 b) then calcAccelRec ts (addAccel acc (accel t1 b))
		    else calcAccelRec (t1sub ++ ts) acc

addAccel (Accel ax1 ay1 az1) (Accel ax2 ay2 az2) = Accel (ax1+ax2) (ay1+ay2) (az1+az2)

f acc (Body x y z vx vy vz m) = acc+x+y+z+vx+vy+vz+m

doSteps 0 bs = bs
doSteps s bs = doSteps (s-1) new_bs
	where						
		box = findBounds bs
		tree = buildTree box bs
#if defined (PAR)
		chunksize = (length bs) `quot` (numCapabilities * 4)
		new_bs = map f bs `using` parListChunk chunksize rdeepseq
		--new_bs = parMap rdeepseq f bs
#else
		new_bs = map f bs
#endif
		f = (updatePos . updateVel)
		updatePos (Body x y z vx vy vz m) = Body (x + timeStep * vx) (y + timeStep * vy) (z + timeStep * vz) vx vy vz m

		updateVel b@(Body x y z vx vy vz m) = Body x y z (vx-ax) (vy-ay) (vz-az) m
			where
				Accel ax ay az = calcAccel tree b

-- func to calc time taken between t0 and t1 in sec
secDiff::ClockTime -> ClockTime -> Float
secDiff (TOD secs1 psecs1) (TOD secs2 psecs2) = fromInteger (psecs2 - psecs1) / 1e12 + fromInteger (secs2 - secs1)

main =
	do
		args <- getArgs
		let 
			n = read (args!!0) :: Int -- num of bodies
			steps = read (args!!1) :: Int -- num of steps (iterations)
		bs <- getBodies n True
		let res1 = foldl' f 0 bs -- do meaningless sum...
		print res1
		t1 <- getClockTime
		let bs' = doSteps steps bs
		print $ foldl' f 0 bs'
		t2 <- getClockTime
		putStrLn $ printf "time taken: %.2fs" $ secDiff t1 t2
