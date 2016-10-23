// Time-stamp: <Sat May 26 03:16:32 BST 2012 pt114>

import scala.actors.Futures
import scala.io.Source
import scala.util.Random
import scala.annotation.tailrec

object bh {
	class Body (val mass: Double, val x: Double, val y: Double, val z: Double, val vx: Double, val vy: Double, val vz: Double)
	
	class Bbox (val minx: Double, val miny: Double, val minz: Double,
						val maxx: Double, val maxy: Double, val maxz: Double)

	class BHTree (val size: Double, val mass: Double,
							val x: Double, val y: Double, val z: Double,
							val nodes: List[BHTree])
	
	val timeStep = 0.001
	val eps = 0.01
	val threshold = 0.25
	
	//  Builds the Barnes-Hut tree.
	//@tailrec
	def buildTree(bb: Bbox, bs: List[Body]): BHTree = {
		val size = List(math.abs(bb.maxx - bb.minx), math.abs(bb.maxy - bb.miny), math.abs(bb.maxz - bb.minz)).min
		val (mass,x,y,z) = calcCentroid(bs)
		
		bs.length match {
			case 0 | 1 =>
				new BHTree(size,mass,x,y,z,List())
			case _ =>
				new BHTree(size,mass,x,y,z,
					splitPoints(bb, bs).map(sp => buildTree(sp._1,sp._2)).toList)
		}
	}
	
	// Calculates the centroid of bodies.
	@inline def calcCentroid(bs: List[Body]): (Double, Double, Double, Double) = {
		val (mass,xs,ys,zs) =
			bs.foldLeft(0.0, 0.0, 0.0, 0.0) ((acc, b: Body) =>
				(acc._1 + b.mass, acc._2 + (b.x * b.mass),
				acc._3 + (b.y * b.mass), acc._4 + (b.z * b.mass)))
		
		(mass,xs/mass,ys/mass,zs/mass)
	}
	
	// Finds the coordinates of the bounding box that contains the bodies.
	@inline def findBounds(bs: List[Body]): Bbox = {
		val bb = bs.foldLeft(0.0, 0.0, 0.0, 0.0, 0.0, 0.0) ((acc, b: Body) =>
			(math.min(acc._1, b.x), math.min(acc._2, b.y), math.min(acc._3, b.z),
			math.max(acc._4, b.x), math.max(acc._5, b.y), math.max(acc._6, b.z)))
		
		new Bbox(bb._1, bb._2, bb._3, bb._4, bb._5, bb._6)
	}
	
	// Checks if the body is inside the given bounding box.
	@inline def inBox(bb: Bbox, b: Body): Boolean =
		(b.x > bb.minx) && (b.x <= bb.maxx) &&
		(b.y > bb.miny) && (b.y <= bb.maxy) &&
		(b.z > bb.minz) && (b.z <= bb.maxz)
	
	// Splits points according to their locations in the box.
	@inline def splitPoints(bb: Bbox, bs: List[Body]): Map[Bbox, List[Body]] = bs.length match {
		case 0 => Map((bb, List()))
		case 1 => Map((bb, bs))
		case _ =>
			val (minx, miny, minz) = (bb.minx, bb.miny, bb.minz)
			val (maxx, maxy, maxz) = (bb.maxx, bb.maxy, bb.maxz)
			val (midx, midy, midz) = ((minx + maxx) / 2.0, (miny + maxy) / 2.0, (minz + maxz) / 2.0)

			val boxes =
				List(new Bbox(minx, miny, minz, midx, midy, midz),
					new Bbox(minx, midy, minz, midx, maxy, midz),
					new Bbox(midx, miny, minz, maxx, midy, midz),
					new Bbox(midx, midy, minz, maxx, maxy, midz),
					new Bbox(minx, miny, midz, midx, midy, maxz),
					new Bbox(minx, midy, midz, midx, maxy, maxz),
					new Bbox(midx, miny, midz, maxx, midy, maxz),
					new Bbox(midx, midy, midz, maxx, maxy, maxz))
			
			bs.groupBy {
				case b if inBox(boxes(0), b) => boxes(0)
				case b if inBox(boxes(1), b) => boxes(1)
				case b if inBox(boxes(2), b) => boxes(2)
				case b if inBox(boxes(3), b) => boxes(3)
				case b if inBox(boxes(4), b) => boxes(4)
				case b if inBox(boxes(5), b) => boxes(5)
				case b if inBox(boxes(6), b) => boxes(6)
				case _ => boxes(7)
			}
	}
	
	@inline def addAccel(acc1:(Double, Double, Double),acc2:(Double, Double, Double)): (Double, Double, Double) = ((acc1._1+acc2._1),(acc1._2+acc2._2),(acc1._3+acc2._3))
	
	@inline def isFar(t: BHTree, b: Body) = {
		val dx = t.x - b.x
		val dy = t.y - b.y
		val dz = t.z - b.z
		val dist = math.sqrt(dx * dx + dy * dy + dz * dz)
		
		((t.size / dist) < threshold)
	}
	
	@inline def accel(t: BHTree, b: Body) = {
		val dx   = b.x - t.x
		val dy   = b.y - t.y
		val dz   = b.z - t.z
		val dsqr = (dx * dx) + (dy * dy) + (dz * dz) + eps
		val d = math.sqrt(dsqr)
		val mag = timeStep / (dsqr * d)
		val mmag = t.mass * mag
		((dx * mmag),(dy * mmag),(dz * mmag))
	}
	
	//original calcAccel (translated from haskell/f#)
	def calcAccel(t: BHTree, b: Body): (Double, Double, Double) = t.nodes match {
		case Nil => accel(t,b)
		case _ => 
			if (isFar(t,b)) accel(t,b)
			else  {
            t.nodes.map(st => 
            	calcAccel(st,b)).foldLeft(0.0,0.0,0.0)(addAccel)
        }
	}
	
	// this version merges the fold and map
	def calcAccel1(t: BHTree, b: Body): (Double, Double, Double) = t.nodes match {
		case Nil => accel(t,b)
		case _ => 
			if (isFar(t,b)) accel(t,b)
			else  t.nodes.foldLeft(0.0,0.0,0.0) ((acc,st) => addAccel(acc,calcAccel1(st,b)))
	}
	
	def calcAccel11(t: BHTree, b: Body): (Double, Double, Double) = {
	
		val dx   = b.x - t.x
		val dy   = b.y - t.y
		val dz   = b.z - t.z
		
		val dsqr = math.sqrt(dx * dx + dy * dy + dz * dz)
	
		t.nodes match {
			case Nil => 
				val d = math.sqrt(dsqr+eps)
				val mag = timeStep / ((dsqr+eps) * d)
				val mmag = t.mass * mag
				((dx * mmag),(dy * mmag),(dz * mmag))
			case _ => 			
				if ((t.size / dsqr) < threshold) {
					val d = math.sqrt(dsqr)
					val mag = timeStep / (dsqr * d)
					val mmag = t.mass * mag
					((dx * mmag),(dy * mmag),(dz * mmag))
				}
				else  {
					t.nodes.foldLeft(0.0,0.0,0.0) ((acc1,st) => {
						val acc2 = calcAccel11(st,b)
						((acc1._1+acc2._1),(acc1._2+acc2._2),(acc1._3+acc2._3))
					})
				}
			}
	}
	
	// calcAccelrec
	def calcAccel2(t: BHTree, b: Body): (Double, Double, Double) = {
		
		def calcAccelRec(subtrees:List[BHTree],acc:(Double, Double, Double)): (Double, Double, Double) = subtrees match {
			case Nil => acc
			case t1::ts => {
				if (t1.nodes.isEmpty || isFar(t1,b)) calcAccelRec(ts,addAccel(acc,accel(t1,b)))
				else calcAccelRec((t1.nodes ++ ts),acc)
			}
		}
		calcAccelRec(List(t), (0,0,0))
	}       
	
	def doSteps(s: Int, bs: List[Body]): List[Body] = s match {
		case 0 => bs
		case _ => {
			val box = findBounds(bs)
			val tree = buildTree(box, bs)
			
			def updateVel(b: Body) : Body = {
				val (ax,ay,az) = calcAccel1(tree, b)

				new Body(b.mass,b.x,b.y,b.z,(b.vx - ax), (b.vy - ay), (b.vz - az))
			}
			
			def updatePos(b: Body) : Body =
				new Body(b.mass,(b.x + timeStep * b.vx) , (b.y + timeStep * b.vy), (b.z + timeStep * b.vz),b.vx,b.vy,b.vz)
			
			val new_bs = bs.map(updateVel _ andThen updatePos)
			
			doSteps(s-1, new_bs)
		}
	}
	
	// Returns nbodies either from generator or input file.
	def getBodies(num: Int, generate: Boolean): List[Body] = {		
		if (generate) {
			def newBody(seed: Int) = {
				val rnd = new Random(seed)
				new Body(rnd.nextDouble(),rnd.nextDouble(), rnd.nextDouble(), rnd.nextDouble(),rnd.nextDouble(), rnd.nextDouble(), rnd.nextDouble())
			}
			
			for (i <- List.range(0, (num-1))) yield (newBody(i))
		}
		else {
			Source.fromFile("input.txt").getLines().toList
				.map((x: String) => x.split(" ").map((y: String) => y.toDouble))
				.map((arr: Array[Double]) => new Body(arr(0),arr(1), arr(2), arr(3),arr(4), arr(5), arr(6)))
		}
	}

	
	def f(bs: List[Body]): Double =
		bs.foldLeft(0.0) ((acc: Double, b: Body) =>
			acc + b.x + b.y + b.z + b.vx + b.vy + b.vz + b.mass)
			
	def main(args: Array[String]) = {
		val n = args(0).toInt
		val steps = args(1).toInt
		val bs = getBodies(n, true)		
		println(f(bs))		
		val startTime = System.currentTimeMillis
		val new_bs = doSteps(steps, bs)
		println(f(new_bs))
		val endTime = System.currentTimeMillis-startTime
		println("\ntime taken: (" + endTime / 1000.0 + " s)")
    }
}
