// Time-stamp: <Sat May 26 03:16:32 BST 2012 pt114>
import scala.actors.Futures
import scala.io.Source
import scala.util.Random
import scala.annotation.tailrec

object allpairs {
	class Body (val mass: Double, val x: Double, val y: Double, val z: Double, val vx: Double, val vy: Double, val vz: Double)
		
	val timeStep = 0.001
	val eps = 0.01

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
    
	def doSteps(s: Int, bs: List[Body]): List[Body] = s match {
		case 0 => bs
		case _ => {
		
			def updateVel(b: Body) : Body = {
				/* //45.09s
				bs.map((b2: Body) => accel(b,b2)).foldLeft(b) ((acc,a) =>
					new Body (acc.mass,acc.x,acc.y,acc.z,(acc.vx-a._1),(acc.vy-a._2),(acc.vz-a._3))
				)*/
				//fold/map merge
				bs.foldLeft(b) ((state,b2) => {
					val (ax,ay,az) = accel(b,b2)
					new Body (state.mass,state.x,state.y,state.z,(state.vx-ax),(state.vy-ay),(state.vz-az))
				})				
			}
			
			def accel(bi:Body,bj:Body):(Double,Double,Double) = {
				if (bi == bj) (0.0,0.0,0.0)
				else {
					val dx = bi.x - bj.x
					val dy = bi.y - bj.y
					val dz = bi.z - bj.z
					val dSquared = dx*dx + dy*dy + dz*dz + eps					
					val distance = math.sqrt (dSquared)
					val mag = timeStep / (dSquared *  distance)
					val mmag = bj.mass*mag
					((dx*mmag),(dy*mmag),(dz*mmag))
				}
			}
			
			def updatePos(b: Body) : Body =
				new Body(b.mass,(b.x + timeStep * b.vx) , (b.y + timeStep * b.vy), (b.z + timeStep * b.vz),b.vx,b.vy,b.vz)
			
			val new_bs = bs.map(updateVel _ andThen updatePos)
			doSteps(s - 1, new_bs)
		}
	}
}
