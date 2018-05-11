// Implementation of the 3D fast marching method to solve the Eikonal equation
// Purpose is now to compute extension velocities and/or reinitialize a level
// set function

#include "stdafx.h"
#include "Interface3D.h"

using namespace std;
//const double INFTY = std::numeric_limits<double>::infinity();

// Requires the priority queue be filled with the initial conditions
// Works with arrays workGrid (for signed distance function) and velocities
void Interface3D::doFastMarch3D(Grid3<double> & T, Grid3<double> & V, EikonalHeap3D & initialConditions, bool reinitialize, double bandRadius) {

   EikonalHeap3D & priorityQueue = initialConditions;

   int nbSize = 0;

   //cerr << "Priority queue size: " << priorityQueue.currentSize() << endl;
   while (!priorityQueue.isEmpty()) {
      Node3D next = priorityQueue.remove();
	  int fineIndex = next.Index();

      //cerr << "*** Next = " << i << " " << j << " " << k << endl

      // Outdated: easier to do this than to remove all instances from the queue, etc.
      //if ((T(fineIndex) != INFTY && !(reinitialize && !narrowBand(fineIndex))) || next.getU() < -1.0) {
	  //  discarded++;
	  //  continue;
	  //}

	  if (next.getU() > bandRadius) { break; } // Stopping condition

	  if (reinitialize) {
		  narrowBand(Material, fineIndex) = true;
		  narrowBandIterator[Material][nbSize++] = next;
          // Put this outside the reinitialize part? Or no
	      signs(fineIndex) = interfaceGrid(Material, fineIndex) > 0.0 ? 1.0 : -1.0;
	  }
      
      T(fineIndex) = next.getU();
      V(fineIndex) = next.getV();
      
      // TODO: Check narrow band
      // Calculate T for all 6 neighbors and insert
	  // If we are reinitializing, want either that T is unassigned, which can be either because it is infinity or not in the band
	  // If we are not reinitializing, we need it to be both in the band and unassigned

	  int leftNeighbor = fineIndex + T.left,
		  rightNeighbor = fineIndex + T.right,
		  topNeighbor = fineIndex + T.top,
		  bottomNeighbor = fineIndex + T.bottom,
		  forwardNeighbor = fineIndex + T.forward,
		  backwardNeighbor = fineIndex + T.backward;

	  int i, j, k;
	  fineFromIndex(fineIndex, i, j, k);

      if (i > 0 && (reinitialize ? (T(leftNeighbor) == INFTY || !narrowBand(Material, leftNeighbor)) 
	                             : (T(leftNeighbor) == INFTY && narrowBand(Material, leftNeighbor)))) {
		  Node3D possible = TrialNode(T, V, leftNeighbor, i-1, j, k);
		  if (possible.getU() <= bandRadius) {
              priorityQueue.insert( possible );
		  }
      }
      if (i < m && (reinitialize ? (T(rightNeighbor) == INFTY || !narrowBand(Material, rightNeighbor)) 
	                             : (T(rightNeighbor) == INFTY && narrowBand(Material, rightNeighbor)))) {
		  Node3D possible = TrialNode(T, V, rightNeighbor, i+1, j, k);
		  if (possible.getU() <= bandRadius) {
              priorityQueue.insert( possible );
		  }          
	  }
      if (j > 0 && (reinitialize ? (T(bottomNeighbor) == INFTY || !narrowBand(Material, bottomNeighbor)) 
	                             : (T(bottomNeighbor) == INFTY && narrowBand(Material, bottomNeighbor)))) {
          Node3D possible = TrialNode(T, V, bottomNeighbor, i, j-1, k);
		  if (possible.getU() <= bandRadius) {
              priorityQueue.insert( possible );
		  }
	  }
      if (j < n && (reinitialize ? (T(topNeighbor) == INFTY || !narrowBand(Material, topNeighbor)) 
	                             : (T(topNeighbor) == INFTY && narrowBand(Material, topNeighbor)))) {
          Node3D possible = TrialNode(T, V, topNeighbor, i, j+1, k);
		  if (possible.getU() <= bandRadius) {
              priorityQueue.insert( possible );
		  }
      }
      if (k > 0 && (reinitialize ? (T(backwardNeighbor) == INFTY || !narrowBand(Material, backwardNeighbor)) 
	                             : (T(backwardNeighbor) == INFTY && narrowBand(Material, backwardNeighbor)))) {
          Node3D possible = TrialNode(T, V, backwardNeighbor, i, j, k-1);
		  if (possible.getU() <= bandRadius) {
              priorityQueue.insert( possible );
		  }
      }
      if (k < p && (reinitialize ? (T(forwardNeighbor) == INFTY || !narrowBand(Material, forwardNeighbor)) 
	                             : (T(forwardNeighbor) == INFTY && narrowBand(Material, forwardNeighbor)))) {
          Node3D possible = TrialNode(T, V, forwardNeighbor, i, j, k+1);
		  if (possible.getU() <= bandRadius) {
              priorityQueue.insert( possible );
		  }
      }
   }

   if (reinitialize) {
	   narrowBandSize[Material] = nbSize;
   }
   // Done!
}

Node3D Interface3D::TrialNode(Grid3<double> & T, Grid3<double> & V, int center, int i, int j, int k) {
   // Solve a quadratic equation in Tijk.
   // Uses the scheme |Grad T|^2 = max(D-x T, -D+x T, 0)^2
   //                            + max(D-y T, -D+y T, 0)^2
   //                            + max(D-z T, -D+z T, 0)2^www
   // and causality, e.g. that T_ijk is bigger than all known surrounding values   
   double xNeighbor = -1.0, yNeighbor = -1.0, zNeighbor = -1.0;
   double xNeighborV = 0.0, yNeighborV = 0.0, zNeighborV = 0.0;
   double xDir = 0.0, yDir = 0.0, zDir = 0.0;

   const int left = center + T.left, right = center + T.right,
	         top = center + T.top, bottom = center + T.bottom,
			 forward = center + T.forward, backward = center + T.backward;

   // To use the value it must be both in the narrow band and assigned
   if (i > 0 && T(left) < INFTY && narrowBand(Material, left)) {
       xNeighbor = T(left);
       xNeighborV = V(left);
       xDir = 1.0;
   }

   if (i < m && T(right) < INFTY && narrowBand(Material, right)) {
      if (xNeighbor < -0.5 || T(right) < xNeighbor) {
         xNeighbor = T(right);
         xNeighborV = V(right);
         xDir = -1.0;
      }
   }
   
   if (j > 0 && T(bottom) < INFTY && narrowBand(Material, bottom)) {
      yNeighbor = T(bottom);
      yNeighborV = V(bottom);
      yDir = 1.0;
   }
   
   if (j < n && T(top) < INFTY && narrowBand(Material, top)) {
      if (yNeighbor < -0.5 || T(top) < yNeighbor) {
         yNeighbor = T(top);
         yNeighborV = V(top);
         yDir = -1.0;
      }
   }

   if (k > 0 && T(backward) < INFTY && narrowBand(Material, backward)) {
      zNeighbor = T(backward);
      zNeighborV = V(backward);
      zDir = 1.0;
   }
   
   if (k < p && T(forward) < INFTY && narrowBand(Material, forward)) {
      if (zNeighbor < -0.5 || T(forward) < zNeighbor) {
         zNeighbor = T(forward);
         zNeighborV = V(forward);
         zDir = -1.0;
      }
   }

   double xExists = xDir * xDir, yExists = yDir * yDir, zExists = zDir * zDir;
   // Solve (if x) (T_ijk - xNeighbor)^2 + (if y) (T_ijk - yNeighbor)^2 
   //     + (if z) (T_ijk - zNeighbor)^2 = h^2/c^2 
   //
   // = (if x) [T_ijk^2 - 2(T_ijk)xNeighbor + xNeighbor^2] + 
   //   (if y) [T_ijk^2 - 2(T_ijk)yNeighbor + yNeighbor^2] 
   //   (if z) [T_ijk^2 - 2(T_ijk)zNeighbor + zNeighbor^2] = h^2/c^2
   // so (if x + if y + if z) T_ijk^2 - 2(T_ijk)((if x) xNeighbor + (if y) yNeighbor + (if z) zNeighbor) 
   //  + (if x) xNeighbor^2 + (if y) yNeighbor^2 + (if z) zNeighbor^2 - h^2/c^2;
   double h = hx; // For now
   double a = xExists + yExists + zExists, 
          b = -2.0 * (xExists * xNeighbor + yExists * yNeighbor + zExists * zNeighbor), 
          c = xExists * xNeighbor*xNeighbor + yExists * yNeighbor*yNeighbor 
            + zExists * zNeighbor*zNeighbor - h*h;
   double disc = b * b - 4 * a * c;
   // + is only acceptable solution because of the causality assumptions we made.
   double u;
   if (disc < 0.0) {
       u = -b / (2.0 * a);
       cerr << "Disc = " << disc << " with xNeighbor =, yNeighbor = " << xNeighbor << ", " << yNeighbor << endl;
   } else {
       u = (-b + sqrt(disc)) / (2.0 * a);
   }
   
   // Calculate v.
   double gradx = xDir * (u - xNeighbor) / h,
          grady = yDir * (u - yNeighbor) / h,
          gradz = zDir * (u - zNeighbor) / h;
   
   double v = (gradx * xNeighborV * xDir + grady * yNeighborV * yDir + gradz * zNeighborV * zDir) /
              (gradx * xDir + grady * yDir + gradz * zDir);   
   
   if (v < 0 && iter > 24) {
// 	   cerr << "Negative v\n";
   }
//   cout << "XDir and YDir " << xDir << ", " << yDir << endl;
//   cout << "Processing (temp) node: " << i << ", " << j << "; v = " << v << endl;
//   cout << "Neighbors (u): " << xNeighbor << "; " << yNeighbor << endl;

   return Node3D(u, v, center);
}
