#include <math.h>
#include <iostream>
#include <fstream>
#include "InterfaceFromVolumeFraction.h"
#include "f2c.h"
#include "clapack.h"
#include <sstream>

using namespace std;

const int Interface3D::oppositeEdge[6] = {5, 4, 3, 2, 1, 0};
const double Interface3D::vertices[4][3] = { {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };
const int Interface3D::edgesTouchingVertices[4][3] = {{0,1,2}, {0,3,4}, {1,3,5}, {2,4,5}};

// Updated.
Interface3D::Interface3D(int gm, int gn, int gp, int _ns, 
    double _h, ifstream & input): 
    coarseM(gm), m(gm * _ns), coarseN(gn), n(gn * _ns), coarseP(gp), p(gp * _ns), 
    coarseSize((gm+2) * (gn+2) * (gp+2)), fineSize((m+1) * (n+1) * (p+1)), ns(_ns), hx(_h),
    hy(_h), hz(_h),
	grid(coarseM+2, coarseN+2, coarseP+2),
    vActs(coarseM+2, coarseN+2, coarseP+2), 
    surfaceAreas(coarseM+2, coarseN+2, coarseP+2), 
    coarseNarrowBand(coarseM+2, coarseN+2, coarseP+2), 
    vActsBack(coarseM+2, coarseN+2, coarseP+2), 
    lambdas(coarseM+2, coarseN+2, coarseP+2), 
    residual(coarseM+2, coarseN+2, coarseP+2), 
    gradientBand(coarseM+2, coarseN+2, coarseP+2), 
    // Fine initializations
    interfaceGrid(m+1, n+1, p+1),
    velocities(m+1, n+1, p+1),
    workGrid(m+1, n+1, p+1),
    oldInterface(m+1, n+1, p+1),
    signs(m+1, n+1, p+1),
    narrowBand(m+1, n+1, p+1), 
    volumeGradients(m+1, n+1, p+1),
	volumeOfGradients(coarseM+2, coarseN+2, coarseP+2)
{       
       allocate();
            
       for (int K = 1; K < coarseP+1; K++) {
           for (int I = 1; I < coarseM+1; I++) {               
               for (int J = 1; J < coarseN+1; J++) {
                   input >> grid(I, J, K); // Read in volume fractions
               }
           }
       }
}

Interface3D::Interface3D(int gm, int gn, int gp, int _ns):
    coarseM(gm), m(gm * _ns), coarseN(gn), n(gn * _ns), coarseP(gp), p(gp * _ns), ns(_ns),
    coarseSize((gm+2) * (gn+2) * (gp+2)), fineSize((gm*_ns+1) * (gn*_ns+1) * (gp*_ns+1)),
    grid(gm+1, gn+1, gp+1),
    vActs(gm+1, gn+1, gp+1), 
    surfaceAreas(gm+1, gn+1, gp+1), 
    coarseNarrowBand(gm+1, gn+1, gp+1), 
    vActsBack(gm+1, gn+1, gp+1), 
    lambdas(gm+1, gn+1, gp+1), 
    residual(gm+1, gn+1, gp+1), 
    gradientBand(gm+1, gn+1, gp+1), 
    // Fine initializations
    interfaceGrid(gm*_ns, gn*_ns, gp*_ns),
    velocities(gm*_ns, gn*_ns, gp*_ns),
    workGrid(gm*_ns, gn*_ns, gp*_ns),
    oldInterface(gm*_ns, gn*_ns, gp*_ns),
    signs(gm*_ns, gn*_ns, gp*_ns),
    narrowBand(gm*_ns, gn*_ns, gp*_ns), 
    volumeGradients(gm*_ns, gn*_ns, gp*_ns),
	volumeOfGradients(gm+1, gn+1, gp+1),
	hx(1.0/_ns), hy(1.0/_ns), hz(1.0/_ns) {
    allocate();
}

// Allocates all data, given that we know how big to make this thing.
void Interface3D::allocate()
{             
       // Initialize volume fraction array
       // Empty front & back face
       for (int i = 0; i < coarseM+2; i++) {
          for (int j = 0; j < coarseN+2; j++) {
              grid(i, j, 0) = 0.0;
              grid(i, j, coarseP+1) = 0.0;
              
              vActs(i, j, 0) = 0.0; 
              vActs(i, j, coarseP+1) = 0.0;
          } 
       }             
            
       for (int k = 1; k < coarseP+1; k++) {
           for (int j = 0; j < coarseN+2; j++) {
               grid(0, j, k) = 0.0; // empty left
               grid(coarseM+1, j, k) = 0.0; // empty right
               
               vActs(0, j, k) = 0.0; // empty left
               vActs(coarseM+1, j, k) = 0.0; // empty right               
           }
           for (int i = 1; i < coarseM+1; i++) { 
               grid(i, 0, k) = 0.0; // Empty top
               vActs(i, 0, k) = 0.0;
               
               grid(i, coarseN+1, k) = 0.0; // Empty bottom
               vActs(i, coarseN+1, k) = 0.0;
           }
       }

       // Seems reasonable to just allocate to max capacity (but only use a small part of it)
       narrowBandIterator = new Node3D[(m+1)*(n+1)*(p+1)];
       narrowBandSize = 0;

	   for (Index3 ind(m, n, p); ind.valid(); ++ind) {
           volumeGradients(ind).setXYZ(ind.i, ind.j, ind.k);
       }

       iter = 0;
}
// Allocate

// Updates the listing of which coarse grid cells are cared about
void Interface3D::updateCoarseNarrowBand() {

	for (Index3 inds(1, 1, 1, coarseM, coarseN, coarseP, coarseM+1, coarseN+1, coarseP+1); inds.valid(); ++inds) {
        // Check neighbors
        bool closeToInterface = false;
        double vrCenter = grid(inds);

		for (Index3 neighbors(inds.i - 1, inds.j - 1, inds.k - 1, inds.i + 1, inds.j + 1, inds.k + 1, coarseM+1, coarseN+1, coarseP+1); neighbors.valid(); ++neighbors) {
            double vr = grid(neighbors), va = vActs(neighbors);
            if ((vr > 1e-8 && vr < 1 - 1e-8) || (va > 1e-8 && va < 1 - 1e-8) ||
                fabs(vr - vrCenter) > 1e-8) {
                closeToInterface = true;
                break;
            }
		}

        coarseNarrowBand(inds) = closeToInterface;
     }         
}

// Old seeding routine. Does a level set of the VoF data
void Interface3D::seedGrid(double levelSet) {

    for (int i = 0; i < fineSize; i++) {
        narrowBand(i) = false;
    }

    int index = 0;
    // Let's seed the interface with the given levelSet.   
    for (Index3 inds(m, n, p); inds.valid(); ++inds) {
        int i = inds.i, j = inds.j, k = inds.k;

        int imod = i % ns, jmod = j % ns, kmod = k % ns, xDir, yDir, zDir;
        if (imod >= ns / 2) { xDir = 1; } else { xDir = -1; }
        if (jmod >= ns / 2) { yDir = 1; } else { yDir = -1; }
        if (kmod >= ns / 2) { zDir = 1; } else { zDir = -1; }
                                  
        double xPos = (double) (imod) / ns - 0.5, yPos = (double) (jmod) / ns - 0.5,
               zPos = (double) (kmod) / ns - 0.5;
               xPos *= xDir; yPos *= yDir; zPos *= zDir;

        int I = (i/ns)+1, J = (j/ns)+1, K = (k/ns)+1;
        int iNext = I+xDir, jNext = J+yDir, kNext = K+zDir, i2Next = I+2*xDir, j2Next = J+2*yDir, k2Next = K+2*zDir, 
            iPrev = I-xDir, jPrev = J-yDir, kPrev = K-zDir;

        if (iNext > coarseM) { iNext = coarseM; }  if (i2Next > coarseM) { i2Next = coarseM; }
        if (jNext > coarseN) { jNext = coarseN; }  if (j2Next > coarseN) { j2Next = coarseN; }
        if (kNext > coarseP) { kNext = coarseP; }  if (k2Next > coarseP) { k2Next = coarseP; }
        if (iNext < 0) { iNext = 0; }  if (i2Next < 0) { i2Next = 0; }
        if (jNext < 0) { jNext = 0; }  if (j2Next < 0) { j2Next = 0; }
        if (kNext < 0) { kNext = 0; }  if (k2Next < 0) { k2Next = 0; }

        if (iPrev < 0) { iPrev = 0; }  if (iPrev > coarseM) { iPrev = coarseM; }
        if (jPrev < 0) { jPrev = 0; }  if (jPrev > coarseN) { jPrev = coarseN; }
        if (kPrev < 0) { kPrev = 0; }  if (kPrev > coarseP) { kPrev = coarseP; }       

        double stencil[64];
        int xcoords[4], ycoords[4], zcoords[4];
        xcoords[0] = iPrev; xcoords[1] = I; xcoords[2] = iNext; xcoords[3] = i2Next;
        ycoords[0] = jPrev; ycoords[1] = J; ycoords[2] = jNext; ycoords[3] = j2Next;
        zcoords[0] = kPrev; zcoords[1] = K; zcoords[2] = kNext; zcoords[3] = k2Next;
		                 
		for (int c = 0; c < 64; c++) {
            int xc = c % 4, yc = (c / 4) % 4, zc = c / 16;
            int II = xcoords[xc], JJ = ycoords[yc], KK = zcoords[zc];
            stencil[c] = grid(II, JJ, KK);
        }
                
        double volfrac = tricubicInterpolation(xPos, yPos, zPos, stencil);

        interfaceGrid(inds) = levelSet - volfrac;
        signs(inds) = (levelSet - volfrac) > 0.0 ? 1.0 : -1.0;

        // Init narrow band
        const int scansize = 3;
        if (i < scansize || j < scansize || k < scansize) { continue; }

        int nMinus = 0, nPlus = 0;
		for (Index3 iprimes(i - scansize, j - scansize, k - scansize, i, j, k, m, n, p); iprimes.valid(); ++iprimes) {
            if (interfaceGrid(iprimes) <= 0.0) {
                nMinus++;
            } else {
				nPlus++;
			}
		}

        if (nMinus > 0 && nPlus > 0) {
            for (Index3 iprimes(i - scansize, j - scansize, k - scansize, i, j, k, m, n, p); iprimes.valid(); ++iprimes) {
                if (!narrowBand(iprimes)) {
                    narrowBandIterator[narrowBandSize++] = Node3D(0.0, 0.0, iprimes);
                    narrowBand(iprimes) = true;
				}
			}
		}
	}
     
     computeVActs();
     updateCoarseNarrowBand();
     checkResult();
     std::cerr << "Narrow band size initial = " << narrowBandSize << endl;
     //reinitialize();
}

// Updated
double Interface3D::iterateLevelSetMethod(double epsilon, double timeStep, bool precise) {

    // Now, must choose the correct time step.
    //double timeStep = calculateBestTimeStep(epsilon, timeStep, currentMag);
    //cout << "Time step chosen was " << timeStep << " (k = " << k << ")\n";
      
 	for (int inds = 0; inds < coarseSize; inds++) {
 		if (coarseNarrowBand(inds)) {
			vActsBack(inds) = vActs(inds);
		}
 	}

    if (precise) {   
        computeCurvatureFlow(epsilon);
        updateInterface(timeStep);
       
        for (int volIter = 0; volIter < 1; volIter++) {
            double volError = computeVActs();

            std::cerr << "Vol iter " << volIter << "; error " << volError << endl;
            //checkResult();
            if (volError < 1e-4) { break; }

            // timeStep is not important in this place
            extensions(false, timeStep, true);
            computeVolumeFlow(0.0, timeStep, true);
            // This time step is suppose to be around 1. The fact that it isn't means the linearization of the
            // change in area isn't working well... Check again with new volGradients!
            updateInterface(1.0);
        }       
    } else {
        computeVActs();
        extensions(false, timeStep, false);
        computeVolumeFlow(epsilon, timeStep, false);
        updateInterface(timeStep);
    }       
                  
    double change = 0.0;
              
    computeVActs(); 
       // Check total change in volume
	for (int inds = 0; inds < coarseSize; inds++) {
        if (coarseNarrowBand(inds)) {
            change += fabs(vActsBack(inds) - vActs(inds));
		}
	}
       
    return change;
}

double Interface3D::computeVActs() {
     // Precompute all the volumes in each cell
     // Contains ghost layers.
     double changeInFraction = 0.0;
     double error = 0.0;
     totalSurfaceArea = 0.0;

	 for (Index3 inds(1,1,1, coarseM, coarseN, coarseP, coarseM+1, coarseN+1, coarseP+1); inds.valid(); ++inds) {
		 if (coarseNarrowBand(inds)) {
             double change = calculateVolumeImmersed(inds.i, inds.j, inds.k);
             changeInFraction += change;
             error += fabs(vActs(inds) - grid(inds));
         }
	 }

      // Ghost layer should be unchanged; remains zero.
      double vol = 0.0;

      updateCoarseNarrowBand();
      
      return error;//changeInFraction;
}

double Interface3D::sigmoid(double x, double L)
{
       return 1.0 / (1.0 + exp(-L*(x-0.5)));
}

void Interface3D::updateInterface(double timeStep) {
       for (int nbInd = 0; nbInd < narrowBandSize; nbInd++) {
           Node3D nbElt = narrowBandIterator[nbInd];
		   int ind = nbElt.Index();
           
           interfaceGrid(ind) -= timeStep * workGrid(ind);
       }
}

// Computes the delta to the next iteration and puts it into workGrid.
double Interface3D::computeVolumeFlow(double epsilon, double timeStep, bool precise) {
       //double changeInFraction = computeVActs();       
       iter++; // just for internal debugging.
              
       //double maxF = 0.0;
       //int maxFi, maxFj;
       double magnitude = 0.0;
             
       // For each node in the narrow band
       for (int nbInd = 0; nbInd < narrowBandSize; nbInd++) {
           Node3D nbElt = narrowBandIterator[nbInd];
		   int centerIndex = nbElt.Index();
                                                                   
           // Need neighbording points for derivatives
           double center, right, left, up, down, forward, back;
           double topright, topleft, botright, botleft;
           double forwardright, forwardleft, backright, backleft;
           double forwardtop, forwardbot, backtop, backbot;

           center = interfaceGrid(centerIndex);
           right = interfaceGrid(centerIndex + interfaceGrid.right);
           left = interfaceGrid(centerIndex + interfaceGrid.left);
           up = interfaceGrid(centerIndex + interfaceGrid.top);
           down = interfaceGrid(centerIndex + interfaceGrid.bottom);
		   forward = interfaceGrid(centerIndex + interfaceGrid.forward);
		   back = interfaceGrid(centerIndex + interfaceGrid.backward);
                              
		   topright = interfaceGrid(centerIndex + interfaceGrid.top + interfaceGrid.right);
		   topleft = interfaceGrid(centerIndex + interfaceGrid.top + interfaceGrid.left);
           botright = interfaceGrid(centerIndex + interfaceGrid.bottom + interfaceGrid.right);
           botleft = interfaceGrid(centerIndex + interfaceGrid.bottom + interfaceGrid.left);  
                             
           forwardright = interfaceGrid(centerIndex + interfaceGrid.forward + interfaceGrid.right);
           forwardleft = interfaceGrid(centerIndex + interfaceGrid.forward + interfaceGrid.left);
           backright = interfaceGrid(centerIndex + interfaceGrid.backward + interfaceGrid.right);
           backleft = interfaceGrid(centerIndex + interfaceGrid.backward + interfaceGrid.left);
                             
           forwardtop = interfaceGrid(centerIndex + interfaceGrid.top + interfaceGrid.forward);                           
           forwardbot = interfaceGrid(centerIndex + interfaceGrid.bottom + interfaceGrid.forward);    
           backtop = interfaceGrid(centerIndex + interfaceGrid.backward + interfaceGrid.top);      
           backbot = interfaceGrid(centerIndex + interfaceGrid.backward + interfaceGrid.bottom);
                      
           // Centered derivatives
           double uxx = (right - 2 * center + left)/(hx*hx),
                  uyy = (up - 2 * center + down)/(hy*hy),
                  uzz = (forward - 2 * center + back)/(hz*hz);
           double uxy = ((topright - topleft)/(2*hx) - (botright - botleft)/(2*hx)) / (2*hy),
                  uxz = ((forwardright - forwardleft)/(2*hx) - (backright - backleft)/(2*hx)) / (2*hz),
                  uyz = ((forwardtop - forwardbot)/(2*hy) - (backtop - backbot)/(2*hy)) / (2*hz);
           double ux  = (right - left)/(2*hx),
                  uy  = (up - down)/(2*hy),
                  uz  = (forward - back)/(2*hz);          
                      
           double gradLen = ux*ux + uy*uy + uz*uz, curvature;
           curvature = ( (uyy + uzz) * ux * ux + (uxx + uzz) * uy * uy + (uxx + uyy) * uz * uz
                       - 2 * ux * uy * uxy - 2 * ux * uz * uxz - 2 * uy * uz * uyz) 
                       / (gradLen * sqrt(gradLen) + hx*hz);
                                          
           //if (fabs(curvature) > 1.0 / hx) {
           //   std::cerr << "Too much curvature\n";
           //}
                              
           // Got to here.
                       
           double dxp = (right - center)/hx,
                  dxm = (center - left)/hx,
                  dyp = (up - center)/hy,
                  dym = (center - down)/hy,
                  dzp = (forward - center)/hz,
                  dzm = (center - back)/hz;
                              
           double gradPlus = 0.0, gradMinus = 0.0;
                       
           double gPlusScratch = max(dxm*fabs(dxm), 0.0) - min(dxp*fabs(dxp), 0.0)
                               + max(dym*fabs(dym), 0.0) - min(dyp*fabs(dyp), 0.0)
                               + max(dzm*fabs(dzm), 0.0) - min(dzp*fabs(dzp), 0.0);                                                 
                             
           if (gPlusScratch >= 0) { gradPlus = sqrt(gPlusScratch); }
           //else { std::cerr << gPlusScratch << endl; system("PAUSE");}
           double gMinusScratch = -min(dxm*fabs(dxm), 0.0) + max(dxp*fabs(dxp), 0.0)
                               - min(dym*fabs(dym), 0.0) + max(dyp*fabs(dyp), 0.0)
                               - min(dzm*fabs(dzm), 0.0) + max(dzp*fabs(dzp), 0.0);
           if (gMinusScratch >= 0) { gradMinus = sqrt(gMinusScratch); }
           //else { std::cerr << gMinusScratch << endl; system("PAUSE");}
                           
           double F = //smoothedSpeedFunc(i, j, k, timeStep);
                      velocities(centerIndex);
                             
		   double workValue = (gradPlus * max(F, 0.0) + gradMinus * min(F, 0.0)) 
                               - epsilon * curvature * sqrt(gradLen);
		   workGrid(centerIndex) = workValue;
                       
           magnitude += fabs(workValue);

       } // For narrow band iterator
       
       //magnitude = changeInFraction;
       return magnitude; // pop pop!
}


// Also puts it into workGrid.
double Interface3D::computeCurvatureFlow(double epsilon) {

       double magnitude = 0.0;

       // For each node in the narrow band
       for (int nbInd = 0; nbInd < narrowBandSize; nbInd++) {
            Node3D nbElt = narrowBandIterator[nbInd];
			int centerIndex = nbElt.Index();

            // Next step: implement curvature in 3D.               
            double right, left, up, down, forward, back;
            double topright, topleft, botright, botleft, center;
            double forwardright, forwardleft, backright, backleft;
            double forwardtop, forwardbot, backtop, backbot;

            center = interfaceGrid(centerIndex);
            right = interfaceGrid(centerIndex + interfaceGrid.right);
            left = interfaceGrid(centerIndex + interfaceGrid.left);
            up = interfaceGrid(centerIndex + interfaceGrid.top);
            down = interfaceGrid(centerIndex + interfaceGrid.bottom);
		    forward = interfaceGrid(centerIndex + interfaceGrid.forward);
		    back = interfaceGrid(centerIndex + interfaceGrid.backward);
                              
		    topright = interfaceGrid(centerIndex + interfaceGrid.top + interfaceGrid.right);
		    topleft = interfaceGrid(centerIndex + interfaceGrid.top + interfaceGrid.left);
            botright = interfaceGrid(centerIndex + interfaceGrid.bottom + interfaceGrid.right);
            botleft = interfaceGrid(centerIndex + interfaceGrid.bottom + interfaceGrid.left);  
                             
            forwardright = interfaceGrid(centerIndex + interfaceGrid.forward + interfaceGrid.right);
            forwardleft = interfaceGrid(centerIndex + interfaceGrid.forward + interfaceGrid.left);
            backright = interfaceGrid(centerIndex + interfaceGrid.backward + interfaceGrid.right);
            backleft = interfaceGrid(centerIndex + interfaceGrid.backward + interfaceGrid.left);
                             
            forwardtop = interfaceGrid(centerIndex + interfaceGrid.top + interfaceGrid.forward);                           
            forwardbot = interfaceGrid(centerIndex + interfaceGrid.bottom + interfaceGrid.forward);    
            backtop = interfaceGrid(centerIndex + interfaceGrid.backward + interfaceGrid.top);      
            backbot = interfaceGrid(centerIndex + interfaceGrid.backward + interfaceGrid.bottom);
                       
            double uxx = (right - 2 * center + left)/(hx*hx),
                    uyy = (up - 2 * center + down)/(hy*hy),
                    uzz = (forward - 2 * center + back)/(hz*hz);
            double uxy = ((topright - topleft)/(2*hx) - (botright - botleft)/(2*hx)) / (2*hy),
                    uxz = ((forwardright - forwardleft)/(2*hx) - (backright - backleft)/(2*hx)) / (2*hz),
                    uyz = ((forwardtop - forwardbot)/(2*hy) - (backtop - backbot)/(2*hy)) / (2*hz);
            double ux  = (right - left)/(2*hx),
                    uy  = (up - down)/(2*hy),
                    uz  = (forward - back)/(2*hz);              
                                        
                       
            double gradLen = ux*ux + uy*uy + uz*uz, curvature;
            curvature = ( (uyy + uzz) * ux * ux + (uxx + uzz) * uy * uy + (uxx + uyy) * uz * uz
                        - 2 * ux * uy * uxy - 2 * ux * uz * uxz - 2 * uy * uz * uyz) 
                        / (gradLen * sqrt(gradLen) + hx*hz);
                                          
            workGrid(centerIndex) = -epsilon * curvature * sqrt(gradLen);
       }
       
      return magnitude; // pop pop!
}

// Reinitializes the level set function to a signed distance function.
void Interface3D::extensions(bool reinitialize, double timeStep, bool precise) {
     
     EikonalHeap3D setup((m+1)*(n+1), m+1, n+1, p+1);
     
     // Find interface and add to queue.
     for (int nbInd = 0; nbInd < narrowBandSize; nbInd++) {
         Node3D nbElt = narrowBandIterator[nbInd];
		 int ind = nbElt.Index(), i, j, k;
		 fineFromIndex(ind, i, j, k);

         if (i == m || j == n || k == p) {
             continue; 
         }

		 // Skip if no intersection
		 int nMinus = 0, nPlus = 0;
		 for (Index3 neighbors(i, j, k, i+1, j+1, k+1, m, n, p); neighbors.valid(); ++neighbors) {
			 if (interfaceGrid(neighbors) > 0.0) { nPlus++; }
			 else { nMinus++; }
		 }

		 if (nMinus == 0 || nPlus == 0) { continue; }
     
         int indices[4][3], transform[3];
         // Each cube is divided into six tetrahedra
         // (0,0,0), (0,1,0), (0,1,1), (0,0,1)
         fillIndices(indices, i, j, k, i, j+1, k, i+1, j+1, k, i, j, k+1);
         // This tetrahedron is not rotated or flipped. (Identity transformation)
         transform[0] = 1;  transform[1] = 2; transform[2] = 3;
         interfaceDistanceSetup(indices, transform, timeStep, setup);
                 
         // (0,1,1), (0,1,0), (1,1,0), (0,0,1)
         fillIndices(indices, i, j+1, k+1, i, j+1, k, i+1, j+1, k, i, j, k+1);
         // x coord of this tetrahedron corresponds to real x coordinate
         // y coord corresponds to negative real z coordinate
         // z coord corresponds to negative real y coordinate
         transform[0] = 1;  transform[1] = -3; transform[2] = -2;
         interfaceDistanceSetup(indices, transform, timeStep, setup);
                  
         // (1,1,1), (0,1,1), (0,0,1), (1,1,0)
         fillIndices(indices, i+1, j+1, k+1, i, j+1, k+1, i, j, k+1, i+1, j+1, k);
         // x coordinate is really negative y
         // y coordinate is really negative x
         // z coordinate is really negative z
         transform[0] = -2;  transform[1] = -1; transform[2] = -3;
         interfaceDistanceSetup(indices, transform, timeStep, setup);
                 
         // (1,1,1), (1,0,1), (0,0,1), (1,1,0)
         fillIndices(indices, i+1, j+1, k+1, i+1, j, k+1, i, j, k+1, i+1, j+1, k);
         // x coordinate is really negative x
         // y coordinate is really negative y
         // z coordinate is really negative z
         transform[0] = -1;  transform[1] = -2;  transform[2] = -3;
         interfaceDistanceSetup(indices, transform, timeStep, setup);
                 
         // (1,0,0), (1,0,1), (0,0,1), (1,1,0)
         fillIndices(indices, i+1, j, k, i+1, j, k+1, i, j, k+1, i+1, j+1, k);
         // x is really -x
         // y is really +z
         // z is really +y
         transform[0] = -1;  transform[1] = 3;  transform[2] = 2;
         interfaceDistanceSetup(indices, transform, timeStep, setup);
                 
         // (0,0,0), (1,0,0), (1,1,0), (0,0,1)
         fillIndices(indices, i, j, k, i+1, j, k, i+1, j+1, k, i, j, k+1);            
         // x is really +y
         // y is really +x
         // z is really +z
         transform[0] = 2;  transform[1] = 1;  transform[2] = 3;
         interfaceDistanceSetup(indices, transform, timeStep, setup);
     } 

     // First, let's get the narrow banding to work on !precise mode.
     if (precise) {
         // setup now contains a heap of interface points; grab all items and throw onto our gradientElement list
         EikonalHeap3D workHeap(m*n, m+1, n+1, p+1);

         // Clear previous volume gradient info
         for (int nbInd = 0; nbInd < narrowBandSize; nbInd++) {
             Node3D nbElt = narrowBandIterator[nbInd];         
			 volumeGradients(nbElt.Index()).setU(-1.0);
         }					 
         allGradients.clear();

		 for (int ind = 0; ind < coarseSize; ind++) {
			 gradientBand(ind).clear();
		 }

         numGradElts = 0;
         while (!setup.isEmpty()) {
             Node3D next = setup.remove();
			 int ind = next.Index();
             //int i = next.X(), j = next.Y(), k = next.Z();
             //int I = 1 + (i / ns), J = 1 + (j / ns), K = 1 + (k / ns);
             //int iStart = (I-1)*ns, jStart = (J-1)*ns, kStart = (K-1)*ns;
			 gradientElement & ge = volumeGradients(ind); //volumeGradients(i,j,k);
             ge.setU(next.getU());
             ge.clearContributions();

        //	 //volumeGradients[I][J][K].push_back(ge);

             //volumeGradients[numGradElts++] = ge;
             //if (numGradElts == capacityGradElts) {
                // std::cerr << "Reallocating interface size\n";
                // // reallocate
                // gradientElement * resized = new gradientElement[capacityGradElts*2];
                // for (int cn = 0; cn < numGradElts; cn++) {
                //	 resized[cn] = volumeGradients[cn];
                // }
                // delete[] volumeGradients;
                // volumeGradients = resized;
                // capacityGradElts *= 2;
             //}

             allGradients.push_back(next);
         }

         // The very painful fmm gradient calculations
         for (std::vector<Node3D>::iterator it = allGradients.begin(); it != allGradients.end(); it++) {
             Node3D & next = *it;
             int ind = next.Index(), i, j, k;
			 fineFromIndex(ind, i, j, k);
             int I = (i/ns)+1, J = (j/ns)+1, K = (k/ns)+1;
             // EXPERIMENTING
             gradientElement & ge = volumeGradients(ind);

             updateGradients(ge, workHeap);
             // Let's try this silliness for now
             //int I2, J2, K2;
             //cellFromIndex(I2, J2, K2, next.getcell());
             //int Ip = I2 - I, Jp = J2 - J, Kp = K2 - K;
             //ge.addContribution(offsetFromDeltas(Ip, Jp, Kp), 1.0);

             // Try this for now -- just add all around you
             // This code is to put grad elements in buckets
             // Adding around is necessary for now to ensure that 0's are extended
             // on adjacent cells to interface too
             // POSSIBLY SLOW
             for (int offset = 0; offset < 27; offset++) {
                 double contribution = ge.getContribution(offset);
                 //if (contribution > 0.0 || contribution < 0.0) {
                     int Ip, Jp, Kp;
                     IndexDeltasFromOffset(Ip, Jp, Kp, offset);
                     gradientBand(I+Ip,J+Jp,K+Kp).push_back(&ge);
                 //}
             }
         }

         workHeap.discardAll();
         computeVolumeOfGradients(workHeap);
         // by this point it's already messed up
         computeLambdas();
         //gaussSeidelLambdas(3);

         for (std::vector<Node3D>::iterator it = allGradients.begin(); it != allGradients.end(); it++) {
             Node3D & next = *it;
			 int ind = next.Index(), i, j, k;
			 fineFromIndex(ind, i, j, k);
             int I = (i/ns)+1, J = (j/ns)+1, K = (k/ns)+1;

             double speed = 0.0;
             for (int offset = 0; offset < 27; offset++) {
                 double contribution = volumeGradients(ind).getContribution(offset);
                 if (contribution < 0.0 || contribution > 0.0) {
                     int Ip, Jp, Kp;				 			 
                     IndexDeltasFromOffset(Ip, Jp, Kp, offset);	
                     speed += lambdas(I+Ip, J+Jp, K+Kp) * contribution;
                 }
             }

             next.setV(speed);
             setup.insert(next);
         }
     } // precise

     // For each node in the narrow band
     for (int nbInd = 0; nbInd < narrowBandSize; nbInd++) {
         Node3D nbElt = narrowBandIterator[nbInd];
		 int ind = nbElt.Index();

         // Old comment. If there's a bug maybe it will offer insight.
         // Set some surrounding entries to infinity as well (Don't want FMM to check non narrow band elements)

         signs(ind) = interfaceGrid(ind) > 0.0 ? 1.0 : -1.0;
         workGrid(ind) = INFTY;
     }

	 if (reinitialize) {
		 for (int nbInd = 0; nbInd < narrowBandSize; nbInd++) {
			 narrowBand(narrowBandIterator[nbInd].Index()) = false;
		 }
     }
     
     // Issue: Assumes square grid for now
     // WARNING: h = hx only works for isotropic grid right now
     double radius = hx * ns * 2.0;
     doFastMarch3D(workGrid, velocities, setup, reinitialize, radius);

     if (!reinitialize) { return; }
     
     // Make signed distance function
     // For each node in the narrow band
     for (int nbInd = 0; nbInd < narrowBandSize; nbInd++) {
         Node3D nbElt = narrowBandIterator[nbInd];
		 int ind = nbElt.Index();
         
         double fmmResult = workGrid(ind);
         if (fmmResult < INFTY) { // in theory should not need this
             interfaceGrid(ind) = fmmResult * signs(ind);
         }					 
     }

     std::cerr << "Narrow band size after extensions " << narrowBandSize << endl;
}

void Interface3D::interfaceDistanceSetup(int indices[4][3], int transform[3], double timeStep, EikonalHeap3D & setup) {
    // The reason we keep track of indices is because I assume all points are
    // a-> (0,0,0), b-> (0,1,0), c-> (1,1,0), d->(0,0,1)
    // and to remember which point was the preimage, I supply the indices
    
    double a = interfaceGrid(indices[0][0], indices[0][1], indices[0][2]), 
           b = interfaceGrid(indices[1][0], indices[1][1], indices[1][2]), 
           c = interfaceGrid(indices[2][0], indices[2][1], indices[2][2]), 
           d = interfaceGrid(indices[3][0], indices[3][1], indices[3][2]);
           
    double nMinus = 0;
    if (a <= 0) { nMinus++; }
    if (b <= 0) { nMinus++; }
    if (c <= 0) { nMinus++; }
    if (d <= 0) { nMinus++; }
           
    if (nMinus == 0 || nMinus == 4) { return; }

    double v[6][3];
    bool edges[6];
           
    if (a * b  <= 0) { 
        v[0][0] = 0.0;
        v[0][1] = -a/(b-a);
        v[0][2] = 0.0;
        edges[0] = true;
    } else { edges[0] = false; }
           
    if (a * c <= 0) {
        v[1][0] = -a/(c-a);
        v[1][1] = v[1][0];
        v[1][2] = 0.0;
        edges[1] = true;
    } else { edges[1] = false; }
       
    if (a * d <= 0) {
        v[2][0] = 0.0;
        v[2][1] = 0.0;
        v[2][2] = -a/(d-a);
        edges[2] = true;
    } else { edges[2] = false; }
        
    if (b * c <= 0) {
        v[3][0] = -b/(c-b);
        v[3][1] = 1.0;
        v[3][2] = 0.0;
        edges[3] = true;
    } else { edges[3] = false; }
       
    if (b * d <= 0) {
        v[4][0] = 0.0;
        v[4][1] = d/(d-b);
        v[4][2] = -b/(d-b);
        edges[4] = true;
    } else { edges[4] = false; }
       
    if (c * d <= 0) {
        v[5][0] = d/(d-c);
        v[5][1] = v[5][0];
        v[5][2] = -c/(d-c);
        edges[5] = true;
    } else { edges[5] = false; }
    
    // Get the polygon that determines our interface.
    double base[4][3]; int e0;
    // Acquire starting vertex
    for (e0 = 0; e0 < 3; e0++) {
        if (edges[e0]) { 
            base[1][0] = v[e0][0];
            base[1][1] = v[e0][1];
            base[1][2] = v[e0][2];
            break;
        }
    }
    int j = 0;
    // Acquire adjacent vertices
    for (int e = e0+1; e < 6; e++) {
        if (edges[e] && oppositeEdge[e] != e0 && j < 2) {
            base[2*j][0] = v[e][0];
            base[2*j][1] = v[e][1];
            base[2*j][2] = v[e][2];
            j++;
        }
    }
    
    // Acquire final vertex, if applicable.
    int nBase = 3;
    if (edges[oppositeEdge[e0]]) {                 
        nBase = 4;
        base[3][0] = v[oppositeEdge[e0]][0];
        base[3][1] = v[oppositeEdge[e0]][1];
        base[3][2] = v[oppositeEdge[e0]][2];    
    }
    
    //std::cerr << "Base: ";
    //for (int r = 0; r < nBase; r++) {
    //    std::cerr << base[r][0] << " " << base[r][1] << " " << base[r][2] << "\n";
    //}
        
    double referenceCoord[3], finalCoord[3], velocity;
    // For each vertex of the tetrahedron
    for (int i = 0; i < 4; i++) {
        //std::cerr << "For vertex " << i << endl;
        // First, let us compute the distance to the plane f(x,y,z)=0
        // f(x,y,z) = a + (b-a)y + (c-b)x + (d-a)z
        double c0 = a, c1 = c-b, c2 = b-a, c3 = d-a;
        
        double n2 = c1*c1 + c2*c2 + c3*c3,
               t = -(c0 + c1 * vertices[i][0] + c2 * vertices[i][1] 
                   + c3 * vertices[i][2]) / n2;
        
        double xint = vertices[i][0] + c1 * t, yint = vertices[i][1] + c2 * t, zint = vertices[i][2] + c3 * t;
        // Check if this point is inside the tetrahedron. If not, proceed to the next step...
        if (yint <= 1.0 - zint && yint >= xint && zint >= 0 && xint >= 0) {
            referenceCoord[0] = xint;  referenceCoord[1] = yint;  referenceCoord[2] = zint;
            computeCoord(indices[0], transform, referenceCoord, finalCoord);
            velocity = smoothedSpeedFunc(finalCoord[0], finalCoord[1], finalCoord[2], timeStep);
            // distance = t * sqrt(n2);
            // WARNING: h is not correct
            // Also this fine cells is wtf.
            setup.insert(Node3D(hx * fabs(t) * sqrt(n2), velocity, indexFromFine(indices[i][0], indices[i][1], indices[i][2])));//,
                                //indexFromFineCells(finalCoord[0], finalCoord[1], finalCoord[2])));
            
            //std::cerr << "Distance is (to plane)" << (t*sqrt(n2)) << endl;
            
            continue;
        }        
        
        // Assuming that failed, next thing to try: Edge intersection.
        double minDist = 3.0;
        for (int j = 0; j < nBase; j++) {
            // Segment is y0 -> z0. x0 is vertices[i].
            // y0, z0 is base[j], base[j+1%nBase]
            // Closest point on segment is y0 + (z0-y0)t,
            // where t = (z0 - y0) dot (y0 - x0) / |z0 - y0|^2
            // want t in [0,1]
            t = 0.0;
            double L = 0.0;
            for (int q = 0; q < 3; q++) {
                double temp = base[j][q] - base[(j+1)%nBase][q];
                t += temp * (base[j][q] - vertices[i][q]);
                L += temp * temp;
            }
            t /= L;
            double dist = 0.0;
            for (int q = 0; q < 3; q++) {
                referenceCoord[q] = base[j][q] * (1 - t) + base[(j+1)%nBase][q] * t;
                double temp = referenceCoord[q] - vertices[i][q];
                dist += temp*temp;
            }
            dist = sqrt(dist);
            
            //std::cerr << "Distance to edge " << j << " is " << dist;
            //if (t < 0 || t > 1) { std::cerr << " (fail) "; }
            //std::cerr << endl;
            
            if (t >= 0 && t <= 1 && dist < minDist) {
                 minDist = dist;
                 computeCoord(indices[0], transform, referenceCoord, finalCoord);
            }
        } // for each edge of base
        
        // Also check vertices
        for (int j = 0; j < nBase; j++) {
            double dist = 0.0;
            for (int q = 0; q < 3; q++) {
                double temp = base[j][q] - vertices[i][q];
                dist += temp * temp;
            } // for each coordinate
            dist = sqrt(dist);
            if (dist < minDist) { 
                minDist = dist; 
                referenceCoord[0] = base[j][0]; referenceCoord[1] = base[j][1];  referenceCoord[2] = base[j][2];
                computeCoord(indices[0], transform, referenceCoord, finalCoord);
            }
        } // for each vertex of base
        
        //std::cerr << "Plane intersect failed, final dist found was " << minDist << endl;
        
        // WARNING: Fix h
        velocity = smoothedSpeedFunc(finalCoord[0], finalCoord[1], finalCoord[2], timeStep);
        setup.insert(Node3D(hx * minDist, velocity, indexFromFine(indices[i][0], indices[i][1], indices[i][2])));//,
                                //indexFromFineCells(finalCoord[0], finalCoord[1], finalCoord[2])));         
    } // for each vertex of tetrahedron i
}

double Interface3D::tricubicSpeedFunc(double i, double j, double k, double timeStep) {
        int I = (int)(i/ns + 0.5), J = (int)(j/ns + 0.5), K = (int)(k/ns + 0.5);
        int iStart = (I-1)*ns, jStart = (J-1)*ns, kStart = (K-1)*ns;
          
		double xPos = (i - iStart - ns / 2.0) / ((double)ns), yPos = (j - jStart - ns / 2.0) / ((double)ns), 
			   zPos = (k - kStart - ns / 2.0) / ((double)ns);

        double rhs[64];
        int xcoords[4], ycoords[4], zcoords[4];
        
		for (int c = 0; c < 4; c++) {
			xcoords[c] = I - 1 + c;
			ycoords[c] = J - 1 + c;
			zcoords[c] = K - 1 + c;

			if (xcoords[c] < 0) { xcoords[c] = 0; } if (xcoords[c] > coarseM) { xcoords[c] = coarseM; }
			if (ycoords[c] < 0) { ycoords[c] = 0; } if (ycoords[c] > coarseN) { ycoords[c] = coarseN; }
			if (zcoords[c] < 0) { zcoords[c] = 0; } if (zcoords[c] > coarseP) { zcoords[c] = coarseP; }
		}

        for (int c = 0; c < 64; c++) {
            int xc = c % 4, yc = (c / 4) % 4, zc = c / 16;
            int II = xcoords[xc], JJ = ycoords[yc], KK = zcoords[zc];
			int coarseCoord = indexFromCoarse(II, JJ, KK);
            rhs[c] = isolatedSpeedFunc(grid(coarseCoord), vActs(coarseCoord), surfaceAreas(coarseCoord), timeStep);
        }
        
        return tricubicInterpolation(xPos, yPos, zPos, rhs);
}

double Interface3D::tricubicInterpolation(double x, double y, double z, double * __restrict rhs) {
    double coeff[64];
    
    /*int offset = 0;
    for (int i = 0; i < 64; i++) {
        double coeffTmp = 0.0;
        for (int j = 0; j < 64; j++) {
            coeffTmp += tricubicMatrix[offset] * rhs[j];
            offset++;
        }
        coeff[i] = coeffTmp;
        //offset += 64;
    }*/

    /*char trans = 'n';
    integer sixtyfour = 64, one = 1;
    double onepointoh = 1.0, zeropointoh = 0.0;

    dgemv_(&trans, &sixtyfour, &sixtyfour, &onepointoh, const_cast<double *>(tricubicMatrix), &sixtyfour, rhs, &one, &zeropointoh, coeff, &one);*/

    int c = 0;
    for (int i = 0; i < 64; i++) {
        double coeffTmp = 0.0;
        for ( ; c < tricubicRowOffsets[i]; c++) {
            coeffTmp += sparseCoeffs[c] * rhs[tricubicSparseIndices[c]];
        }
        coeff[i] = coeffTmp;
    }

    //std::cerr << "Coeff vs Coeff2 " << coeff[13] << " " << coeff2[13] << endl;

    double sum = 0.0;
    for (int i = 3; i >= 0; i--) {
        //s *= x;
        double ysum = 0.0;
        for (int j = 3; j >= 0; j--) {
            //s *= ;
            double zsum = 0.0;
            for (int k = 3; k >= 0; k--) {
                zsum = z*zsum + coeff[i+4*j+16*k];
            }
            ysum = ysum*y + zsum;
        }
        sum = sum*x + ysum;
    }

    return sum;
}

void Interface3D::updateGradients(gradientElement & ge, EikonalHeap3D & workHeap) {
    // Algorithm: Compute volume flow at point (i,j,k) (store in workGrid) with some small tau
    // Compute tetrahedon changes in surrounding tetra
    // or maybe just surrounding cubes if we are lazy
    // then do get.setContribution(offsets, dvol);
    int i = ge.X(), j = ge.Y(), k = ge.Z();
    int I = i / ns + 1, J = j / ns + 1, K = k / ns + 1;
	int centerIndex = indexFromFine(i, j, k);

    double backupVal = interfaceGrid(centerIndex);
    const double tau = 1e-3;

    double right, left, up, down, forward, back, center;

	center  = interfaceGrid(centerIndex);
    right   = interfaceGrid(centerIndex + interfaceGrid.right);
    left    = interfaceGrid(centerIndex + interfaceGrid.left);
    up      = interfaceGrid(centerIndex + interfaceGrid.top);
	down    = interfaceGrid(centerIndex + interfaceGrid.bottom);
	forward = interfaceGrid(centerIndex + interfaceGrid.forward);
	back    = interfaceGrid(centerIndex + interfaceGrid.backward);                              
                       
    double ux  = (right - left)/(2*hx),
           uy  = (up - down)/(2*hy),
           uz  = (forward - back)/(2*hz);       

    double gradLen = ux*ux + uy*uy + uz*uz;
                       
    double dxp = (right - center)/hx,
           dxm = (center - left)/hx,
           dyp = (up - center)/hy,
           dym = (center - down)/hy,
           dzp = (forward - center)/hz,
           dzm = (center - back)/hz;
    double gradPlus = 0.0;
    double gPlusScratch = max(dxm*fabs(dxm), 0.0) - min(dxp*fabs(dxp), 0.0)
                        + max(dym*fabs(dym), 0.0) - min(dyp*fabs(dyp), 0.0)
                        + max(dzm*fabs(dzm), 0.0) - min(dzp*fabs(dzp), 0.0);
    if (gPlusScratch >= 0) { gradPlus = sqrt(gPlusScratch); }
    //else { std::cerr << gPlusScratch << endl; system("PAUSE");}
    double newVal = interfaceGrid(centerIndex) - tau * gradPlus;

    // Find change in volume in those boxes sharing a vertex in common with (i,j,k)

	for (Index3 inds(i-1, j-1, k-1, i, j, k, m, n, p); inds.valid(); ++inds) {

	    int Ip = 1 + inds.i / ns - I, Jp = 1 + inds.j / ns - J, Kp = 1 + inds.k / ns - K;
        int offset = offsetFromDeltas(Ip, Jp, Kp);
        
		double newVolume, oldVolume, dummyArea;
        // Change to outward flow at that vertex
        interfaceGrid(centerIndex) = newVal;
		
		newVolume = fluidAndAreaBox(interfaceGrid(inds), 
			                        interfaceGrid(inds + interfaceGrid.right),
                                    interfaceGrid(inds + interfaceGrid.top), 
									interfaceGrid(inds + interfaceGrid.right + interfaceGrid.top),                                    
									interfaceGrid(inds + interfaceGrid.forward), 
			                        interfaceGrid(inds + interfaceGrid.right + interfaceGrid.forward),
                                    interfaceGrid(inds + interfaceGrid.top + interfaceGrid.forward), 
									interfaceGrid(inds + interfaceGrid.right + interfaceGrid.top + interfaceGrid.forward), dummyArea);

        // Check original volume in cell
        interfaceGrid(centerIndex) = backupVal;
        oldVolume = fluidAndAreaBox(interfaceGrid(inds), 
			                        interfaceGrid(inds + interfaceGrid.right),
                                    interfaceGrid(inds + interfaceGrid.top), 
									interfaceGrid(inds + interfaceGrid.right + interfaceGrid.top),                                    
									interfaceGrid(inds + interfaceGrid.forward), 
			                        interfaceGrid(inds + interfaceGrid.right + interfaceGrid.forward),
                                    interfaceGrid(inds + interfaceGrid.top + interfaceGrid.forward), 
									interfaceGrid(inds + interfaceGrid.right + interfaceGrid.top + interfaceGrid.forward), dummyArea);

        ge.addContribution(offset, fabs(newVolume - oldVolume) / (tau * ns * ns * ns));
    }

    ge.normalize();
}

void Interface3D::computeVolumeOfGradients(EikonalHeap3D & workHeap) {
    int num = 0;

    // For each gradient we must do a FMM
    // but makes sense not to initialize each time
    // so that is why we input workHeap

    // We actually require a for-each-cell calculation here
	for (Index3 coarseInds(1, 1, 1, coarseM, coarseN, coarseP, coarseM+1, coarseN+1, coarseP+1); coarseInds.valid(); ++coarseInds) {
        if (surfaceAreas(coarseInds) <= 1e-14) continue;
        // Update interface in this cell by adding a constant tau = 0.01 of our volume gradients
        // Store result in workGrid
        // (f(x + tau v) - f(x)) / tau
        double tau = 0.005; //min(max(fabs(grid[I][J][K] - vActs[I][J][K]), 0.0001), 0.2);
		
		int I = coarseInds.i, J = coarseInds.j, K = coarseInds.k;
		int cellStartI = (I-1)*ns, cellStartJ = (J-1)*ns, cellStartK = (K-1)*ns;
		int clearWindowStartI = cellStartI - ns - 4, clearWindowStartJ = cellStartJ - ns - 4, clearWindowStartK = cellStartK - ns - 4,
			clearWindowEndI = cellStartI + 2 * ns + 4, clearWindowEndJ = cellStartJ + 2 * ns + 4, clearWindowEndK = cellStartK + 2 * ns + 4;

		clearWindowStartI = max(0, clearWindowStartI);
		clearWindowStartJ = max(0, clearWindowStartJ);
		clearWindowStartK = max(0, clearWindowStartK);
		clearWindowEndI = min(m, clearWindowEndI);
		clearWindowEndJ = min(n, clearWindowEndJ);
		clearWindowEndK = min(p, clearWindowEndK);

		for (Index3 clearIndices(clearWindowStartI, clearWindowStartJ, clearWindowStartK, clearWindowEndI, clearWindowEndJ, clearWindowEndK, 
			m, n, p); clearIndices.valid(); ++clearIndices) {
				workGrid(clearIndices) = INFTY;
		}

		// Now compute f(x + tau v)                    
		std::vector<gradientElement*> & cellInterfaceElts = gradientBand(coarseInds);
		for (std::vector<gradientElement*>::iterator c = cellInterfaceElts.begin(); c != cellInterfaceElts.end(); c++) {
			gradientElement * ge = (*c);
			int i = ge->X(), j = ge->Y(), k = ge->Z(),
				IStar = i/ns + 1, JStar = j/ns + 1, KStar = k/ns + 1;

			int Ip = I - IStar, Jp = J - JStar, Kp = K - KStar;
			double speed = ge->getContribution(offsetFromDeltas(Ip, Jp, Kp));
			if (isNaN(speed)) {
				cerr << "Come on, then!\n";
			}
			workHeap.insert(Node3D(ge->U(), speed, indexFromFine(i, j, k)));
		}

		doFastMarch3D(workGrid, velocities, workHeap, false, 2.0 * hx);
		workHeap.discardAll();

		int flowStartI = cellStartI-2, flowStartJ = cellStartJ-2, flowStartK = cellStartK-2,
			flowEndI = cellStartI + ns + 2, flowEndJ = cellStartJ + ns + 2, flowEndK = cellStartK + ns + 2;

		flowStartI = max(1, flowStartI);
		flowStartJ = max(1, flowStartJ);
		flowStartK = max(1, flowStartK);
		flowEndI = min(m-1, flowEndI);
		flowEndJ = min(n-1, flowEndJ);
		flowEndK = min(p-1, flowEndK);

		for (Index3 flowIter (flowStartI, flowStartJ, flowStartK, flowEndI, flowEndJ, flowEndK, m, n, p); flowIter.valid(); ++flowIter) {

			if (workGrid(flowIter) == INFTY) {
				velocities(flowIter) = 0.0;				
			}

			double right, left, up, down, forward, back, center;
            center = interfaceGrid(flowIter);
			right = interfaceGrid(flowIter + interfaceGrid.right);
            left  = interfaceGrid(flowIter + interfaceGrid.left);
			up    = interfaceGrid(flowIter + interfaceGrid.top);
            down  = interfaceGrid(flowIter + interfaceGrid.bottom);
            forward = interfaceGrid(flowIter + interfaceGrid.forward);
            back    = interfaceGrid(flowIter + interfaceGrid.backward);                            
                       
            double ux  = (right - left)/(2*hx),
                   uy  = (up - down)/(2*hy),
                   uz  = (forward - back)/(2*hz);       

            double gradLen = ux*ux + uy*uy + uz*uz;
            double dxp = (right - center)/hx,
                   dxm = (center - left)/hx,
                   dyp = (up - center)/hy,
                   dym = (center - down)/hy,
                   dzp = (forward - center)/hz,
                   dzm = (center - back)/hz;
                              
            double gradPlus = 0.0;
                       
            double gPlusScratch = max(dxm*fabs(dxm), 0.0) - min(dxp*fabs(dxp), 0.0)
                                + max(dym*fabs(dym), 0.0) - min(dyp*fabs(dyp), 0.0)
                                + max(dzm*fabs(dzm), 0.0) - min(dzp*fabs(dzp), 0.0);
            if (gPlusScratch >= 0) { gradPlus = sqrt(gPlusScratch); }
                    
			//workGrid[i][j][k] -= tau * velocities[i][j][k] * gradPlus;
			workGrid(flowIter) = interfaceGrid(flowIter) - tau * velocities(flowIter) * gradPlus;

			if (isNaN(workGrid(flowIter))) {
				cerr << "Ok dudes\n";
			}

			//interfaceGrid[i][j][k] -= tau * velocities[i][j][k] * gradPlus;
		}

		num++; 

		// Clear entries
		gradientElement & whichElement = volumeOfGradients(coarseInds);
        whichElement.clearContributions();

        // Then compute volume difference between interfaceGrid and workGrid for each tetrahedron
        // startin from iStart-1 to iStart + ns (+1 grid cell around coarse cell)
        // Add changed volume to appropriate volumeOfGradient

		int volStartI = cellStartI-1, volStartJ = cellStartJ-1, volStartK = cellStartK-1,
            volEndI = cellStartI+ns, volEndJ = cellStartJ+ns, volEndK = cellStartK+ns;

		volStartI = max(0, volStartI);
		volStartJ = max(0, volStartJ);
		volStartK = max(0, volStartK);
		volEndI = min(m-1, volEndI);
		volEndJ = min(n-1, volEndJ);
		volEndK = min(p-1, volEndK);

		for (Index3 deltaVIndices(volStartI, volStartJ, volStartK, volEndI, volEndJ, volEndK, m, n, p); deltaVIndices.valid(); ++deltaVIndices) {
			int Ip = -1 + (deltaVIndices.i - cellStartI + ns) / ns, 
				Jp = -1 + (deltaVIndices.j - cellStartJ + ns) / ns, 
				Kp = -1 + (deltaVIndices.k - cellStartK + ns) / ns;
			int index = offsetFromDeltas(Ip, Jp, Kp);

			double dummyarea = 0.0;
			double oldVolume = fluidAndAreaBox(interfaceGrid(deltaVIndices),     
				                               interfaceGrid(deltaVIndices + interfaceGrid.right), 
								               interfaceGrid(deltaVIndices + interfaceGrid.top),   
											   interfaceGrid(deltaVIndices + interfaceGrid.right + interfaceGrid.top),
							                   interfaceGrid(deltaVIndices + interfaceGrid.forward),   
											   interfaceGrid(deltaVIndices + interfaceGrid.forward + interfaceGrid.right), 
											   interfaceGrid(deltaVIndices + interfaceGrid.forward + interfaceGrid.top), 
											   interfaceGrid(deltaVIndices + interfaceGrid.forward + interfaceGrid.right + interfaceGrid.top), dummyarea);

			double newVolume = fluidAndAreaBox(workGrid(deltaVIndices),     
				                               workGrid(deltaVIndices + interfaceGrid.right), 
								               workGrid(deltaVIndices + interfaceGrid.top),   
											   workGrid(deltaVIndices + interfaceGrid.right + interfaceGrid.top),
							                   workGrid(deltaVIndices + interfaceGrid.forward),   
											   workGrid(deltaVIndices + interfaceGrid.forward + interfaceGrid.right), 
											   workGrid(deltaVIndices + interfaceGrid.forward + interfaceGrid.top), 
											   workGrid(deltaVIndices + interfaceGrid.forward + interfaceGrid.right + interfaceGrid.top), dummyarea);
			double changeInVol = (newVolume - oldVolume) / (tau * ns * ns * ns);


			whichElement.addContribution(index, changeInVol);
		} // Done with changed volume

		// not quite the right fix
		if (whichElement.getContribution(middleCell) < 0.0) {
			std::cerr << "Negative: " << num << endl;
			double fix = surfaceAreas(coarseInds) / 6.0 - whichElement.getContribution(middleCell);
			whichElement.addContribution(middleCell, fix);
		}

		for (Index3 clearIndices(clearWindowStartI, clearWindowStartJ, clearWindowStartK, clearWindowEndI, clearWindowEndJ, clearWindowEndK, 
			m, n, p); clearIndices.valid(); ++clearIndices) {
				workGrid(clearIndices) = INFTY;
		}
	}
}

void Interface3D::computeLambdas() {
    // Should be a sparse matrix, but since it's on the coarse grid,
    // maybe the dense version won't kill us.
    //std::cerr << "Computing \\lambda\n";
    int * matrixRow = new int[(coarseM+2)*(coarseN+2)*(coarseP+2)];
    double * A, *rhs;

    int rows = 0;

	for (int inds = 0; inds < coarseSize; inds++) {
		lambdas(inds) = 0.0;
        residual(inds) = grid(inds) - vActs(inds);
		double diag = fabs(volumeOfGradients(inds).getContribution(middleCell));
		
		if (surfaceAreas(inds) >= 1e-14 && diag > 0) { 
			matrixRow[rows++] = inds; 
		}
     }

     // Create the matrix
     A = new double[rows * rows];
     rhs = new double[rows];

     //ofstream matrixFile("AGMatrix.txt");
     //ofstream residualFile("Residual.txt");
     
     for (int r = 0; r < rows; r++) {
         // ijkID counts Z, then J, then I
         int indR = matrixRow[r],
			 IR = indR / ((coarseN+2) * (coarseP+2)), 
             JR = (indR / (coarseP+2)) % (coarseN+2), 
             KR = indR % (coarseP+2);

         for (int c = 0; c < rows; c++) {
             int indC = matrixRow[c],
				 IC = indC / ((coarseN+2) * (coarseP+2)), 
                 JC = (indC / (coarseP+2)) % (coarseN+2),
                 KC = indC % (coarseP+2);

             double val = 0.0;
             if (abs(IC - IR) <= 1 && abs(JC - JR) <= 1 && abs(KC - KR) <= 1) {
                 int whichCell = middleCell + (IC - IR) + 3 * (JC - JR) + 9 * (KC - KR);
                 if (fabs(volumeOfGradients(indC).getContribution(middleCell)) > 0.1 || whichCell == middleCell) {
                     val = volumeOfGradients(indR).getContribution(whichCell); 
                 }
             }
             A[c * rows + r] = val;
             //matrixFile << val << " ";
         }		 
         //matrixFile << endl;
         rhs[r] = residual(indR);
         //residualFile << rhs[r] << endl;
     }
     
     //matrixFile.close();
     //residualFile.close();
     
     integer numEq = rows, numRhs = 1, ldA = rows, ldB = rows, info;
     integer * pivotArray = new integer[rows];

     dgesv_(&numEq, &numRhs, A, &ldA, pivotArray, rhs, &ldB, &info);

     //ofstream lhs("ComputedX.txt");
     for (int r = 0; r < rows; r++) {
         int indR = matrixRow[r],
			 IR = indR / ((coarseN+2) * (coarseP+2)), 
             JR = (indR / (coarseP+2)) % (coarseN+2), 
             KR = indR % (coarseP+2);
         //const double maxSpeed = 0.1 * hx;//0.02 * hx;
         double maxSpeed = surfaceAreas(indR) * 0.25 * hx;
         if (rhs[r] > maxSpeed) { 
             rhs[r] = maxSpeed; 
         }
         else if (rhs[r] < -maxSpeed) { 
             rhs[r] = -maxSpeed; 
         }
         //lhs << rhs[r] << endl;

         lambdas(indR) = rhs[r];
         //lambdas[IR][JR][KR] = residual[IR][JR][KR] / volumeOfGradients[IR][JR][KR][13];		
     }
     //lhs.close();

     delete[] A;
     delete[] rhs;
     delete[] pivotArray;
     delete[] matrixRow;    

     //std::cerr << "Lambda computed\n";
}

void Interface3D::gaussSeidelLambdas(int numIters) {
    double L1Res = 0.0;
    // Set residuals to Vr - Va, and lambdas to 0
	for (int inds = 0; inds < coarseSize; inds++) {
		double diag = volumeOfGradients(inds).getContribution(middleCell);
        if (diag < 0.0 && surfaceAreas(inds) > 1e-8) { std::cerr << "Diag = " << diag << "\n"; }
                
        lambdas(inds) = 0.0;
        double r = grid(inds) - vActs(inds);
        residual(inds) = r;
        L1Res += fabs(r);
	}

    std::cerr << "Residual starts at " << L1Res << endl;

	const int right = residual.right, top = residual.top, forward = residual.forward;

    for (int step = 0; step < numIters; step++) {

        for (int inds = 0; inds < coarseSize; inds++) {
			double diag = volumeOfGradients(inds).getContribution(middleCell);
            if (surfaceAreas(inds) <= 1e-8 || fabs(diag) <= 1e-8) {
                continue;
            }

			const int center = (int)inds;

            double lambdaChange = residual(inds) / diag;
            const double maxSpeed = 0.1 * hx;
            if (lambdaChange > maxSpeed) { lambdaChange = maxSpeed; }
            if (lambdaChange < -maxSpeed) { lambdaChange = -maxSpeed; }

            int offset = 0;
            // Update residuals
            for (int offset = 0; offset < 27; offset++) {
                int Ip, Jp, Kp;
                IndexDeltasFromOffset(Ip, Jp, Kp, offset);
                int neighbor = center + Ip * right + Jp * top + Kp * forward;

                if (surfaceAreas(neighbor) <= 1e-8) {
                    continue;
                }
				residual(neighbor) -= volumeOfGradients(neighbor).getContribution(26 - offset) * lambdaChange;
            } // Updating residuals
            lambdas(inds) += lambdaChange;
		}

        // Check residual (debug)
        L1Res = 0.0;
        for (int inds = 0; inds < coarseSize; inds++) {
            L1Res += fabs(residual(inds));
        } // done w residual
        std::cerr << "Residual " << L1Res << endl;
    }

	for (int inds = 0; inds < coarseSize; inds++) {
        double area = surfaceAreas(inds);
        if (area < 1e-8) { continue; }

        double maxSpeed = 0.1 * hx;
                
        if (lambdas(inds) > maxSpeed) { lambdas(inds) = maxSpeed; }
        else if (lambdas(inds) < -maxSpeed) { lambdas(inds) = -maxSpeed; }
    }
}

void Interface3D::checkResult(bool verbose){
      
       double maxError = 0.0, totalRequiredVol = 0.0, totalComputedVol = 0.0,
              sumOfErrors = 0.0;
              
       computeVActs();
       
       int mixed = 0;
       
	   for (Index3 inds(1,1,1, coarseM, coarseN, coarseP, coarseM+1, coarseN+1, coarseP+1); inds.valid(); ++inds) {
           double vAct = vActs(inds);
           double vReq = grid(inds);
               
           totalComputedVol += vAct;
           totalRequiredVol += vReq;
               
           if ((vReq > 1e-8 && vReq < 1 - 1e-8) ||
               (vAct > 1e-8 && vAct < 1 - 1e-8)) {
               mixed++;
               if (verbose) { cout << "Cell (" << inds.i << "," << inds.j << "," << inds.k << "): required = " << vReq << ", actual = " << vAct << endl; }
           }
           double error = fabs(vReq - vAct);
           sumOfErrors += error;
           if (error > maxError) { maxError = error; }
       }        
              
       // Count faces
       int faces = 0;
	   for (Index3 inds(0,0,0, m-1, n-1, p-1, m, n, p); inds.valid(); ++inds) {
		   int center = (int) inds;
		   const int forward = interfaceGrid.forward, right = interfaceGrid.right, top = interfaceGrid.top;
           // (0,0,0), (0,1,0), (1,1,0), (0,0,1)
           double dummyArea = 0.0;
           if (fluidTetrahedron(interfaceGrid(center), interfaceGrid(center + top), 
			                    interfaceGrid(center + right + top), interfaceGrid(center + forward), dummyArea) > 1e-9) {
               faces++;
           }
		   if(fluidTetrahedron(interfaceGrid(center + forward + top), interfaceGrid(center + top),
                               interfaceGrid(center + right + top), interfaceGrid(center + forward), dummyArea) > 1e-9) {
               faces++;
           }
		   if(fluidTetrahedron(interfaceGrid(center + forward + top + right), interfaceGrid(center + forward + top),
                               interfaceGrid(center + forward), interfaceGrid(center + right + top), dummyArea) > 1e-9) {
               faces++;
           }                   
           if(fluidTetrahedron(interfaceGrid(center + right + top + forward), interfaceGrid(center + right + forward), 
                               interfaceGrid(center + forward), interfaceGrid(center + right + top), dummyArea) > 1e-9) {
               faces++;
           }                     
           if (fluidTetrahedron(interfaceGrid(center + right), interfaceGrid(center + forward + right),
                               interfaceGrid(center + forward), interfaceGrid(center + right + top), dummyArea) > 1e-9) {
               faces++;
           }
           if(fluidTetrahedron(interfaceGrid(center), interfaceGrid(center + right),
                               interfaceGrid(center + right + top), interfaceGrid(center + forward), dummyArea) > 1e-9) {
               faces++;
           }
       }
       
       
       cout << "Max error = " << maxError << "; L1 error = " << sumOfErrors << "\nTotal volume in interface = " 
            << totalComputedVol << "; required vol = " << totalRequiredVol << endl;
            
       std::cerr << "Mixed " << ((double)mixed/(coarseM * coarseN * coarseP)) << " | Faces " << faces << endl;
       std::cerr << "Average error = " << (sumOfErrors/mixed) << endl;
}

double Interface3D::calculateVolumeImmersed(
                    int I, int J, int K) {
       double vol = 0.0, area = 0.0;
       int span = ns;
       double unitVol = 1.0 / (span*span*span), unitArea = 1.0 / (span*span);
       
       int startIndexX = (I-1)*ns, startIndexY = (J-1)*ns, startIndexZ = (K-1)*ns;
       
	   for (Index3 inds(startIndexX, startIndexY, startIndexZ, startIndexX + span - 1, startIndexY + span - 1, startIndexZ + span - 1,
		   m, n, p); inds.valid(); ++inds) {

			   const int center = (int)inds, right = interfaceGrid.right, top = interfaceGrid.top, forward = interfaceGrid.forward;
			   // 000, 100, 010, 110, 001, 101, 011, 111
			   vol += fluidAndAreaBox(interfaceGrid(center), interfaceGrid(center + right), 
                                      interfaceGrid(center + top), interfaceGrid(center + right + top),
									  interfaceGrid(center + forward), interfaceGrid(center + right + forward), 
                                      interfaceGrid(center + top + forward), 
									  interfaceGrid(center + right + top + forward), 
                                      area);                     
       }
       
       vol *= unitVol;
       area *= unitArea;

	   int IJK = indexFromCoarse(I, J, K);
              
       double change = fabs(vol - vActs(IJK));
       
       vActs(IJK) = vol;
       surfaceAreas(IJK) = area;

       totalSurfaceArea += area;
       
       return change;
}

Interface3D::~Interface3D() {
    // Delete all that memory we allocated.
	// Handled by destructors
}
 
void Interface3D::normalToInterface(int i, int j, int k, 
                                        double & n1, double & n2, double & n3) {
     int numVectors = 0;
     n1 = 0.0; n2 = 0.0; n3 = 0.0;
     double n1t, n2t, n3t;
           
     double A, B, C, D;
     // (0,0,0), (0,1,0), (1,1,0), (0,0,1)
     A = interfaceGrid(i,j,k); B = interfaceGrid(i,j+1,k); 
     C = interfaceGrid(i+1,j+1,k); D = interfaceGrid(i,j,k+1);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = C-B;
        n2t = B-A;
        n3t = D-A;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        std::cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }           
     // (1,1,1), (0,1,1), (0,0,1), (1,1,0)
     A = interfaceGrid(i+1,j+1,k+1); B = interfaceGrid(i,j+1,k+1); 
     C = interfaceGrid(i,j,k+1); D = interfaceGrid(i+1,j+1,k);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = -B+A;//C-B;
        n2t = B-C;//B;
        n3t = -D+A;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        std::cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";        
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }   
     // (1,1,1), (1,0,1), (0,0,1), (1,1,0)
     A = interfaceGrid(i+1,j+1,k+1); B = interfaceGrid(i+1,j,k+1); 
     C = interfaceGrid(i,j,k+1); D = interfaceGrid(i+1,j+1,k);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = B-C;
        n2t = -B+A;
        n3t = -D+A;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        std::cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";        
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }   
     // (1,0,0), (1,0,1), (0,0,1), (1,1,0)
     A = interfaceGrid(i+1,j,k); B = interfaceGrid(i+1,j,k+1); 
     C = interfaceGrid(i,j,k+1); D = interfaceGrid(i+1,j+1,k);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = B-C;
        n2t = D-A;
        n3t = B-A;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        std::cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";        
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }   
     // (0,1,1), (0,1,0), (1,1,0), (0,0,1)
     A = interfaceGrid(i,j+1,k+1); B = interfaceGrid(i,j+1,k); 
     C = interfaceGrid(i+1,j+1,k); D = interfaceGrid(i,j,k+1);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = C-B;
        n2t = -D+A; //B;
        n3t = -B+A; //D;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        std::cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";        
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }   
     // (0,0,0), (1,0,0), (1,1,0), (0,0,1)
     A = interfaceGrid(i,j,k); B = interfaceGrid(i+1,j,k); 
     C = interfaceGrid(i+1,j+1,k); D = interfaceGrid(i,j,k+1);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = B-A;//C-B;
        n2t = C-B;//B;
        n3t = D-A;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        std::cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";        
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }              
           
     if (numVectors > 0) {
        n1 /= numVectors;
        n2 /= numVectors;
        n3 /= numVectors;                                            
     } else {
        std::cerr << "No normal in this cell...\n";
        //for (int ii = i-2; ii <= i+2; ii++) {
        //    for (int jj = j-2; jj <= j+2; jj++) {
        //        for (int kk = k-2; kk <= k+2; kk++) {
        //            std::cerr << interfaceGrid[ii][jj][kk] << " ";
        //        }
        //        std::cerr << endl;
        //    }
        //    std::cerr << endl;
        //}
     }
}

void Interface3D::saveInterface() {
    // For each node in the narrow band
    for (int nbInd = 0; nbInd < narrowBandSize; nbInd++) {
		int ind = narrowBandIterator[nbInd].Index();
        oldInterface(ind) = interfaceGrid(ind);
    }
}

double Interface3D::interfaceChange() {
    double change = 0.0;
    
    // For each node in the narrow band
    for (int nbInd = 0; nbInd < narrowBandSize; nbInd++) {
        Node3D nbElt = narrowBandIterator[nbInd];
		int ind = nbElt.Index();

        double phi = interfaceGrid(ind), 
               psi =  oldInterface(ind);
        change += fabs(phi - psi) * exp(-0.125 * psi * psi / (hx*hx));
    }

    return change * hx * hy * hz;
}
 