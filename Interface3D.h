#include "stdafx.h"
#include "Heap3D.h"
#include "gradientElement.h"
#include "MultiGrid.h"

using namespace std;

// Commented out functions have not yet been updated.

typedef std::vector<gradientElement*> gradientList;

class Interface3D {

	  int Material;
      int coarseSize, fineSize;

      // (numPhases x) coarseM x coarseN x coarseP
      MultiGrid3<double> grid;
      MultiGrid3<double> vActs, vActsBack, surfaceAreas;
	  Grid3<double> lambdas, residual, curvatureSums;

      // (numPhases x) m x n x p
      MultiGrid3<double> interfaceGrid, oldInterface;
	  Grid3<double> workGrid, velocities, signs;
	  Grid3<gradientElement> volumeGradients; // Volume that each gradient element occupies in own cell and neighbors

	  // coarseM x coarseN x coarseP x 27
	  Grid3<gradientElement> volumeOfGradients; // Matrix. Volume that each volume gradient function occupies within its own cell and neighbors.
	  Grid3<gradientList> gradientBand; // Is a smaller narrow band, for each cell, where gradients will lie

	  std::vector<Node3D> allGradients;
	  int numGradElts, capacityGradElts;

	  // Computes the speed function constants for each coarse grid cell
	  // That is, moves the interface by Lambda[I][J][K] * VolGradient[I][J][K]
	  // For now VolGradient = const function
	  void computeVolumeOfGradients(EikonalHeap3D & workHeap);
	  void computeLambdas();
	  void gaussSeidelLambdas(int numIters);
	  void updateGradients(gradientElement & ge, EikonalHeap3D & workHeap);

	  // Level set methods ////////////////////////
	  // Outward normal speed function plus curvature flow
	  double computeVolumeFlow(double epsilon, double timeStep);
	  // Curvature flow only
      double computeCurvatureFlow(double epsilon);
	  void curvatureAt(int centerIndex, int i, int j, int k, double & curvature, double & gradLen);
	  void curvatureAt(int i, int j, int k, double & curvature, double & gradLen) {
		  curvatureAt(indexFromFine(i, j, k), i, j, k, curvature, gradLen);
	  }
	  // Updates phi = phi + timeStep * workGrid
	  void updateInterface(double timeStep);

	  MultiGrid3<bool> narrowBand, coarseNarrowBand; // One m,n,p, the other coarseM, coarseN, coarseP
	  int * narrowBandSize;
	  Node3D ** narrowBandIterator;
	  void updateCoarseNarrowBand();

	  // Updates so that there is no overlap between multiple level set functions
	  void fixLevelsetOverlap();
	  void fixLevelsetOverlapNoNB();
	  //////////////////////////////////////////////

      // Calculating indices
	  void coarseFromIndex(int ind, int & I, int & J, int & K) { 
		  K = ind % (coarseP+2); J = (ind / (coarseP+2)) % (coarseN+2); I = ind / ((coarseN+2)*(coarseP+2)); 
	  }
	  void fineFromIndex  (int ind, int & i, int & j, int & k) { 
		  k = ind % (p+1); j = (ind / (p+1)) % (n+1); i = ind / ((n+1)*(p+1)); 
	  }
	  int indexFromCoarse(int I, int J, int K) { return K + (coarseP+2) * (J + (coarseN+2)*I); }
	  int indexFromFine  (int i, int j, int k) { return k + (p+1) * (j + (n+1)*i); }

	  int coarseIndexFromFine(int i, int j, int k) {
		  return (k/ns) + 1 + (coarseP+2) * ((j/ns)+1 + (coarseN+2)*((i/ns)+1));
	  }

      // Maps tetrahedron rotated coordinates back onto (i,j,k) space (extended to real numbers)
      void computeCoord(int tetStart[3], int transform[3], double coords[3], double result[3]) {
          for (int i = 0; i < 3; i++) {
              int dir = transform[i] > 0 ? 1 : -1, whichCoord = abs(transform[i]) - 1;
              result[whichCoord] = tetStart[whichCoord] + dir * coords[i];
          }                           
      }

      // Given a tetrahedron of four points (indices), computes curvature over the endpoints of the
	  // polygon defining the interface within that tetrahedron, and takes the maximum.
      double getMaxCurvatureOnInterface(int indA, int indB, int indC, int indD);
	  // faceListIn: An array of polygons (which are each an array of points, which are an array of 3 coordinates)
	  // faceListOut: where the new faces will go
	  // Chop faceListIn by the plane ax + by + cz + d = 0
	  void chopPolytope(double Verts[][3], int numVerts, int ** faces, int nFaceIn, 
	                    double newVerts[][3], int & numNewVerts, int ** newFaces, int & nFaceOut,
	                    double a, double b, double c, double d);
	  double volumeOfPolytope(const double Vertices[][3], int numVerts, int ** faces, int numFaces);

      //double computeVolumeFlow(double epsilon, double timeStep, bool precise = false);
      //double computeCurvatureFlow(double epsilon);
      //double calculateBestTimeStep(double epsilon, double k, double mag1);                                              
      //double objective(double epsilon);     
      double calculateVolumeImmersed(int I, int J, int K);      
      
      double isolatedSpeedFunc(double vReq, double vAct, double area, double timeStep) { 
             const double delta = 1e-3;
             double error = vReq - vAct, sign = error >= 0 ? 1.0 : -1.0;
			 // New: Want to cap maximum flow rate as relative to surface area in cell
			 // Fixes problems with some lumps
             return sign * min( fabs(error) / (area * 5.0 * timeStep + 1e-20), 0.1 * hx * area / timeStep );
			 //return sign * min( fabs(error) / (area * 5.0 * timeStep + 1e-20), 0.25 * hx / timeStep );
			 //return 0;
      }
      
      double smoothedSpeedFunc(double i, double j, double k, double timeStep);
	  double tricubicSpeedFunc(double i, double j, double k, double timeStep);

	  double tricubicInterpolation(double x, double y, double z, double * __restrict rhs);
	  double sigmoid(double, double);
      
	  void tricubicCoeff(double * coeff, double * __restrict rhs);
	  double tricubicEval(double x, double y, double z, double * __restrict coeff);
	  void tricubicDerivatives(double x, double y, double z, double & dx, double & dy, double & dz, double * __restrict coeff);
	  void tricubicDistanceSetup(double x0, double y0, double z0, int i, int j, int k, double * __restrict coeff, double timeStep, EikonalHeap3D & setup);

	  // Hopefully can be removed once reinitialization improves
	  double safeExtrapolate(int mat, int ic, int jc, int kc, int ir, int jr, int kr) {
		  double r = interfaceGrid(mat, ir,jr,kr), c = interfaceGrid(mat, ic, jc, kc);
		  if (narrowBand(mat, ir,jr,kr) && narrowBand(mat, ir,jr,kr) && fabs(r - c) <= 1.5 * hx) {
			  return 2.0 * c - r;
		  } else return c;
	  }

      static const double vertices[4][3];
      static const int oppositeEdge[6];
      static const int edgesTouchingVertices[4][3];
	  static int ** unitTetFaces;

	  int tetrahedraOffsets[6][4];
      
      void fillIndices(int indices[4][3], int a0, int b0, int c0, int a1, int b1, 
                       int c1, int a2, int b2, int c2, int a3, int b3, int c3) {
          indices[0][0] = a0;  indices[0][1] = b0;  indices[0][2] = c0;
          indices[1][0] = a1;  indices[1][1] = b1;  indices[1][2] = c1;
          indices[2][0] = a2;  indices[2][1] = b2;  indices[2][2] = c2;
          indices[3][0] = a3;  indices[3][1] = b3;  indices[3][2] = c3;
      }

      double areaTriangle(double A[3], double B[3], double C[3]) {
		  double a1 = B[0] - A[0], b1 = B[1] - A[1], c1 = B[2] - A[2],
			     a2 = C[0] - A[0], b2 = C[1] - A[1], c2 = C[2] - A[2];
		  double cross1 = b1 * c2 - c1 * b2, cross2 = c1 * a2 - a1 * c2, cross3 = a1 * b2 - a2 * b1;
		  double crossSq = cross1*cross1 + cross2*cross2 + cross3*cross3;
          /*double a = sqrt( (B[0]-A[0])*(B[0]-A[0]) + (B[1]-A[1])*(B[1]-A[1]) + (B[2]-A[2])*(B[2]-A[2]) ),
                 b = sqrt( (C[0]-A[0])*(C[0]-A[0]) + (C[1]-A[1])*(C[1]-A[1]) + (C[2]-A[2])*(C[2]-A[2]) ),
                 c = sqrt( (C[0]-B[0])*(C[0]-B[0]) + (C[1]-B[1])*(C[1]-B[1]) + (C[2]-B[2])*(C[2]-B[2]) ),
                 s = 0.5 * (a+b+c);
          double asq = s*(s-a)*(s-b)*(s-c);
          return asq > 0.0 ? sqrt(asq) : 0.0;*/
		  return crossSq > 0.0 ? 0.5 * sqrt(crossSq) : 0.0;
      }

	  // For the tetrahedron defined by (0,0,0), (0,1,0), (1,1,0), (0,0,1), assuming f(0,0,0) = a, f(0,1,0) = b, etc,
	  // and that f = c_0 + c_1 x + c_2 y + c_3 z, finds the plane surface f(x,y,z) = 0.
	  void getInterfacePolygon(double a, double b, double c, double d, int & nBase, double base[4][3]);
      
      // Updates the interface with the contents of workGrid
      //void updateInterface(double timeStep);
      
      void doFastMarch3D(Grid3<double> & T, Grid3<double> & V, EikonalHeap3D & initialConditions, bool reinitialize, double bandRadius = INFTY);
      // Helper function for fast marching method
      Node3D TrialNode(Grid3<double> & T, Grid3<double> & V, int centerIndex, int i, int j, int k);
                                    
      void normalToInterface(int mat, int i, int j, int k, double & n1, double & n2, double & n3);

	  static const int middleCell = 13;

	  void IndexDeltasFromOffset(int & Ip, int & Jp, int & Kp, int offset) {
		  Ip = -1 + (offset % 3);
		  Jp = -1 + ((offset % 9) / 3);
		  Kp = -1 + offset / 9;
	  }

	  int offsetFromDeltas(int Ip, int Jp, int Kp) {
		  return 9 * (Kp + 1) + 3 * (Jp + 1) + (Ip + 1);
	  }

	  int indexFromFineCells(double i, double j, double k) {
		  int I = ((int) i) / ns + 1,
			  J = ((int) j) / ns + 1,
			  K = ((int) k) / ns + 1;
		  return indexFromCoarse(I,J,K);
	  }
            
      double totalSurfaceArea;
      int iter;

	  // Setup routines ///////////////////
      void allocate();
	  void setupGhostFractions();
	  void updateGhostLayer();
	  void computeVolfracsFromPhi(char * exactFilename = __nullptr);
	  void adjustVolumes(MultiGrid3<double> & whichGrid);
      
public:
	  
      const int coarseM, coarseN, coarseP, m, n, p, ns;
	  const int numPhases;
      double hx, hy, hz;
      
      // Allocates memory and sets up (reads) volume fractions
      Interface3D(int nPhases, int gm, int gn, int gp, int _ns, double _h, 
                                  ifstream & input);
      // Only allocates data
      Interface3D(int nPhases, int gm, int gn, int gp, int _ns);
      //InterfaceFromVolumeFraction(const InterfaceFromVolumeFraction & other);
      ~Interface3D();
      
      void seedGrid(double levelSet);
      void advancedSeeding();              
      
      double iterateLevelSetMethod(double epsilon, double timeStep, bool precise = false);
      
      // Reinitializes the interface to the signed distance function
      void extensions(bool reinitialize, double timeStep, bool precise);
	  void reinitialize() {
		  for (Material = 0; Material < numPhases; Material++) { 
			  extensions(true, 0.01, false); 
			  cerr << "Narrow band size (Mat " << Material << "): " << narrowBandSize[Material] << endl;
		  }
	  }
      
      // Tells you how good your answer is.
      void checkResult(bool verbose = true, bool voronoiVol = false);
      
      // In file output.cpp
      //void drawInterface(const char * filename);
      void outputInterface(const char * filename);     
	  void outputNarrowBand(const char * filename); // just for debugs basically
      void outputInterfaceCoarse(const char * filename);
	  void outputVelocities(const char * filename);

	  void outputInterfaceMesh(const char * filename);
	  void addPolygonsForTetrahedron(std::vector<double> & points, std::vector<int> & polySizes, std::vector<int> & faceTypes,
		                             std::vector<double> & errorVals, int indices[4][3], int transform[3], int coarseInd);
      
      // Test case functions
      void initCube(double length, double x = 0, double y = 0, double z = 0, 
                    double alpha = 0, double beta = 0, double gamma = 0, double graphRange = 6.0);    
	  void initSphere(double radius, double x = 0, double y = 0, double z = 0, double graphRange = 6.0);
      void initDeltoidish(double radius, double x = 0.0, double y = 0.0, double z = 0.0, 
                          double alpha = 0.0, double beta = 0.0, double gamma = 0.0, double graphRange = 6.0);              
	  void initCone(double radius, double height, double x = 0.0, double y = 0.0, double z = 0.0,
		                  double alpha = 0.0, double beta = 0.0, double gamma = 0.0, double graphRange = 6.0);
	  void initTripleLine(double x = 0.0, double y = 0.0, double z = 0.0, double alpha = 0.0, double beta = 0.0, double gamma = 0.0, bool inSphere = false);
	  void initTripleLineFlat(double x = 0.0, double y = 0.0, double z = 0.0, double alpha = 0.0, double beta = 0.0, double gamma = 0.0);
	  void initQuadruplePoint(double x = 0.0, double y = 0.0, double z = 0.0, double alpha = 0.0, double beta = 0.0, double gamma = 0.0);
	  // The following have no parameters because they are from the paper Smooth, Volume-Accurate Material Interface Reconstruction
	  void initHemispheres(double radius, double x, double y, double z, 
                                 double alpha, double beta, double gamma, double graphRange);
	  void initCubeSphereIntersect(double alpha = 0.0, double beta = 0.0, double gamma = 0.0);
	  void initConcentricCircles();
	  void initTetra(double x = 0, double y = 0, double z = 0, double alpha = 0, double beta = 0, double gamma = 0);
	  void initRandomVoronoi(int seed, bool swirl = false);

      static Interface3D * SFBay(int ns);

      // TO DO: Figure out if it's fast enough.
	  // Note that surfaceArea is an accumulator (it will add new surfaceArea to old value)
      double fluidTetrahedron(double a, double b, double c, double d, double & surfaceArea) {
           const double oneSixth = 1.0 / 6.0;
           int nMinus = 0;

		   int sgnA = 1, sgnB = 1, sgnC = 1, sgnD = 1;
           if (a <= 0) { sgnA = -1; nMinus++; }
           if (b <= 0) { sgnB = -1; nMinus++; }
           if (c <= 0) { sgnC = -1; nMinus++; }
           if (d <= 0) { sgnD = -1; nMinus++; }
           
           if (nMinus == 0) { return 0; }
           else if (nMinus == 4) { return oneSixth; }
           else {
                double v[6][3];
                bool edges[6];
           
                if (sgnA * sgnB <= 0) { 
                     v[0][0] = 0.0;
                     v[0][1] = fabs(a) > 1e-20 ? -a/(b-a) : 0.0;
                     v[0][2] = 0.0;
                     edges[0] = true;
                } else { edges[0] = false; }
           
                if (sgnA * sgnC <= 0) {
                     v[1][0] = fabs(a) > 1e-20 ? -a/(c-a) : 0.0;
                     v[1][1] = v[1][0];
                     v[1][2] = 0.0;
                     edges[1] = true;
                } else { edges[1] = false; }
           
                if (sgnA * sgnD <= 0) {
                     v[2][0] = 0.0;
                     v[2][1] = 0.0;
                     v[2][2] = fabs(a) > 1e-20 ? -a/(d-a) : 0.0;
                     edges[2] = true;
                } else { edges[2] = false; }
            
                if (sgnB * sgnC <= 0) {
                     v[3][0] = fabs(b) > 1e-20 ? -b/(c-b) : 0.0;
                     v[3][1] = 1.0;
                     v[3][2] = 0.0;
                     edges[3] = true;
                } else { edges[3] = false; }
           
                if (sgnB * sgnD <= 0) {
                     v[4][0] = 0.0;
                     v[4][1] = fabs(d) > 1e-20 ? d/(d-b) : 0.0;
                     v[4][2] = fabs(b) > 1e-20 ? -b/(d-b) : 0.0;
                     edges[4] = true;
                } else { edges[4] = false; }
           
                if (sgnC * sgnD <= 0) {
                     v[5][0] = fabs(d) > 1e-20 ? d/(d-c) : 0.0;
                     v[5][1] = v[5][0];
                     v[5][2] = fabs(c) > 1e-20 ? -c/(d-c) : 0.0;
                     edges[5] = true;
                } else { edges[5] = false; }
                
                if (nMinus == 2) { // Volume is of a pyramid thingie + a tetrahedron
                     // Six relevant vertices are the four intersect vertices,
                     // plus the two that are minus.
                     int fifth, sixth;
                     if (a <= 0) { fifth = 0; }
                     else if (b <= 0) { fifth = 1; }
                     else { fifth = 2; }
                     
                     if (b <= 0 && fifth != 1) { sixth = 1; }
                     else if (c <= 0 && fifth != 2) { sixth = 2; }
                     else { sixth = 3; }
                     
                     // The pyramid is the four intersect vertices + vertices[fifth]
                     // Volume is 1/3 * area(base) * height
                     // height = vertices[fifth] - any of the vertices dot normal
                     double base[4][3]; int e0;
                     for (e0 = 0; e0 < 3; e0++) {
                         if (edges[e0]) { 
                            base[1][0] = v[e0][0];
                            base[1][1] = v[e0][1];
                            base[1][2] = v[e0][2];
                            break;
                         }
                     }
                     int j = 0;
                     for (int e = e0+1; e < 6; e++) {
                         if (edges[e] && oppositeEdge[e] != e0 && j < 2) {
							 if (2*j == 4) {
								 cerr << "Base corrupted\n";
							 }
                              base[2*j][0] = v[e][0];
                              base[2*j][1] = v[e][1];
                              base[2*j][2] = v[e][2];
                              j++;
                         }
                     }
                     
                     base[3][0] = v[oppositeEdge[e0]][0];
                     base[3][1] = v[oppositeEdge[e0]][1];
                     base[3][2] = v[oppositeEdge[e0]][2];                     
                     
                     double areaBase = areaTriangle(base[0], base[1], base[2]) + 
                                       areaTriangle(base[2], base[3], base[0]);
                                       
                     // Return second argument of surface area of interface.
                     surfaceArea += areaBase;
                     
                     // Recall the interpolant is
                     // f(x,y,z) = a + (c-b)x + (b-a)y + (d-a)z
                     // and so the normal vector is ((c-b), (b-a), (d-a))                                       
                     double normalLength = sqrt( (b-a)*(b-a) + (c-b)*(c-b) + (d-a)*(d-a) );
                     if (normalLength < 1e-20) { normalLength = 1.0; } // for safety
                     double height = ((c-b) * (vertices[fifth][0] - base[0][0]) +
                                      (b-a) * (vertices[fifth][1] - base[0][1]) +
                                      (d-a) * (vertices[fifth][2] - base[0][2])) / normalLength;
                     height = fabs(height);
                                          
                     // The tetrahedron vertices[sixth] + vertices[fifth] + 
                     // the two vertices (out of v0 to v3) that are on an edge
                     // touching the sixth vertex
                     j = 0;
                     for (int i = 0; i < 3; i++) {
                         int t = edgesTouchingVertices[sixth][i];
                         if (edges[t]) {
							 if (j >= 4) {
								 cerr << "Base corrupted (2nd place)\n";
							 }
                              base[j][0] = v[t][0];
                              base[j][1] = v[t][1];
                              base[j][2] = v[t][2];
                              j++;
                         }
                     }
                     base[2][0] = vertices[fifth][0];
                     base[2][1] = vertices[fifth][1];
                     base[2][2] = vertices[fifth][2];
                     base[3][0] = vertices[sixth][0];
                     base[3][1] = vertices[sixth][1];
                     base[3][2] = vertices[sixth][2];
                     
                     // Compute volume of this crap
                     for (int i = 1; i < 4; i++) {
                         for (j = 0; j < 3; j++) {
                             base[i][j] -= base[0][j];
                         }                                          
                     }
                     
                     // triple product
                     double tripleProd = base[1][0] * (base[2][1] * base[3][2] - base[2][2] * base[3][1])
                                       - base[1][1] * (base[2][0] * base[3][2] - base[2][2] * base[3][0])
                                       + base[1][2] * (base[2][0] * base[3][1] - base[2][1] * base[3][0]);
                                       
                     return oneSixth * (fabs(tripleProd) + 2.0 * areaBase * height);
                     
                } else {
                     double tet[3][3];
                     int fourth;
                     if (nMinus == 1) {
                        if (a <= 0) { fourth = 0; }
                        else if (b <= 0) { fourth = 1; }
                        else if (c <= 0) { fourth = 2; }
                        else if (d <= 0) { fourth = 3; }
                     } else {
                        if (a > 0) { fourth = 0; }
                        else if (b > 0) { fourth = 1; }
                        else if (c > 0) { fourth = 2; }
                        else if (d > 0) { fourth = 3; }                            
                     }
                     
                     int j = 0;
                     for (int i = 0; i < 6; i++) {
                         if (edges[i]) {
                              tet[j][0] = v[i][0] - vertices[fourth][0];
                              tet[j][1] = v[i][1] - vertices[fourth][1];
                              tet[j][2] = v[i][2] - vertices[fourth][2];
                              j++;
                         }
                         if (j > 2) break;
                     }
                     
                     // Compute area (if this turns out to be a good idea, will design two output argument function later)
                     surfaceArea += areaTriangle(tet[0], tet[1], tet[2]);                     
                     
                     double tripleProd = tet[0][0] * (tet[1][1] * tet[2][2] - tet[1][2] * tet[2][1])
                                       - tet[0][1] * (tet[1][0] * tet[2][2] - tet[1][2] * tet[2][0])
                                       + tet[0][2] * (tet[1][0] * tet[2][1] - tet[1][1] * tet[2][0]);
                                                                              
                     if (nMinus == 1) { 
                          return oneSixth * fabs(tripleProd);
                     } else {
                          return oneSixth * (1 - fabs(tripleProd));
                     }
                }
           }       
      }

	  // Note that area is an accumulation
	  double fluidAndAreaBox(double p000, double p100, double p010, double p110,
		                     double p001, double p101, double p011, double p111,
							 double & area) {
	      double vol = 0.0;
          // (0,0,0), (0,1,0), (1,1,0), (0,0,1)
		  vol += fluidTetrahedron(p000, p010, p110, p001, area);

          // (0,1,1), (0,1,0), (1,1,0), (0,0,1)
          // Mapping it back onto the first tetrahedron, we see that
          // the z coordinate is really the y coordinate of the original (negative sign),
          // the x coordinate is really the x coordiante of the original (positive sign),
          // the y coordinate is really the z coordinate of the original (negative sign)                   
          vol += fluidTetrahedron(p011, p010, p110, p001, area);
		  // (1,1,1), (0,1,1), (0,0,1), (1,1,0)
		  vol += fluidTetrahedron(p111, p011, p001, p110, area);
		  // (1,1,1), (1,0,1), (0,0,1), (1,1,0)
		  vol += fluidTetrahedron(p111, p101, p001, p110, area);              
          // (1,0,0), (1,0,1), (0,0,1), (1,1,0)
		  vol += fluidTetrahedron(p100, p101, p001, p110, area);
          // (0,0,0), (1,0,0), (1,1,0), (0,0,1)
          vol += fluidTetrahedron(p000, p100, p110, p001, area);		

		  return vol;
	  }
                        
      double computeVActs(); // Returns change in vActs.
	  double computeAllVActs();

	  void computeVActsExactMulti(); // Also computes vActs, but computes the exact triangularization and finds its volume.
	  void computeSingleExactVAct(int mat); // Computes for specified material, and tries to be somewhat efficient at least.
      
      static bool isNaN(double x) { return !(x > 0.0 || x < 1.0); }      
      
      // Checks against exact answer
      void checkSphereR2Normal() {
           double x = sqrt(1.0/3.0), y = sqrt(1.0/3.0), z = sqrt(1.0/3.0);
           double gridCx = (m-1)/2.0, gridCy = (n-1)/2.0, gridCz = (p-1)/2.0;
           
           int i = (int)(2*x/hx + gridCx), j = (int)(2*y/hy + gridCy), k = (int)(2*z/hz + gridCz);
           if (i > m || j > n || k > p) { return; }
           
           cerr << "Grid point " << i << " " << j << " " << k << endl;
           
           double n1, n2, n3;
           normalToInterface(0, i, j, k, n1, n2, n3);
           
           cerr << "Normal obtained (" << n1 << "," << n2 << "," << n3 << "); Desired (" << x << "," << y << "," << z << ")\n";
           
           double error = sqrt( (n1-x)*(n1-x) + (n2-y)*(n2-y) + (n3-z)*(n3-z) );
           
           cerr << "Normal error = " << error << endl;
      }

	  void saveInterface();
	  double interfaceChange();
};
