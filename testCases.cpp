#include "stdafx.h"
#include "Interface3D.h"

using namespace std;

void Interface3D::computeVolfracsFromPhi(char * exactFilename) {
	// Now compute volume
	for (Material = 0; Material < numPhases; Material++) {
		for (Index3 inds(1, 1, 1, coarseM, coarseN, coarseP, coarseM + 1, coarseN + 1, coarseP + 1); inds.valid(); ++inds) {
			// Init to nonzero / nonone value so NB starts at max
			grid(Material, inds) = 0.5;
		}
		updateCoarseNarrowBand();
	}

	if (exactFilename != __nullptr) {
		outputInterfaceMesh(exactFilename);
	}
	computeVActsExactMulti();

	for (Material = 0; Material < numPhases; Material++) {
		for (Index3 inds(1, 1, 1, coarseM, coarseN, coarseP, coarseM + 1, coarseN + 1, coarseP + 1); inds.valid(); ++inds) {
			//cerr << "Computed volume fraction " << grid[I][J][K] << endl;
			grid(Material, inds) = vActs(Material, inds);
			vActs(Material, inds) = 0.0;
		}
		updateCoarseNarrowBand();
	}

	//
	//adjustVolumes(grid);
}

void Interface3D::initCube(double length, double x, double y, double z,
	double alpha, double beta, double gamma, double graphRange)
{
	hx = hy = hz = graphRange / m;
	double h = hx;

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = -graphRange / 2 + h * (i - 0.5),
					yp = -graphRange / 2 + h * (j - 0.5),
					zp = -graphRange / 2 + h * (k - 0.5);

				// Rotate!
				double xt, yt, zt;

				// First in alpha (z)
				xt = xp * cos(alpha) - yp * sin(alpha);
				yt = xp * sin(alpha) + yp * cos(alpha);
				xp = xt; yp = yt;
				// Next in beta (y)
				xt = xp * cos(beta) - zp * sin(beta);
				zt = xp * sin(beta) + zp * cos(beta);
				xp = xt; zp = zt;
				// Finally in gamma (x)
				yt = yp * cos(gamma) - zp * sin(gamma);
				zt = yp * sin(gamma) + zp * cos(gamma);
				yp = yt; zp = zt;

				// Translate.
				xp += x;  yp += y; zp += z;

				//cerr << "xp yp zp = " << xp << " " << yp << " " << zp << endl;

				// Finally, compute the function
				interfaceGrid(0, inds) = max(-length / 2.0 + xp,
					max(-length / 2.0 - xp,
						max(-length / 2.0 + yp,
							max(-length / 2.0 - yp,
								max(-length / 2.0 + zp,
									-length / 2.0 - zp)))));

				if (numPhases > 1) {
					interfaceGrid(1, inds) = -interfaceGrid(0, inds);
				}
				inds++;
			}
		}
	}

	computeVolfracsFromPhi();
}

void Interface3D::initTetra(double x, double y, double z, double alpha, double beta, double gamma)
{
	hx = hy = hz = 2.0 / m;
	double h = hx;

	if (numPhases < 3) {
		cerr << "Boom\n";
		exit(1);
	}

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = -0.5 + h * i,
					yp = -0.5 + h * j,
					zp = -0.5 + h * k;

				// Rotate!
				double xt, yt, zt;

				// First in alpha (z)
				xt = xp * cos(alpha) - yp * sin(alpha);
				yt = xp * sin(alpha) + yp * cos(alpha);
				xp = xt; yp = yt;
				// Next in beta (y)
				xt = xp * cos(beta) - zp * sin(beta);
				zt = xp * sin(beta) + zp * cos(beta);
				xp = xt; zp = zt;
				// Finally in gamma (x)
				yt = yp * cos(gamma) - zp * sin(gamma);
				zt = yp * sin(gamma) + zp * cos(gamma);
				yp = yt; zp = zt;

				// Translate.
				xp += x;  yp += y; zp += z;

				//cerr << "xp yp zp = " << xp << " " << yp << " " << zp << endl;

				// Finally, compute the function
				interfaceGrid(0, inds) = -max(-zp,
					max(-xp,
						max(xp - yp,
							zp + yp - 1)));
				interfaceGrid(1, inds) = max(-zp,
					max(-xp,
						max(xp - yp,
							max(zp + yp - 1,
								max(0.5 * xp - yp - zp + 0.6, // base
									-0.08 * xp + 0.48 * yp + 0.8 * zp - 0.48)))));
				interfaceGrid(2, inds) = max(-interfaceGrid(2, inds), max(-interfaceGrid(0, inds),
					0.5 * xp - yp - zp + 0.6));
				inds++;
			}
		}
	}

	outputInterface("Tetra.vti");

	computeVolfracsFromPhi();
}

void Interface3D::initSphere(double radius, double x, double y, double z, double graphRange)
{
	hx = hy = hz = graphRange / m;
	double h = hx;

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = -graphRange / 2 + h * (i - 0.5),
					yp = -graphRange / 2 + h * (j - 0.5),
					zp = -graphRange / 2 + h * (k - 0.5);

				// Translate.
				xp += x;  yp += y; zp += z;

				// Put a sphere in there
				interfaceGrid(0, inds) = sqrt(xp*xp + yp * yp + zp * zp) - radius;
				if (numPhases > 1) {
					interfaceGrid(1, inds) = -sqrt(xp*xp + yp * yp + zp * zp) + radius;
				}
				inds++;
			}
		}
	}

	computeVolfracsFromPhi();
}

void Interface3D::initDeltoidish(double radius, double x, double y, double z,
	double alpha, double beta, double gamma, double graphRange)
{
	hx = hy = hz = graphRange / m;
	double h = hx;

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = -graphRange / 2 + h * (i - 0.5),
					yp = -graphRange / 2 + h * (j - 0.5),
					zp = -graphRange / 2 + h * (k - 0.5);

				// Rotate!
				double xt, yt, zt;

				// First in alpha (z)
				xt = xp * cos(alpha) - yp * sin(alpha);
				yt = xp * sin(alpha) + yp * cos(alpha);
				xp = xt; yp = yt;
				// Next in beta (y)
				xt = xp * cos(beta) - zp * sin(beta);
				zt = xp * sin(beta) + zp * cos(beta);
				xp = xt; zp = zt;
				// Finally in gamma (x)
				yt = yp * cos(gamma) - zp * sin(gamma);
				zt = yp * sin(gamma) + zp * cos(gamma);
				yp = yt; zp = zt;

				// Translate.
				xp += x;  yp += y; zp += z;

				//cerr << "xp yp zp = " << xp << " " << yp << " " << zp << endl;

				// Finally, compute the function
				double r2 = xp * xp + yp * yp + zp * zp, a = 0.5;
				interfaceGrid(0, inds) = r2 * r2 + 18.0 * a * a * r2 - 27 * a * a * a * a -
					8.0 * a * (xp*xp*xp + yp * yp*yp - 3 * xp * zp * zp /*- 3 * yp * zp * zp
							   - 3 * xp * zp * zp*/) - 1.0;
				inds++;
			}
		}
	}

	computeVolfracsFromPhi();
}

void Interface3D::initCone(double radius, double height, double x, double y, double z,
	double alpha, double beta, double gamma, double graphRange) {
	hx = hy = hz = graphRange / m;
	double h = hx;

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = -graphRange / 2 + h * (i - 0.5),
					yp = -graphRange / 2 + h * (j - 0.5),
					zp = -graphRange / 2 + h * (k - 0.5);

				// Rotate!
				double xt, yt, zt;

				// First in alpha (z)
				xt = xp * cos(alpha) - yp * sin(alpha);
				yt = xp * sin(alpha) + yp * cos(alpha);
				xp = xt; yp = yt;
				// Next in beta (y)
				xt = xp * cos(beta) - zp * sin(beta);
				zt = xp * sin(beta) + zp * cos(beta);
				xp = xt; zp = zt;
				// Finally in gamma (x)
				yt = yp * cos(gamma) - zp * sin(gamma);
				zt = yp * sin(gamma) + zp * cos(gamma);
				yp = yt; zp = zt;

				// Translate.
				xp += x;  yp += y; zp += z;

				//cerr << "xp yp zp = " << xp << " " << yp << " " << zp << endl;

				// Finally, compute the function
				double r2xy = xp * xp + yp * yp, a2 = height * height / (radius*radius);
				interfaceGrid(0, inds) = max(a2 * r2xy - (zp*zp + height * zp + height * height / 4.0),
					max(zp - height / 2.0, -zp - height / 2.0));
				inds++;
			}
		}
	}

	computeVolfracsFromPhi();
}

Interface3D * Interface3D::SFBay(int ns) {

	// each (x,y) cell will be 100 x 100 m
	double scale = 100.0;
	Interface3D * theBay = new Interface3D(1, 250, 250, 10, ns);
	theBay->hx = theBay->hy = 100.0 / (ns * scale);
	double maxDepth = 15000.0; // in cm
	theBay->hz = maxDepth / (10 * ns * 100.0 * scale); // maxDepth / 100 (to meters) / (10 * ns) # of grid points

	ifstream infile("sfbay.txt");

	int rows, cols;
	infile >> rows >> cols;

	int nread = 0;

	double mymax = 0.0;
	double cellHeight = maxDepth / theBay->coarseP;

	for (int I = 0; I < rows; I++) {
		if (I > 550) break;
		for (int J = 0; J < cols; J++) {
			int depth;
			infile >> depth;
			if (J <= 250 && I >= 300) {
				for (int K = 0; K < theBay->coarseP + 2; K++) {
					if (depth > mymax) mymax = depth;
					double dist = -((K - 1) * cellHeight - depth) / cellHeight,
						vol;
					if (dist < 0) {
						vol = 0.0;
					}
					else if (dist < 1) {
						vol = dist;
					}
					else vol = 1.0;

					//cerr << vol << " ";

					theBay->grid(0, 1 + I - 300, 1 + J, K) = vol;

					/*int iStart = (I-300)*ns, jStart = J * ns, kStart = (K-1)*ns;
					if (K < 1 || K == theBay->coarseP+1 ||
						J >= 250 || (I-300)>=250) { continue; }
					//cerr << iStart << " " << jStart << " " << kStart << endl;
					for (int i = 0; i < ns; i++) {
						for (int j = 0; j < ns; j++) {
							for (int k = 0; k < ns; k++) {
								theBay->interfaceGrid[iStart+i][jStart+j][kStart+k] =
										(double) k / theBay->ns - vol;
							}
						}
					} // test*/
				}
			}
		}
	}

	Grid3<double> tmp(250, 250, 10);
	for (int I = 1; I <= 250; I++) {
		for (int J = 1; J <= 250; J++) {
			for (int K = 1; K <= 10; K++) {
				double sum = 0.0;
				for (int Ip = -1; Ip <= 1; Ip++) {
					for (int Jp = -1; Jp <= 1; Jp++) {
						for (int Kp = -1; Kp <= 1; Kp++) {
							if (Ip == 0 && Jp == 0 && Kp == 0) continue;
							sum += theBay->grid(0, I + Ip, J + Jp, K + Kp);
						}
					}
				}
				tmp(I - 1, J - 1, K - 1) = theBay->grid(0, I, J, K) * 0.5 + 0.5 * sum / 26.0;
			}
		}
	}
	for (int I = 1; I <= 250; I++) {
		for (int J = 1; J <= 250; J++) {
			for (int K = 1; K <= 10; K++) {
				theBay->grid(0, I, J, K) = tmp(I - 1, J - 1, K - 1);
			}
		}
	}

	theBay->updateCoarseNarrowBand();

	cerr << "Max depth found was " << mymax << endl;
	cerr << "Read input\n";
	system("PAUSE");

	return theBay;
}

void Interface3D::adjustVolumes(MultiGrid3<double> & whichGrid) {
	for (Index3 inds(1, 1, 1, coarseM, coarseN, coarseP, coarseM + 1, coarseN + 1, coarseP + 1); inds.valid(); ++inds) {
		double totalVol = 0.0;
		for (int mat = 0; mat < numPhases; mat++) {
			totalVol += whichGrid(mat, inds);
		}

		double error = 1.0 - totalVol;

		for (int mat = 0; mat < numPhases; mat++) {
			whichGrid(mat, inds) *= (1.0 + error / totalVol);
		}
	}
}


void Interface3D::initTripleLine(double x, double y, double z, double alpha, double beta, double gamma, bool inSphere) {
	hx = hy = hz = 2.0 / m;
	double h = hx;
	const double graphRange = 2.0;
	if (numPhases < 3 || (numPhases < 4 && inSphere)) {
		cerr << "Error, cannot initialize a triple point without at least three phases\n";
		exit(1);
	}

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = -graphRange / 2 + h * i,
					yp = -graphRange / 2 + h * j,
					zp = -graphRange / 2 + h * k;

				// Rotate!
				double xt, yt, zt;

				// Let's do a little deformation
				double alphaStar = alpha + 0.25 * zp;// + 0.5 * sqrt(xp*xp + yp*yp);

				// First in alpha (z)
				xt = xp * cos(alphaStar) - yp * sin(alphaStar);
				yt = xp * sin(alphaStar) + yp * cos(alphaStar);
				xp = xt; yp = yt;
				// Next in beta (y)
				xt = xp * cos(beta) - zp * sin(beta);
				zt = xp * sin(beta) + zp * cos(beta);
				xp = xt; zp = zt;
				// Finally in gamma (x)
				yt = yp * cos(gamma) - zp * sin(gamma);
				zt = yp * sin(gamma) + zp * cos(gamma);
				yp = yt; zp = zt;

				// Translate.
				xp += x;  yp += y; zp += z;

				double val1 = sqrt(xp*xp + (yp - 0.5)*(yp - 0.5) + zp * zp) - 0.6,
					val2 = sqrt((xp - sqrt(3.0) / 4.0)*(xp - sqrt(3.0) / 4.0) + (yp + 0.25)*(yp + 0.25) + zp * zp) - 0.6,
					val3 = sqrt((xp + sqrt(3.0) / 4.0)*(xp + sqrt(3.0) / 4.0) + (yp + 0.25)*(yp + 0.25) + zp * zp) - 0.6;

				double //val1p = 0.5 * (val1 - min(val2, val3)),
					   //val2p = 0.5 * (val2 - min(val1, val3)),
					   //val3p = 0.5 * (val3 - min(val1, val2)),
					valSphere = 0.75 - sqrt(xp*xp + yp * yp + zp * zp); // Magic number for now. Can add variable later.

				if (inSphere) {
					interfaceGrid(0, inds) = val1,//max(val1p, -valSphere);
						interfaceGrid(1, inds) = val2,//max(val2p, -valSphere);
						interfaceGrid(2, inds) = val3,//max(val3p, -valSphere);
						interfaceGrid(3, inds) = valSphere;
				}
				else {
					interfaceGrid(0, inds) = val1;
					interfaceGrid(1, inds) = val2;
					interfaceGrid(2, inds) = val3;
				}
				inds++;
			}
		}
	}

	computeVolfracsFromPhi("TripleLineMeshInit.vtp");

	//outputInterface("TripleLineTest.vti");
}

void Interface3D::initTripleLineFlat(double x, double y, double z, double alpha, double beta, double gamma) {
	hx = hy = hz = 2.0 / m;
	double h = hx;
	const double graphRange = 2.0;
	if (numPhases < 3) {
		cerr << "Error, cannot initialize a triple point without at least three phases\n";
		exit(1);
	}

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = -graphRange / 2 + h * i,
					yp = -graphRange / 2 + h * j,
					zp = -graphRange / 2 + h * k;

				// Rotate!
				double xt, yt, zt;

				// Let's do a little deformation
				double alphaStar = alpha;// + 0.25 * zp;// + 0.5 * sqrt(xp*xp + yp*yp);

				// First in alpha (z)
				xt = xp * cos(alphaStar) - yp * sin(alphaStar);
				yt = xp * sin(alphaStar) + yp * cos(alphaStar);
				xp = xt; yp = yt;
				// Next in beta (y)
				xt = xp * cos(beta) - zp * sin(beta);
				zt = xp * sin(beta) + zp * cos(beta);
				xp = xt; zp = zt;
				// Finally in gamma (x)
				yt = yp * cos(gamma) - zp * sin(gamma);
				zt = yp * sin(gamma) + zp * cos(gamma);
				yp = yt; zp = zt;

				// Translate.
				xp += x;  yp += y; zp += z;

				double val1 = yp, val2 = -yp, val3 = 100 * zp;

				interfaceGrid(0, inds) = val1;
				interfaceGrid(1, inds) = val2;
				interfaceGrid(2, inds) = val3;

				inds++;
			}
		}
	}

	computeVolfracsFromPhi("FlatTripleLineMeshInit.vtp");

	//outputInterface("TripleLineTest.vti");
}

void Interface3D::initConcentricCircles() {
	// Radii of spheres: 6/13, 4.75/13, 3.5/13, 2.25/13, 1/13
	// Create volume fractions as a set of annuli	
	hx = hy = hz = 1.0 / m;
	double h = hx;
	const double graphRange = 1.0;
	if (numPhases < 6) {
		cerr << "Error, need six phases for this Concentric Spheres test\n";
		exit(1);
	}

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = h * i,
					yp = h * j,
					zp = h * k;
				// Radii of spheres: 6/13, 4.75/13, 3.5/13, 2.25/13, 1/13                      
				double sphere = sqrt((xp - 0.5)*(xp - 0.5) + (yp - 0.5)*(yp - 0.5) + (zp - 0.5)*(zp - 0.5)),
					sphere1 = sphere - 1.0 / 13.0,
					sphere2 = sphere - 2.25 / 13.0,
					sphere3 = sphere - 3.5 / 13.0,
					sphere4 = sphere - 4.75 / 13.0,
					sphere5 = sphere - 6.0 / 13.0;

				double val1 = sphere1,
					val2 = max(-sphere1, sphere2),
					val3 = max(-sphere2, sphere3),
					val4 = max(-sphere3, sphere4),
					val5 = max(-sphere4, sphere5),
					val6 = -sphere5;

				// Fix overlap (just in case)
				double val1p = 0.5 * (val1 - min(val2, min(val3, min(val4, min(val5, val6))))),
					val2p = 0.5 * (val2 - min(val1, min(val3, min(val4, min(val5, val6))))),
					val3p = 0.5 * (val3 - min(val1, min(val2, min(val4, min(val5, val6))))),
					val4p = 0.5 * (val4 - min(val1, min(val2, min(val3, min(val5, val6))))),
					val5p = 0.5 * (val5 - min(val1, min(val2, min(val3, min(val4, val6))))),
					val6p = 0.5 * (val6 - min(val1, min(val2, min(val3, min(val4, val5)))));

				interfaceGrid(0, inds) = val1p;
				interfaceGrid(1, inds) = val2p;
				interfaceGrid(2, inds) = val3p;
				interfaceGrid(3, inds) = val4p;
				interfaceGrid(4, inds) = val5p;
				interfaceGrid(5, inds) = val6p;
				inds++;
			}
		}
	}

	computeVolfracsFromPhi("ConcentricSpheresTest.vti");
}

void Interface3D::initCubeSphereIntersect(double alpha, double beta, double gamma) {
	// Material 1 = Sphere, Center (0.6, 0.6, 0.6) R 0.2
	// Material 2 = Box, (0.2, 0.2, 0.2) to (0.6, 0.6, 0.6), subtracted out the parts in the sphere
	// Material 3 = Outside

	// Sphere volfracs created easily
	// Box: Max(Cube, -Sphere)
	// Outside: Min(Cube, Sphere)
	hx = hy = hz = 1.0 / m;
	double h = hx;
	const double graphRange = 1.0;
	if (numPhases < 3) {
		cerr << "Error, need three phases for this Sphere/Box test\n";
		exit(1);
	}

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = h * i,
					yp = h * j,
					zp = h * k;

				// Rotate!
				double xt, yt, zt;

				// Let's do a little deformation
				double alphaStar = alpha;// + 0.25 * zp;// + 0.5 * sqrt(xp*xp + yp*yp);

				// First in alpha (z)
				xt = xp * cos(alphaStar) - yp * sin(alphaStar);
				yt = xp * sin(alphaStar) + yp * cos(alphaStar);
				xp = xt; yp = yt;
				// Next in beta (y)
				xt = xp * cos(beta) - zp * sin(beta);
				zt = xp * sin(beta) + zp * cos(beta);
				xp = xt; zp = zt;
				// Finally in gamma (x)
				yt = yp * cos(gamma) - zp * sin(gamma);
				zt = yp * sin(gamma) + zp * cos(gamma);
				yp = yt; zp = zt;

				// val1 = sphere
				// val2 = 
				double sphere = sqrt((xp - 0.6)*(xp - 0.6) + (yp - 0.6)*(yp - 0.6) + (zp - 0.6)*(zp - 0.6)) - 0.2,
					box = max(-(xp - 0.2),
						max(xp - 0.6,
							max(-(yp - 0.2),
								max(yp - 0.6,
									max(-(zp - 0.2),
										zp - 0.6)))));

				double val1 = sphere, val2 = max(box, -sphere), val3 = -min(box, sphere);

				// Fix overlap (just in case)
				double val1p = 0.5 * (val1 - min(val2, val3)),
					val2p = 0.5 * (val2 - min(val1, val3)),
					val3p = 0.5 * (val3 - min(val1, val2));

				interfaceGrid(0, inds) = val1; //p
				interfaceGrid(1, inds) = val2;
				interfaceGrid(2, inds) = val3;
				inds++;
			}
		}
	}

	computeVolfracsFromPhi("BoxSphereTest.vtp");
}

void Interface3D::initQuadruplePoint(double x, double y, double z, double alpha, double beta, double gamma) {
	const double graphRange = 1.0;
	hx = hy = hz = graphRange / m;
	double h = hx;

	if (numPhases < 4) {
		cerr << "Error, cannot initialize a quadruple point without at least four phases\n";
		exit(1);
	}

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = h * i,
					yp = h * j,
					zp = h * k;

				// Rotate!
				double xt, yt, zt;

				// Let's do a little deformation
				double alphaStar = alpha + 0.25 * (zp - 0.5) + 0.5 * sqrt((xp - 0.5)*(xp - 0.5) + (yp - 0.5)*(yp - 0.5));

				// First in alpha (z)
				xt = xp * cos(alphaStar) - yp * sin(alphaStar);
				yt = xp * sin(alphaStar) + yp * cos(alphaStar);
				xp = xt; yp = yt;
				// Next in beta (y)
				xt = xp * cos(beta) - zp * sin(beta);
				zt = xp * sin(beta) + zp * cos(beta);
				xp = xt; zp = zt;
				// Finally in gamma (x)
				yt = yp * cos(gamma) - zp * sin(gamma);
				zt = yp * sin(gamma) + zp * cos(gamma);
				yp = yt; zp = zt;

				// Translate.
				xp += x;  yp += y; zp += z;

				double dist1 = sqrt((xp - 0.0)*(xp - 0.0) + (yp - 0.0)*(yp - 0.0) + (zp - 0.0)*(zp - 0.0)),
					dist2 = sqrt((xp - 1.0)*(xp - 1.0) + (yp - 1.0)*(yp - 1.0) + (zp - 0.0)*(zp - 0.0)),
					dist3 = sqrt((xp - 0.0)*(xp - 0.0) + (yp - 1.0)*(yp - 1.0) + (zp - 1.0)*(zp - 1.0)),
					dist4 = sqrt((xp - 1.0)*(xp - 1.0) + (yp - 0.0)*(yp - 0.0) + (zp - 1.0)*(zp - 1.0));

				/*double val1 = min(dist1 - dist2, min(dist1 - dist3, dist1 - dist4)),
					   val2 = min(dist2 - dist1, min(dist2 - dist3, dist2 - dist4)),
					   val3 = min(dist3 - dist1, min(dist3 - dist2, dist3 - dist4)),
					   val4 = min(dist4 - dist1, min(dist4 - dist2, dist4 - dist3));

				double val1p = 0.5 * (val1 - min(val2, min(val3, val4))),
					   val2p = 0.5 * (val2 - min(val1, min(val3, val4))),
					   val3p = 0.5 * (val3 - min(val1, min(val2, val4))),
					   val4p = 0.5 * (val4 - min(val1, min(val2, val3)));*/

				interfaceGrid(0, inds) = dist1;
				interfaceGrid(1, inds) = dist2;
				interfaceGrid(2, inds) = dist3;
				interfaceGrid(3, inds) = dist4;

				inds++;
			}
		}
	}

	computeVolfracsFromPhi("QuadruplePointTest.vtp");
}

void Interface3D::initRandomVoronoi(int seed, bool swirl) {
	srand(seed);
	hx = hy = hz = 6.0 / m;

	for (int ind = 0; ind < fineSize; ind++) {
		for (int mat = 0; mat < numPhases; mat++) {
			interfaceGrid(mat, ind) = INFTY;
			narrowBand(mat, ind) = true;
		}
	}

	int maxInt = RAND_MAX;
	for (Material = 0; Material < numPhases; Material++) {
		double x = ((double)rand()) / maxInt * m,
			y = ((double)rand()) / maxInt * n,
			z = ((double)rand()) / maxInt * p;

		//cerr << "Points: " << x << " " << y << " " << z << endl;

		int ind = 0;
		for (int i = 0; i <= m; i++) {
			for (int j = 0; j <= n; j++) {
				for (int k = 0; k <= p; k++) {
					double coordX = (double)i, coordY = (double)j, coordZ = (double)k;
					if (swirl) {
						double cx = m / 2, cy = n / 2, cz = p / 2, r = sqrt((coordX - cx)*(coordX - cx) + (coordY - cy)*(coordY - cy) + (coordZ - cz)*(coordZ - cz));
						double alpha = M_PI * r * exp(-r / m) / m, gamma = M_PI * r * exp(-1.2*r / m) / m;
						double xt, yt, zt;

						xt = cos(alpha) * (coordX - cx) - sin(alpha) * (coordY - cy);
						yt = sin(alpha) * (coordX - cx) + cos(alpha) * (coordY - cy);
						coordX = xt + cx; coordY = yt + cy;

						xt = cos(gamma) * (coordX - cx) - sin(gamma) * (coordZ - cz);
						zt = sin(gamma) * (coordX - cx) + cos(gamma) * (coordZ - cz);
						coordX = xt + cx; coordZ = zt + cz;
					}
					double dx = x - coordX, dy = y - coordY, dz = z - coordZ;
					double dist = sqrt(dx*dx + dy * dy + dz * dz);

					interfaceGrid(Material, ind) = dist / ns;
					ind++;
					//cerr << "Inserting " << dist << " at " << ip << " " << jp << " " << kp << endl;
				}
			}
		}
	}

	//fixLevelsetOverlapNoNB();

	//reinitialize();

	computeVolfracsFromPhi("PerfectVoronoiMesh.vtp");
	//outputInterface("RandomVoronoiTest.vti");
}

void Interface3D::initHemispheres(double radius, double x, double y, double z,
	double alpha, double beta, double gamma, double graphRange)
{
	hx = hy = hz = graphRange / m;
	double h = hx;

	if (numPhases < 3) { exit(1); }

	// Initialize interfaceGrid (phi) to be the signed distance function to the cube.
	int inds = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double xp = -graphRange / 2 + h * (i - 0.5),
					yp = -graphRange / 2 + h * (j - 0.5),
					zp = -graphRange / 2 + h * (k - 0.5);

				// Rotate!
				double xt, yt, zt;

				// First in alpha (z)
				xt = xp * cos(alpha) - yp * sin(alpha);
				yt = xp * sin(alpha) + yp * cos(alpha);
				xp = xt; yp = yt;
				// Next in beta (y)
				xt = xp * cos(beta) - zp * sin(beta);
				zt = xp * sin(beta) + zp * cos(beta);
				xp = xt; zp = zt;
				// Finally in gamma (x)
				yt = yp * cos(gamma) - zp * sin(gamma);
				zt = yp * sin(gamma) + zp * cos(gamma);
				yp = yt; zp = zt;

				// Translate.
				xp += x;  yp += y; zp += z;

				//cerr << "xp yp zp = " << xp << " " << yp << " " << zp << endl;

				// Finally, compute the function
				double r = sqrt(xp*xp + yp * yp + zp * zp);
				double top = max(r - radius, zp), bottom = max(r - radius, -zp), outside = radius - r;

				interfaceGrid(0, inds) = top;
				interfaceGrid(1, inds) = bottom;
				interfaceGrid(2, inds) = outside;
				inds++;
			}
		}
	}

	computeVolfracsFromPhi("HemispheresTest.vtp");
}