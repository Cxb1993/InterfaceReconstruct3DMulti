#include "stdafx.h"
#include "Interface3D.h"

using namespace std;

const int Interface3D::oppositeEdge[6] = { 5, 4, 3, 2, 1, 0 };
const double Interface3D::vertices[4][3] = { {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };
const int Interface3D::edgesTouchingVertices[4][3] = { {0,1,2}, {0,3,4}, {1,3,5}, {2,4,5} };
int ** Interface3D::unitTetFaces = new int*[4];

// Updated.
Interface3D::Interface3D(int nPhases, int gm, int gn, int gp, int _ns,
	double _h, ifstream & input) : numPhases(nPhases),
	coarseM(gm), m(gm * _ns), coarseN(gn), n(gn * _ns), coarseP(gp), p(gp * _ns), ns(_ns), hx(_h),
	hy(_h), hz(_h), coarseSize((gm + 2) * (gn + 2) * (gp + 2)), fineSize((gm*_ns + 1) * (gn*_ns + 1) * (gp*_ns + 1)),
	grid(nPhases, gm + 1, gn + 1, gp + 1),
	vActs(nPhases, gm + 1, gn + 1, gp + 1),
	surfaceAreas(nPhases, gm + 1, gn + 1, gp + 1),
	coarseNarrowBand(nPhases, gm + 1, gn + 1, gp + 1),
	vActsBack(nPhases, gm + 1, gn + 1, gp + 1),
	lambdas(gm + 1, gn + 1, gp + 1),
	residual(gm + 1, gn + 1, gp + 1),
	gradientBand(gm + 1, gn + 1, gp + 1),
	curvatureSums(gm + 1, gn + 1, gp + 1),
	// Fine initializations
	interfaceGrid(nPhases, gm*_ns, gn*_ns, gp*_ns),
	velocities(gm*_ns, gn*_ns, gp*_ns),
	workGrid(gm*_ns, gn*_ns, gp*_ns),
	oldInterface(nPhases, gm*_ns, gn*_ns, gp*_ns),
	signs(gm*_ns, gn*_ns, gp*_ns),
	narrowBand(nPhases, gm*_ns, gn*_ns, gp*_ns),
	volumeGradients(gm*_ns, gn*_ns, gp*_ns),
	volumeOfGradients(gm + 1, gn + 1, gp + 1)
{
	allocate();

	for (Material = 0; Material < numPhases; Material++) {
		for (int k = 1; k < coarseP + 1; k++) {
			for (int i = 1; i < coarseM + 1; i++) {
				for (int j = 1; j < coarseN + 1; j++) {
					input >> grid(Material, i, j, k); // Read in volume fractions
				}
			}
		}
	}
}

Interface3D::Interface3D(int nPhases, int gm, int gn, int gp, int _ns) :
	numPhases(nPhases), coarseM(gm), m(gm * _ns), coarseN(gn), n(gn * _ns), coarseP(gp), p(gp * _ns), ns(_ns),
	coarseSize((gm + 2) * (gn + 2) * (gp + 2)), fineSize((gm*_ns + 1) * (gn*_ns + 1) * (gp*_ns + 1)),
	grid(nPhases, gm + 1, gn + 1, gp + 1),
	vActs(nPhases, gm + 1, gn + 1, gp + 1),
	surfaceAreas(nPhases, gm + 1, gn + 1, gp + 1),
	coarseNarrowBand(nPhases, gm + 1, gn + 1, gp + 1),
	vActsBack(nPhases, gm + 1, gn + 1, gp + 1),
	lambdas(gm + 1, gn + 1, gp + 1),
	residual(gm + 1, gn + 1, gp + 1),
	gradientBand(gm + 1, gn + 1, gp + 1),
	curvatureSums(gm + 1, gn + 1, gp + 1),
	// Fine initializations
	interfaceGrid(nPhases, gm*_ns, gn*_ns, gp*_ns),
	velocities(gm*_ns, gn*_ns, gp*_ns),
	workGrid(gm*_ns, gn*_ns, gp*_ns),
	oldInterface(nPhases, gm*_ns, gn*_ns, gp*_ns),
	signs(gm*_ns, gn*_ns, gp*_ns),
	narrowBand(nPhases, gm*_ns, gn*_ns, gp*_ns),
	volumeGradients(gm*_ns, gn*_ns, gp*_ns),
	volumeOfGradients(gm + 1, gn + 1, gp + 1),
	hx(1.0 / _ns), hy(1.0 / _ns), hz(1.0 / _ns) {
	allocate();
}

// Allocates all data, given that we know how big to make this thing.
void Interface3D::allocate()
{
	// Initialize volume fraction array
	// Empty front & back face
	for (Material = 0; Material < numPhases; Material++) {
		for (int i = 0; i < coarseM + 2; i++) {
			for (int j = 0; j < coarseN + 2; j++) {
				grid(Material, i, j, 0) = 0.0;
				grid(Material, i, j, coarseP + 1) = 0.0;

				vActs(Material, i, j, 0) = 0.0;
				vActs(Material, i, j, coarseP + 1) = 0.0;
			}
		}


		for (int k = 1; k < coarseP + 1; k++) {
			for (int j = 0; j < coarseN + 2; j++) {
				grid(Material, 0, j, k) = 0.0; // empty left
				grid(Material, coarseM + 1, j, k) = 0.0; // empty right

				vActs(Material, 0, j, k) = 0.0; // empty left
				vActs(Material, coarseM + 1, j, k) = 0.0; // empty right               
			}
			for (int i = 1; i < coarseM + 1; i++) {
				grid(Material, i, 0, k) = 0.0; // Empty top
				vActs(Material, i, 0, k) = 0.0;

				grid(Material, i, coarseN + 1, k) = 0.0; // Empty bottom
				vActs(Material, i, coarseN + 1, k) = 0.0;
			}
		}
	}

	// Initialize volume fraction array
	//setupGhostFractions();
	// Seems reasonable to just allocate to max capacity (but only use a small part of it)
	narrowBandSize = new int[numPhases];
	narrowBandIterator = new Node3D*[numPhases];
	for (int mat = 0; mat < numPhases; mat++) {
		narrowBandIterator[mat] = new Node3D[(m + 1)*(n + 1)*(p + 1)];
		narrowBandSize[mat] = 0;
	}

	for (Index3 ind(m, n, p); ind.valid(); ++ind) {
		volumeGradients(ind).setXYZ(ind.i, ind.j, ind.k);
	}

	unitTetFaces[0] = new int[5]; unitTetFaces[1] = new int[5]; unitTetFaces[2] = new int[5]; unitTetFaces[3] = new int[5];
	unitTetFaces[0][0] = 0; unitTetFaces[0][1] = 1; unitTetFaces[0][2] = 2; unitTetFaces[0][3] = 0; unitTetFaces[0][4] = -1;
	unitTetFaces[1][0] = 0; unitTetFaces[1][1] = 3; unitTetFaces[1][2] = 1; unitTetFaces[1][3] = 0; unitTetFaces[1][4] = -1;
	unitTetFaces[2][0] = 0; unitTetFaces[2][1] = 2; unitTetFaces[2][2] = 3; unitTetFaces[2][3] = 0; unitTetFaces[2][4] = -1;
	unitTetFaces[3][0] = 1; unitTetFaces[3][1] = 3; unitTetFaces[3][2] = 2; unitTetFaces[3][3] = 1; unitTetFaces[3][4] = -1;

	const int right = interfaceGrid.right, top = interfaceGrid.top, forward = interfaceGrid.forward;
	// (0,0,0), (0,1,0), (1,1,0), (0,0,1)
	tetrahedraOffsets[0][0] = 0;
	tetrahedraOffsets[0][1] = top;
	tetrahedraOffsets[0][2] = right + top;
	tetrahedraOffsets[0][3] = forward;
	// (0,1,1), (0,1,0), (1,1,0), (0,0,1)
	tetrahedraOffsets[1][0] = top + forward;
	tetrahedraOffsets[1][1] = top;
	tetrahedraOffsets[1][2] = right + top;
	tetrahedraOffsets[1][3] = forward;
	// (1,1,1), (0,1,1), (0,0,1), (1,1,0)
	tetrahedraOffsets[2][0] = right + top + forward;
	tetrahedraOffsets[2][1] = top + forward;
	tetrahedraOffsets[2][2] = forward;
	tetrahedraOffsets[2][3] = right + top;
	// (1,1,1), (1,0,1), (0,0,1), (1,1,0)
	tetrahedraOffsets[3][0] = right + top + forward;
	tetrahedraOffsets[3][1] = right + forward;
	tetrahedraOffsets[3][2] = forward;
	tetrahedraOffsets[3][3] = right + top;
	// (1,0,0), (1,0,1), (0,0,1), (1,1,0)
	tetrahedraOffsets[4][0] = right;
	tetrahedraOffsets[4][1] = right + forward;
	tetrahedraOffsets[4][2] = forward;
	tetrahedraOffsets[4][3] = right + top;
	// (0,0,0), (1,0,0), (1,1,0), (0,0,1)
	tetrahedraOffsets[5][0] = 0;
	tetrahedraOffsets[5][1] = right;
	tetrahedraOffsets[5][2] = right + top;
	tetrahedraOffsets[5][3] = forward;

	iter = 0;
}
// Allocate

/*void Interface3D::setupGhostFractions() {
	   for (int mat = 0; mat < numPhases; mat++) {
		   for (int i = 1; i < coarseM+1; i++) {
			   for (int j = 1; j < coarseN+1; j++) {
				   grid(mat,i,j,0) = grid(mat,i,j,1); // copy front face
				   grid(mat,i,j,coarseP+1) = grid(mat,i,j,coarseP); // copy back face
			   }
			   for (int k = 1; k < coarseP+1; k++) {
				   grid[mat][i][0][k] = grid[mat][i][1][k]; // copy bottom face
				   grid[mat][i][coarseN+1][k] = grid[mat][i][coarseN][k]; // copy top face
			   }

			   // edges w/ variable i (4)
			   grid[mat][i][0][0] = grid[mat][i][1][1];
			   grid[mat][i][0][coarseP+1] = grid[mat][i][1][coarseP];
			   grid[mat][i][coarseN+1][0] = grid[mat][i][coarseN][1];
			   grid[mat][i][coarseN+1][coarseP+1] = grid[mat][i][coarseN][coarseP];
		   }

		   for (int j = 1; j < coarseN+1; j++) {
			   // edges w/ variable j (4)
			   grid[mat][0][j][0] = grid[mat][1][j][1];
			   grid[mat][0][j][coarseP+1] = grid[mat][1][j][coarseP];
			   grid[mat][coarseM+1][j][0] = grid[mat][coarseM][j][1];
			   grid[mat][coarseM+1][j][coarseP+1] = grid[mat][coarseM][j][coarseP];

			   for (int k = 1; k < coarseP+1; k++) {
				   grid[mat][0][j][k] = grid[mat][1][j][k]; // copy left face
				   grid[mat][coarseM+1][j][k] = grid[mat][coarseM][j][k]; // copy right face
			   }
		   }

		   for (int k = 1; k < coarseP+1; k++) {
			   // edges w/ variable k(4)
			   grid[mat][0][0][k] = grid[mat][1][1][k];
			   grid[mat][0][coarseN+1][k] = grid[mat][1][coarseN][k];
			   grid[mat][coarseM+1][0][k] = grid[mat][coarseM][1][k];
			   grid[mat][coarseM+1][coarseN+1][k] = grid[mat][coarseM][coarseN][k];
		   }
		   // Do 8 corners
		   grid[mat][0][0][0] = grid[mat][1][1][1];
		   grid[mat][coarseM+1][0][0] = grid[mat][coarseM][1][1];
		   grid[mat][0][coarseN+1][0] = grid[mat][1][coarseN][1];
		   grid[mat][coarseM+1][coarseN+1][0] = grid[mat][coarseM][coarseN][1];

		   grid[mat][0][0][coarseP+1] = grid[mat][1][1][coarseP];
		   grid[mat][coarseM+1][0][coarseP+1] = grid[mat][coarseM][1][coarseP];
		   grid[mat][0][coarseN+1][coarseP+1] = grid[mat][1][coarseN][coarseP];
		   grid[mat][coarseM+1][coarseN+1][coarseP+1] = grid[mat][coarseM][coarseN][coarseP];

		   for (int I = 0; I < coarseM+2; I++) {
			   for (int J = 0; J < coarseN+2; J++) {
				   for (int K = 0; K < coarseP+2; K++) {
					   surfaceAreas[mat][I][J][K] = 0.0;
					   for (int mat = 0; mat < numPhases; mat++) { vActs[mat][I][J][K] = 0.0; }
				   }
			   }
		   }
	   }
}*/

void Interface3D::updateGhostLayer() {
	for (int mat = 0; mat < numPhases; mat++) {
		for (int i = 1; i < m; i++) {
			for (int j = 1; j < n; j++) {
				interfaceGrid(mat, i, j, 0) = safeExtrapolate(mat, i, j, 1, i, j, 2);
				interfaceGrid(mat, i, j, p) = safeExtrapolate(mat, i, j, p - 1, i, j, p - 2);
			}
			for (int k = 1; k < p; k++) {
				interfaceGrid(mat, i, 0, k) = safeExtrapolate(mat, i, 1, k, i, 2, k);
				interfaceGrid(mat, i, n, k) = safeExtrapolate(mat, i, n - 1, k, i, n - 2, k);
			}

			// edges w/ variable i (4)
			interfaceGrid(mat, i, 0, 0) = safeExtrapolate(mat, i, 1, 0, i, 2, 0);
			interfaceGrid(mat, i, 0, p) = safeExtrapolate(mat, i, 1, p, i, 2, p);
			interfaceGrid(mat, i, n, 0) = safeExtrapolate(mat, i, n - 1, 0, i, n - 2, 0);
			interfaceGrid(mat, i, n, p) = safeExtrapolate(mat, i, n - 1, p, i, n - 2, p);
		}

		for (int j = 1; j < n; j++) {
			for (int k = 1; k < p; k++) {
				interfaceGrid(mat, 0, j, k) = safeExtrapolate(mat, 1, j, k, 2, j, k);
				interfaceGrid(mat, m, j, k) = safeExtrapolate(mat, m - 1, j, k, m - 2, j, k);
			}
			// edges w/ variable j (4)
			interfaceGrid(mat, 0, j, 0) = safeExtrapolate(mat, 1, j, 0, 2, j, 0);
			interfaceGrid(mat, 0, j, p) = safeExtrapolate(mat, 1, j, p, 2, j, p);
			interfaceGrid(mat, m, j, 0) = safeExtrapolate(mat, m - 1, j, 0, m - 2, j, 0);
			interfaceGrid(mat, m, j, p) = safeExtrapolate(mat, m - 1, j, p, m - 2, j, p);
		}

		for (int k = 1; k < p; k++) {
			// edges w/ variable k(4)
			interfaceGrid(mat, 0, 0, k) = safeExtrapolate(mat, 1, 0, k, 2, 0, k);
			interfaceGrid(mat, 0, n, k) = safeExtrapolate(mat, 1, n, k, 2, n, k);
			interfaceGrid(mat, m, 0, k) = safeExtrapolate(mat, m - 1, 0, k, m - 2, 0, k);
			interfaceGrid(mat, m, n, k) = safeExtrapolate(mat, m - 1, n, k, m - 2, n, k);
		}
		// Do 8 corners
		interfaceGrid(mat, 0, 0, 0) = safeExtrapolate(mat, 1, 0, 0, 2, 0, 0);
		interfaceGrid(mat, m, 0, 0) = safeExtrapolate(mat, m - 1, 0, 0, m - 2, 0, 0);
		interfaceGrid(mat, 0, n, 0) = safeExtrapolate(mat, 0, n - 1, 0, 0, n - 2, 0);
		interfaceGrid(mat, m, n, 0) = safeExtrapolate(mat, m - 1, n, 0, m - 2, n, 0);

		interfaceGrid(mat, 0, 0, p) = safeExtrapolate(mat, 1, 0, p, 2, 0, p);
		interfaceGrid(mat, m, 0, p) = safeExtrapolate(mat, m - 1, 0, p, m - 2, 0, p);
		interfaceGrid(mat, 0, n, p) = safeExtrapolate(mat, 0, n - 1, p, 0, n - 2, p);
		interfaceGrid(mat, m, n, p) = safeExtrapolate(mat, m - 1, n, p, m - 2, n, p);
	}
}

// Updates the listing of which coarse grid cells are cared about
// Depends on [material]
void Interface3D::updateCoarseNarrowBand() {
	for (Index3 inds(1, 1, 1, coarseM, coarseN, coarseP, coarseM + 1, coarseN + 1, coarseP + 1); inds.valid(); ++inds) {
		// Check neighbors
		bool closeToInterface = false;
		double vrCenter = grid(Material, inds);

		for (Index3 neighbors(inds.i - 1, inds.j - 1, inds.k - 1, inds.i + 1, inds.j + 1, inds.k + 1, coarseM + 1, coarseN + 1, coarseP + 1); neighbors.valid(); ++neighbors) {
			double vr = grid(Material, neighbors), va = vActs(Material, neighbors);
			if ((vr > 1e-8 && vr < 1 - 1e-8) || (va > 1e-8 && va < 1 - 1e-8) ||
				fabs(vr - vrCenter) > 1e-8) {
				closeToInterface = true;
				break;
			}
		}

		coarseNarrowBand(Material, inds) = closeToInterface;
	}
}

// Original seeding routine, and the current one too! Does a level set of the VoF data
void Interface3D::seedGrid(double levelSet) {

	// TEST THE MESHER
	/*std::vector<double> pp;
	std::vector<int> ss;
	std::vector<int> faceType;
	int indices[4][3], transform[3];
	transform[0] = 1; transform[1] = 2; transform[2] = 3;
	fillIndices(indices, 10, 10, 10, 10, 11, 10, 11, 11, 10, 10, 10, 11);

	// Make some fake data for the crapapapapa
	if (numPhases < 3) { cerr << "Blam\n"; system("pause"); exit(1); }
	// phi_v(x) = (v - p) dot (p - x)
	// p = (0.25, 0.5, 0.25). v = endpoint corresponding to material. x = location
	interfaceGrid[0][10][10][10] = (0.0 - 0.25) * (0.25 - 0.0) + (0.0 - 0.5) * (0.5 - 0.0) + (0.0 - 0.25) * (0.25 - 0.0);
	interfaceGrid[1][10][10][10] = (0.0 - 0.25) * (0.25 - 0.0) + (1.0 - 0.5) * (0.5 - 0.0) + (0.0 - 0.25) * (0.25 - 0.0);
	interfaceGrid[2][10][10][10] = (0.0 - 0.25) * (0.25 - 0.0) + (0.0 - 0.5) * (0.5 - 0.0) + (1.0 - 0.25) * (0.25 - 0.0);

	interfaceGrid[0][10][11][10] = (0.0 - 0.25) * (0.25 - 0.0) + (0.0 - 0.5) * (0.5 - 1.0) + (0.0 - 0.25) * (0.25 - 0.0);
	interfaceGrid[1][10][11][10] = (0.0 - 0.25) * (0.25 - 0.0) + (1.0 - 0.5) * (0.5 - 1.0) + (0.0 - 0.25) * (0.25 - 0.0);
	interfaceGrid[2][10][11][10] = (0.0 - 0.25) * (0.25 - 0.0) + (0.0 - 0.5) * (0.5 - 1.0) + (1.0 - 0.25) * (0.25 - 0.0);

	interfaceGrid[0][11][11][10] = (0.0 - 0.25) * (0.25 - 1.0) + (0.0 - 0.5) * (0.5 - 1.0) + (0.0 - 0.25) * (0.25 - 0.0);
	interfaceGrid[1][11][11][10] = (0.0 - 0.25) * (0.25 - 1.0) + (1.0 - 0.5) * (0.5 - 1.0) + (0.0 - 0.25) * (0.25 - 0.0);
	interfaceGrid[2][11][11][10] = (0.0 - 0.25) * (0.25 - 1.0) + (0.0 - 0.5) * (0.5 - 1.0) + (1.0 - 0.25) * (0.25 - 0.0);

	interfaceGrid[0][10][10][11] = (0.0 - 0.25) * (0.25 - 0.0) + (0.0 - 0.5) * (0.5 - 0.0) + (0.0 - 0.25) * (0.25 - 1.0);
	interfaceGrid[1][10][10][11] = (0.0 - 0.25) * (0.25 - 0.0) + (1.0 - 0.5) * (0.5 - 0.0) + (0.0 - 0.25) * (0.25 - 1.0);
	interfaceGrid[2][10][10][11] = (0.0 - 0.25) * (0.25 - 0.0) + (0.0 - 0.5) * (0.5 - 0.0) + (1.0 - 0.25) * (0.25 - 1.0);

	addPolygonsForTetrahedron(pp, ss, faceType, indices, transform);*/
	// END

	//setupGhostFractions();

	for (int ind = 0; ind < fineSize; ind++) {
		for (int mat = 0; mat < numPhases; mat++) { narrowBand(mat, ind) = false; }
		velocities(ind) = 0.0;
	}

	// Let's seed the interface with the given levelSet.
	for (Material = 0; Material < numPhases; Material++) {
		int nbSizeMat = 0;
		for (Index3 inds(m, n, p); inds.valid(); ++inds) {
			int i = inds.i, j = inds.j, k = inds.k;
			double v100, v010, v110, v000, v101, v011, v111, v001;

			int I = i / ns + 1, J = j / ns + 1, K = k / ns + 1;
			if (i == m) { I = coarseM; } if (j == n) { J = coarseN; } if (k == p) { K = coarseP; }
			int IJK = indexFromCoarse(I, J, K);

			int imod = i % ns, jmod = j % ns, kmod = k % ns, xDir, yDir, zDir;
			if (imod >= ns / 2) { xDir = 1; }
			else { xDir = -1; }
			if (jmod >= ns / 2) { yDir = 1; }
			else { yDir = -1; }
			if (kmod >= ns / 2) { zDir = 1; }
			else { zDir = -1; }

			if (I + xDir > coarseM || I + xDir < 1) { xDir = 0; }
			if (J + yDir > coarseN || J + yDir < 1) { yDir = 0; }
			if (K + zDir > coarseP || K + zDir < 1) { zDir = 0; }

			double xPos = (double)(imod) / ns - 0.5, yPos = (double)(jmod) / ns - 0.5,
				zPos = (double)(kmod) / ns - 0.5;
			xPos *= xDir; yPos *= yDir; zPos *= zDir;

			v000 = grid(Material, IJK);
			v100 = grid(Material, IJK + xDir * grid.right);
			v010 = grid(Material, IJK + yDir * grid.top);
			v110 = grid(Material, IJK + xDir * grid.right + yDir * grid.top);
			v001 = grid(Material, IJK + zDir * grid.forward);
			v101 = grid(Material, IJK + xDir * grid.right + zDir * grid.forward);
			v011 = grid(Material, IJK + yDir * grid.top + zDir * grid.forward);
			v111 = grid(Material, IJK + xDir * grid.right + yDir * grid.top + zDir * grid.forward);

			double volfrac = v000 + (v100 - v000) * xPos + // 1, x
				(v010 - v000) * yPos + (v001 - v000) * zPos + // y, z
				(v110 - v010 - v100 + v000) * xPos * yPos + // xy
				(v101 - v001 - v100 + v000) * xPos * zPos + // xz
				(v011 - v001 - v010 + v000) * yPos * zPos + // yz
				(v111 - v110 - v101 - v011 + v100 + v010 + v001 - v000) * xPos * yPos * zPos; // xyz
// Commenting this line lets the initial condition remain whatever it was before this function
// was called. (Often it is the exact distance function)
			interfaceGrid(Material, inds) = levelSet - volfrac;

			if (volfrac > 1e8 || volfrac < -1e8) {
				cerr << "What the hell?\n";
			}

			//narrowBandIterator[Material][nbSizeMat++] = Node3D(0.0, 0.0, i, j, k);
			//narrowBand[Material][i][j][k] = true;
			// Init narrow band
			const int scansize = 3;
			if (i < scansize || j < scansize || k < scansize) { continue; }

			int nMinus = 0, nPlus = 0;
			for (Index3 iprimes(i - scansize, j - scansize, k - scansize, i, j, k, m, n, p); iprimes.valid(); ++iprimes) {
				if (interfaceGrid(Material, iprimes) <= 0.0) {
					nMinus++;
				}
				else {
					nPlus++;
				}
			}

			if (nMinus > 0 && nPlus > 0) {
				for (Index3 iprimes(i - scansize, j - scansize, k - scansize, i, j, k, m, n, p); iprimes.valid(); ++iprimes) {
					if (!narrowBand(Material, iprimes)) {
						narrowBandIterator[Material][nbSizeMat++] = Node3D(0.0, 0.0, iprimes);
						narrowBand(Material, iprimes) = true;
					}
				}
			}
		}

		narrowBandSize[Material] = nbSizeMat;
		computeVActs();
		updateCoarseNarrowBand();

		cerr << "Narrow band size (Mat " << Material << "): " << nbSizeMat << endl;
	}

	fixLevelsetOverlap();
	checkResult();
	//reinitialize();
}

// Updated
double Interface3D::iterateLevelSetMethod(double epsilon, double timeStep, bool precise) {

	// Now, must choose the correct time step.
	//double timeStep = calculateBestTimeStep(epsilon, timeStep, currentMag);
	//cout << "Time step chosen was " << timeStep << " (k = " << k << ")\n";

	if (precise) {
		// Save old volume fractions (just to see the total change)
		for (Material = 0; Material < numPhases; Material++) {
			for (int inds = 0; inds < coarseSize; inds++) {
				if (narrowBand(Material, inds)) {
					vActsBack(Material, inds) = vActs(Material, inds);
					curvatureSums(inds) = 0.0;
				}
			}
		}

		updateGhostLayer();

		for (Material = 0; Material < numPhases; Material++) {
			// Flow kappa
			computeCurvatureFlow(epsilon);
			updateInterface(timeStep);

			// Needs to be updated.
			for (int volIter = 0; volIter < 1; volIter++) {
				computeSingleExactVAct(Material);
				//double volError = computeVActs();

				//cerr << "Vol iter " << volIter << "; error " << volError << endl;
				//checkResult();
				//if (volError < 1e-4) { break; }

				// timeStep is not important in this place
				extensions(false, timeStep, true);

				computeVolumeFlow(0.0, timeStep);
				// This time step is suppose to be around 1. The fact that it isn't means the linearization of the
				// change in area isn't working well... Check again with new volGradients!
				updateInterface(1.0);
			}
		}

		updateGhostLayer();
		//fixLevelsetOverlap();
		//reinitialize();

		double change = 0.0;
		computeAllVActs();
		// Check total change in volume
		for (Material = 0; Material < numPhases; Material++) {
			for (int inds = 0; inds < coarseSize; inds++) {
				if (narrowBand(Material, inds)) {
					change += fabs(vActsBack(Material, inds) - vActs(Material, inds));
				}
			}
		}

		return change;
	}
	else {
		updateGhostLayer();

		for (Material = 0; Material < numPhases; Material++) {
			for (int inds = 0; inds < coarseSize; inds++) {
				curvatureSums(inds) = 0.0;
			}
			//computeVActs();
			extensions(false, timeStep, false); // turned reinitialize on
			computeVolumeFlow(epsilon, timeStep);
			updateInterface(timeStep);
			//outputVelocities("V.vti");
		}

		updateGhostLayer();
		//fixLevelsetOverlap();		

		iter++;

		return computeAllVActs();
	}
}

void Interface3D::fixLevelsetOverlap() {
	double * phiTmp = new double[numPhases];

	for (Material = 0; Material < numPhases; Material++) {
		int nbSize = narrowBandSize[Material];
		for (int nbInd = 0; nbInd < nbSize; nbInd++) {
			Node3D nbElt = narrowBandIterator[Material][nbInd];
			int inds = nbElt.Index();

			bool skip = false;
			for (int mat = 0; mat < numPhases; mat++) {
				phiTmp[mat] = 200.0;
				if (mat < Material && narrowBand(mat, inds)) {
					skip = true; break;

				}
			}

			if (skip) { continue; }

			for (int mat = 0; mat < numPhases; mat++) {
				for (int othermat = 0; othermat < numPhases; othermat++) {
					if (othermat != mat && narrowBand(othermat, inds) && interfaceGrid(othermat, inds) < phiTmp[mat]) {
						phiTmp[mat] = interfaceGrid(othermat, inds);
					}
				}
			}

			for (int mat = 0; mat < numPhases; mat++) {
				if (narrowBand(mat, inds) && fabs(phiTmp[mat]) < 100.0) {
					interfaceGrid(mat, inds) = 0.5 * (interfaceGrid(mat, inds) - phiTmp[mat]);
				}
			}
		}
	}
}

void Interface3D::fixLevelsetOverlapNoNB() {
	// Need to fix level set overlap, except on full grid, no narrow band
	double * phiTmp = new double[numPhases];
	int ind = 0;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				for (int mat = 0; mat < numPhases; mat++) {
					phiTmp[mat] = 1e8;
					for (int mat2 = 0; mat2 < numPhases; mat2++) {
						if (mat2 == mat) { continue; }
						if (interfaceGrid(mat2, ind) < phiTmp[mat]) { phiTmp[mat] = interfaceGrid(mat2, ind); }
					}
				}

				for (int mat = 0; mat < numPhases; mat++) {
					interfaceGrid(mat, ind) = 0.5 * (interfaceGrid(mat, ind) - phiTmp[mat]);
				}
				ind++;
			}
		}
	}
	delete[] phiTmp;
}

// Depends on internal var Material
double Interface3D::computeVActs() {
	// Precompute all the volumes in each cell
	// Contains ghost layers.
	double changeInFraction = 0.0;
	double error = 0.0;
	totalSurfaceArea = 0.0;
	for (Index3 inds(1, 1, 1, coarseM, coarseN, coarseP, coarseM + 1, coarseN + 1, coarseP + 1); inds.valid(); ++inds) {
		if (coarseNarrowBand(Material, inds)) {
			double change = calculateVolumeImmersed(inds.i, inds.j, inds.k);
			changeInFraction += change;
			error += fabs(vActs(Material, inds) - grid(Material, inds));
		}
	}
	// Ghost layer should be unchanged; remains zero.
	double vol = 0.0;

	updateCoarseNarrowBand();

	return error;//changeInFraction;
}

double Interface3D::computeAllVActs() {
	double totalErr = 0.0;
	for (Material = 0; Material < numPhases; Material++) {
		totalErr += computeVActs();
	}
	//adjustVolumes(vActs);
	return totalErr;
}

double Interface3D::sigmoid(double x, double L)
{
	return 1.0 / (1.0 + exp(-L * (x - 0.5)));
}

// Updates based phi indexed by member variable Material
void Interface3D::updateInterface(double timeStep) {
	int nbSize = narrowBandSize[Material];
	for (int nbInd = 0; nbInd < nbSize; nbInd++) {
		Node3D nbElt = narrowBandIterator[Material][nbInd];
		int index = nbElt.Index();

		interfaceGrid(Material, index) -= timeStep * workGrid(index);
	}
}

// Computes the delta to the next iteration and puts it into workGrid.
// Depends on [Material]
double Interface3D::computeVolumeFlow(double epsilon, double timeStep) {
	//double changeInFraction = computeVActs();       

	//double maxF = 0.0;
	//int maxFi, maxFj;
	double magnitude = 0.0;

	// For each node in the narrow band
	int nbSize = narrowBandSize[Material];
	for (int nbInd = 0; nbInd < nbSize; nbInd++) {
		Node3D nbElt = narrowBandIterator[Material][nbInd];
		int centerIndex = nbElt.Index(), i, j, k;
		fineFromIndex(centerIndex, i, j, k);
		int I = (i / ns) + 1, J = (j / ns) + 1, K = (k / ns) + 1, IJK = indexFromCoarse(I, J, K);

		if (i == 0 || i == m || j == 0 || j == n || k == 0 || k == p) { workGrid(centerIndex) = 0.0; continue; }

		// Possibly temporary: only update nodes fully within the narrow band
		bool skip = false;
		for (Index3 iprimes(i - 1, j - 1, k - 1, i, j, k, m, n, p); iprimes.valid(); ++iprimes) {
			if (!narrowBand(Material, iprimes)) {
				skip = true;
				break;
			}
		}

		if (skip) { workGrid(centerIndex) = 0.0; continue; }

		// Need neighbording points for derivatives
		double right, left, up, down, forward, back;
		double topright, topleft, botright, botleft, center;
		double forwardright, forwardleft, backright, backleft;
		double forwardtop, forwardbot, backtop, backbot;

		center = interfaceGrid(Material, centerIndex);
		right = interfaceGrid(Material, centerIndex + interfaceGrid.right);
		left = interfaceGrid(Material, centerIndex + interfaceGrid.left);
		up = interfaceGrid(Material, centerIndex + interfaceGrid.top);
		down = interfaceGrid(Material, centerIndex + interfaceGrid.bottom);
		forward = interfaceGrid(Material, centerIndex + interfaceGrid.forward);
		back = interfaceGrid(Material, centerIndex + interfaceGrid.backward);

		topright = interfaceGrid(Material, centerIndex + interfaceGrid.top + interfaceGrid.right);
		topleft = interfaceGrid(Material, centerIndex + interfaceGrid.top + interfaceGrid.left);
		botright = interfaceGrid(Material, centerIndex + interfaceGrid.bottom + interfaceGrid.right);
		botleft = interfaceGrid(Material, centerIndex + interfaceGrid.bottom + interfaceGrid.left);

		forwardright = interfaceGrid(Material, centerIndex + interfaceGrid.forward + interfaceGrid.right);
		forwardleft = interfaceGrid(Material, centerIndex + interfaceGrid.forward + interfaceGrid.left);
		backright = interfaceGrid(Material, centerIndex + interfaceGrid.backward + interfaceGrid.right);
		backleft = interfaceGrid(Material, centerIndex + interfaceGrid.backward + interfaceGrid.left);

		forwardtop = interfaceGrid(Material, centerIndex + interfaceGrid.top + interfaceGrid.forward);
		forwardbot = interfaceGrid(Material, centerIndex + interfaceGrid.bottom + interfaceGrid.forward);
		backtop = interfaceGrid(Material, centerIndex + interfaceGrid.backward + interfaceGrid.top);
		backbot = interfaceGrid(Material, centerIndex + interfaceGrid.backward + interfaceGrid.bottom);

		// Centered derivatives
		double uxx = (right - 2 * center + left) / (hx*hx),
			uyy = (up - 2 * center + down) / (hy*hy),
			uzz = (forward - 2 * center + back) / (hz*hz);
		double uxy = ((topright - topleft) / (2 * hx) - (botright - botleft) / (2 * hx)) / (2 * hy),
			uxz = ((forwardright - forwardleft) / (2 * hx) - (backright - backleft) / (2 * hx)) / (2 * hz),
			uyz = ((forwardtop - forwardbot) / (2 * hy) - (backtop - backbot) / (2 * hy)) / (2 * hz);
		double ux = (right - left) / (2 * hx),
			uy = (up - down) / (2 * hy),
			uz = (forward - back) / (2 * hz);

		double gradLen = ux * ux + uy * uy + uz * uz, curvature;
		curvature = ((uyy + uzz) * ux * ux + (uxx + uzz) * uy * uy + (uxx + uyy) * uz * uz
			- 2 * ux * uy * uxy - 2 * ux * uz * uxz - 2 * uy * uz * uyz)
			/ (gradLen * sqrt(gradLen) + hx * hz);

		//if (fabs(curvature) > 1.0 / hx) {
		//   cerr << "Too much curvature\n";
		//}

		// Got to here.

		double dxp = (right - center) / hx,
			dxm = (center - left) / hx,
			dyp = (up - center) / hy,
			dym = (center - down) / hy,
			dzp = (forward - center) / hz,
			dzm = (center - back) / hz;

		double gradPlus = 0.0, gradMinus = 0.0;

		double gPlusScratch = max(dxm*fabs(dxm), 0.0) - min(dxp*fabs(dxp), 0.0)
			+ max(dym*fabs(dym), 0.0) - min(dyp*fabs(dyp), 0.0)
			+ max(dzm*fabs(dzm), 0.0) - min(dzp*fabs(dzp), 0.0);

		if (gPlusScratch >= 0) { gradPlus = sqrt(gPlusScratch); }
		//else { cerr << gPlusScratch << endl; system("PAUSE");}
		double gMinusScratch = -min(dxm*fabs(dxm), 0.0) + max(dxp*fabs(dxp), 0.0)
			- min(dym*fabs(dym), 0.0) + max(dyp*fabs(dyp), 0.0)
			- min(dzm*fabs(dzm), 0.0) + max(dzp*fabs(dzp), 0.0);
		if (gMinusScratch >= 0) { gradMinus = sqrt(gMinusScratch); }
		//else { cerr << gMinusScratch << endl; system("PAUSE");}

		double F = //smoothedSpeedFunc(i, j, k, timeStep);
			velocities(centerIndex);

		double curvatureCorrection = 1.0;
		//if (epsilon * curvatureSums[I][J][K] >  0.1 * hx / timeStep) {
		//    curvatureCorrection =  0.1 * hx / (timeStep * epsilon * curvatureSums[I][J][K]);
		//}
		double avgKappa = curvatureSums(IJK) / surfaceAreas(Material, IJK);
		//if (I == 7 && J == 7 && K == 7) {
	 //	   cerr << "/\n";
		//}
		if (avgKappa * hx * ns > 1.0) {
			curvatureCorrection = 0.25;
		}

		double workVal = (gradPlus * max(F, 0.0) + gradMinus * min(F, 0.0))
			- epsilon * curvature * sqrt(gradLen) * curvatureCorrection;
		workGrid(centerIndex) = workVal;

		if (fabs(workVal) > 1e4 || isNaN(workVal)) {
			cerr << "oh poo\n";
		}

		magnitude += fabs(workVal);

	} // For narrow band iterator

	//magnitude = changeInFraction;
	return magnitude; // pop pop!
}

// Also puts it into workGrid.
// [Material] dependent
double Interface3D::computeCurvatureFlow(double epsilon) {
	double magnitude = 0.0;

	// For each node in the narrow band
	int nbSize = narrowBandSize[Material];
	for (int nbInd = 0; nbInd < nbSize; nbInd++) {
		Node3D nbElt = narrowBandIterator[Material][nbInd];
		int centerIndex = nbElt.Index(), i, j, k;
		fineFromIndex(centerIndex, i, j, k);

		double gradLen, curvature;
		curvatureAt(centerIndex, i, j, k, curvature, gradLen);

		//int IJK = coarseIndexFromFine(i,  j, k);
		//double avgKappa = curvatureSums(IJK) / surfaceAreas(Material,IJK);
	   //if (I == 7 && J == 7 && K == 7) {
	//	   cerr << "/\n";
	   //}
	   //if (avgKappa * hx * ns > 1.0) {
	   //   curvature *= 0.25;
	   //}

		workGrid(centerIndex) = -epsilon * curvature * sqrt(gradLen);
	}

	return magnitude; // pop pop!
}

void Interface3D::curvatureAt(int centerIndex, int i, int j, int k, double & curvature, double & gradLen) {
	double right, left, up, down, forward, back;
	double topright, topleft, botright, botleft, center;
	double forwardright, forwardleft, backright, backleft;
	double forwardtop, forwardbot, backtop, backbot;

	if (i == 0 || j == 0 || k == 0 || i == m || j == n || k == p) { curvature = 0.0; gradLen = 1.0; return; }

	center = interfaceGrid(Material, centerIndex);
	right = interfaceGrid(Material, centerIndex + interfaceGrid.right);
	left = interfaceGrid(Material, centerIndex + interfaceGrid.left);
	up = interfaceGrid(Material, centerIndex + interfaceGrid.top);
	down = interfaceGrid(Material, centerIndex + interfaceGrid.bottom);
	forward = interfaceGrid(Material, centerIndex + interfaceGrid.forward);
	back = interfaceGrid(Material, centerIndex + interfaceGrid.backward);

	topright = interfaceGrid(Material, centerIndex + interfaceGrid.top + interfaceGrid.right);
	topleft = interfaceGrid(Material, centerIndex + interfaceGrid.top + interfaceGrid.left);
	botright = interfaceGrid(Material, centerIndex + interfaceGrid.bottom + interfaceGrid.right);
	botleft = interfaceGrid(Material, centerIndex + interfaceGrid.bottom + interfaceGrid.left);

	forwardright = interfaceGrid(Material, centerIndex + interfaceGrid.forward + interfaceGrid.right);
	forwardleft = interfaceGrid(Material, centerIndex + interfaceGrid.forward + interfaceGrid.left);
	backright = interfaceGrid(Material, centerIndex + interfaceGrid.backward + interfaceGrid.right);
	backleft = interfaceGrid(Material, centerIndex + interfaceGrid.backward + interfaceGrid.left);

	forwardtop = interfaceGrid(Material, centerIndex + interfaceGrid.top + interfaceGrid.forward);
	forwardbot = interfaceGrid(Material, centerIndex + interfaceGrid.bottom + interfaceGrid.forward);
	backtop = interfaceGrid(Material, centerIndex + interfaceGrid.backward + interfaceGrid.top);
	backbot = interfaceGrid(Material, centerIndex + interfaceGrid.backward + interfaceGrid.bottom);

	double uxx = (right - 2 * center + left) / (hx*hx),
		uyy = (up - 2 * center + down) / (hy*hy),
		uzz = (forward - 2 * center + back) / (hz*hz);
	double uxy = ((topright - topleft) / (2 * hx) - (botright - botleft) / (2 * hx)) / (2 * hy),
		uxz = ((forwardright - forwardleft) / (2 * hx) - (backright - backleft) / (2 * hx)) / (2 * hz),
		uyz = ((forwardtop - forwardbot) / (2 * hy) - (backtop - backbot) / (2 * hy)) / (2 * hz);
	double ux = (right - left) / (2 * hx),
		uy = (up - down) / (2 * hy),
		uz = (forward - back) / (2 * hz);

	gradLen = ux * ux + uy * uy + uz * uz, curvature;
	curvature = ((uyy + uzz) * ux * ux + (uxx + uzz) * uy * uy + (uxx + uyy) * uz * uz
		- 2 * ux * uy * uxy - 2 * ux * uz * uxz - 2 * uy * uz * uyz)
		/ (gradLen * sqrt(gradLen) + hx * hz);
}

// Reinitializes the level set function to a signed distance function.
// TEMP: Transitioning to chopp initialization of FMM, so has both currenty
void Interface3D::extensions(bool reinitialize, double timeStep, bool precise) {

	EikonalHeap3D setup((m + 1)*(n + 1), m + 1, n + 1, p + 1);
	EikonalHeap3D burn(50, m + 1, n + 1, p + 1);

	// Find interface and add to queue.
	int nbSize = narrowBandSize[Material];
	for (int nbInd = 0; nbInd < nbSize; nbInd++) {
		Node3D nbElt = narrowBandIterator[Material][nbInd];
		int index = nbElt.Index(), i, j, k;
		fineFromIndex(index, i, j, k);

		// Since little cube goes from i to i+1, we don't want to go out of bounds
		if (i == m || j == n || k == p) {
			continue;
		}

		int nMinus = 0, nPlus = 0;
		for (Index3 iprimes(i, j, k, i + 1, j + 1, k + 1, m, n, p); iprimes.valid(); ++iprimes) {
			if (interfaceGrid(Material, iprimes) > 0.0) {
				nPlus++;
			}
			else {
				nMinus--;
			}
		}
		if (nMinus == 0 || nPlus == 0) { continue; }

		double kappaSums = 0.0, maxKappa = 0.0;

		// Each cube is divided into six tetrahedra
		for (int T = 0; T < 6; T++) {
			double newKappa = getMaxCurvatureOnInterface(index + tetrahedraOffsets[T][0], index + tetrahedraOffsets[T][1],
				index + tetrahedraOffsets[T][2], index + tetrahedraOffsets[T][3]);
			kappaSums += newKappa;
		}

		int I = (i / ns) + 1, J = (j / ns) + 1, K = (k / ns) + 1;
		curvatureSums(I, J, K) += kappaSums;

		// Chopp
		double coeff[64], stencil[64];
		// Init stencil
		int c = 0;
		for (int kp = -1; kp <= 2; kp++) {
			for (int jp = -1; jp <= 2; jp++) {
				for (int ip = -1; ip <= 2; ip++) {
					int ic = i + ip, jc = j + jp, kc = k + kp;
					if (ic <= 0) { ic = 0; } if (ic > m) { ic = m; }
					if (jc <= 0) { jc = 0; } if (jc > n) { jc = n; }
					if (kc <= 0) { kc = 0; } if (kc > p) { kc = p; }

					stencil[c++] = interfaceGrid(Material, ic, jc, kc);
				}
			}
		}

		tricubicCoeff(coeff, stencil);


		for (Index3 iprimes(i, j, k, i + 1, j + 1, k + 1, m, n, p); iprimes.valid(); ++iprimes) {
			double x = (double)iprimes.i - i, y = (double)iprimes.j - j, z = (double)iprimes.k - k;
			tricubicDistanceSetup(x, y, z, iprimes.i, iprimes.j, iprimes.k, coeff, timeStep, setup);
		}

		burn.discardAll();
	}

	// First, let's get the narrow banding to work on !precise mode.
	if (precise) {
		// setup now contains a heap of interface points; grab all items and throw onto our gradientElement list
		EikonalHeap3D workHeap(m*n, m + 1, n + 1, p + 1);

		// Clear previous volume gradient info
		int nbSize = narrowBandSize[Material];
		for (int nbInd = 0; nbInd < nbSize; nbInd++) {
			Node3D nbElt = narrowBandIterator[Material][nbInd];
			volumeGradients(nbElt.Index()).setU(-1.0);
		}
		allGradients.clear();

		for (int ind = 0; ind < coarseSize; ind++) {
			gradientBand(ind).clear();
		}

		numGradElts = 0;
		while (!setup.isEmpty()) {
			Node3D next = setup.remove();
			//int i = next.X(), j = next.Y(), k = next.Z();
			//int I = 1 + (i / ns), J = 1 + (j / ns), K = 1 + (k / ns);
			//int iStart = (I-1)*ns, jStart = (J-1)*ns, kStart = (K-1)*ns;
			gradientElement & ge = volumeGradients(next.Index());
			ge.setU(next.getU());
			ge.clearContributions();

			//	 //volumeGradients[I][J][K].push_back(ge);

				 //volumeGradients[numGradElts++] = ge;
				 //if (numGradElts == capacityGradElts) {
					// cerr << "Reallocating interface size\n";
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
			int index = next.Index(), i, j, k;
			fineFromIndex(index, i, j, k);
			int I = (i / ns) + 1, J = (j / ns) + 1, K = (k / ns) + 1;

			// Boundary nodes should not need to be updated, as they will be extrapolated later anyway.
			if (i == m || j == n || k == p || i == 0 || j == 0 || k == 0) { continue; }
			// EXPERIMENTING
			gradientElement & ge = volumeGradients(index);

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
				gradientBand(I + Ip, J + Jp, K + Kp).push_back(&ge);
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
			int index = next.Index(), i, j, k;
			fineFromIndex(index, i, j, k);
			int I = (i / ns) + 1, J = (j / ns) + 1, K = (k / ns) + 1;

			double speed = 0.0;
			for (int offset = 0; offset < 27; offset++) {
				double contribution = volumeGradients(index).getContribution(offset);
				if (contribution < 0.0 || contribution > 0.0) {
					int Ip, Jp, Kp;
					IndexDeltasFromOffset(Ip, Jp, Kp, offset);
					speed += lambdas(I + Ip, J + Jp, K + Kp) * contribution;
					if (isNaN(contribution) || isNaN(lambdas(I + Ip, J + Jp, K + Kp))) {
						cerr << "Nice day isn't it\n";
					}
				}
			}

			next.setV(speed);
			setup.insert(next);
		}
	} // precise

	// For each node in the narrow band
	for (int nbInd = 0; nbInd < nbSize; nbInd++) {
		Node3D nbElt = narrowBandIterator[Material][nbInd];
		int ind = nbElt.Index();

		// Old comment. If there's a bug maybe it will offer insight.
		// Set some surrounding entries to infinity as well (Don't want FMM to check non narrow band elements)

		//signs[i][j][k] = interfaceGrid[Material][i][j][k] > 0.0 ? 1.0 : -1.0;
		workGrid(ind) = INFTY;
		if (reinitialize) {
			narrowBand(Material, ind) = false;
		}
	}

	// Issue: Signs, if reinitializing, may not have been initialized

	// Issue: Assumes square grid for now
	// WARNING: h = hx only works for isotropic grid right now
	double radius = hx * ns * 10.0;
	doFastMarch3D(workGrid, velocities, setup, reinitialize, radius);

	// Make signed distance function
	// For each node in the narrow band
	if (reinitialize) {
		nbSize = narrowBandSize[Material];
		for (int nbInd = 0; nbInd < nbSize; nbInd++) {
			Node3D nbElt = narrowBandIterator[Material][nbInd];
			int ind = nbElt.Index();

			double fmmResult = workGrid(ind);
			if (fmmResult < INFTY) { // in theory should not need this
				interfaceGrid(Material, ind) = fmmResult * signs(ind); // I think this is causing probs
			}
		}
	}
}

// Depends on value of Material
double Interface3D::getMaxCurvatureOnInterface(int indA, int indB, int indC, int indD) {

	double a = interfaceGrid(Material, indA),
		b = interfaceGrid(Material, indB),
		c = interfaceGrid(Material, indC),
		d = interfaceGrid(Material, indD);

	double nMinus = 0;
	if (a <= 0) { nMinus++; }
	if (b <= 0) { nMinus++; }
	if (c <= 0) { nMinus++; }
	if (d <= 0) { nMinus++; }

	if (nMinus == 0 || nMinus == 4) { return 0.0; }

	int i, j, k;

	double curvA, curvB, curvC, curvD, scratch, curvMax = 0.0;

	fineFromIndex(indA, i, j, k);
	curvatureAt(indA, i, j, k, curvA, scratch);
	fineFromIndex(indB, i, j, k);
	curvatureAt(indB, i, j, k, curvB, scratch);
	fineFromIndex(indC, i, j, k);
	curvatureAt(indC, i, j, k, curvC, scratch);
	fineFromIndex(indD, i, j, k);
	curvatureAt(indD, i, j, k, curvD, scratch);

	if (a * b <= 0) {
		curvMax = max(curvMax, fabs(-a / (b - a) * curvB + b / (b - a) * curvA));
	}

	if (a * c <= 0) {
		curvMax = max(curvMax, fabs(-a / (c - a) * curvC + c / (c - a) * curvA));
	}

	if (a * d <= 0) {
		curvMax = max(curvMax, fabs(-a / (d - a) * curvD + d / (d - a) * curvA));
	}

	if (b * c <= 0) {
		curvMax = max(curvMax, fabs(-b / (c - b) * curvC + c / (c - b) * curvB));
	}

	if (b * d <= 0) {
		curvMax = max(curvMax, fabs(-b / (d - b) * curvD + d / (d - b) * curvB));
	}

	if (c * d <= 0) {
		curvMax = max(curvMax, fabs(-c / (d - c) * curvD + d / (d - c) * curvC));
	}

	return curvMax;
}

double Interface3D::smoothedSpeedFunc(double i, double j, double k, double timeStep) {

	//return tricubicSpeedFunc(i,j,k,timeStep); // what if we replaced it?

	if (i > m) { i = m - 1e-10; } if (i < 0) { i = 0.0; }
	if (j > n) { j = n - 1e-10; } if (j < 0) { j = 0.0; }
	if (k > p) { k = p - 1e-10; } if (k < 0) { k = 0.0; }

	int I = (int)i / ns + 1, J = (int)j / ns + 1, K = (int)k / ns + 1;
	int iStart = (I - 1)*ns, jStart = (J - 1)*ns, kStart = (K - 1)*ns;
	int IJK = indexFromCoarse(I, J, K);

	double F000, F001, F010, F011, F100, F101, F110, F111;
	int iDir, jDir, kDir;
	if (i - iStart <= ns / 2) { iDir = -1; }
	else { iDir = 1; }
	if (j - jStart <= ns / 2) { jDir = -1; }
	else { jDir = 1; }
	if (k - kStart <= ns / 2) { kDir = -1; }
	else { kDir = 1; }

	F000 = isolatedSpeedFunc(grid(Material, IJK), vActs(Material, IJK), surfaceAreas(Material, IJK), timeStep);

	double xPos = 0.0, yPos = 0.0, zPos = 0.0;

	/*if (i - iStart < 1.0) {
		xPos = (1.0 - i + iStart) / 2.0;
	} else if (i - iStart > ns - 1.0) {
		xPos = (i - iStart - ns + 1.0) / 2.0;
	}

	if (j - jStart < 1.0) {
		yPos = (1.0 - j + jStart) / 2.0;
	} else if (j - jStart > ns - 1.0) {
		yPos = (j - jStart - ns + 1.0) / 2.0;
	}

	if (k - kStart < 1.0) {
		zPos = (1.0 - k + kStart) / 2.0;
	} else if (k - kStart > ns - 1.0) {
		zPos = (k - kStart - ns + 1.0) / 2.0;
	}*/

	xPos = (i - iStart - ns / 2.0) / ((double)ns) * iDir;
	yPos = (j - jStart - ns / 2.0) / ((double)ns) * jDir;
	zPos = (k - kStart - ns / 2.0) / ((double)ns) * kDir;

	// Replaces ghost cells
	if (I + iDir == 0 || I + iDir == coarseM + 1) { iDir = 0; }
	if (J + jDir == 0 || J + jDir == coarseN + 1) { jDir = 0; }
	if (K + kDir == 0 || K + kDir == coarseP + 1) { kDir = 0; }

	int rightInd = IJK + grid.right * iDir, topRightInd = rightInd + grid.top * jDir, topInd = IJK + grid.top * jDir,
		forwardInd = IJK + grid.forward * kDir, forwardRightInd = rightInd + grid.forward * kDir,
		forwardTopInd = topInd + grid.forward * kDir, forwardTopRightInd = forwardTopInd + grid.right * iDir;

	F100 = isolatedSpeedFunc(grid(Material, rightInd), vActs(Material, rightInd), surfaceAreas(Material, rightInd), timeStep);
	F010 = isolatedSpeedFunc(grid(Material, topInd), vActs(Material, topInd), surfaceAreas(Material, topInd), timeStep);
	F110 = isolatedSpeedFunc(grid(Material, topRightInd), vActs(Material, topRightInd), surfaceAreas(Material, topRightInd), timeStep);
	F001 = isolatedSpeedFunc(grid(Material, forwardInd), vActs(Material, forwardInd), surfaceAreas(Material, forwardInd), timeStep);
	F101 = isolatedSpeedFunc(grid(Material, forwardRightInd), vActs(Material, forwardRightInd), surfaceAreas(Material, forwardRightInd), timeStep);
	F011 = isolatedSpeedFunc(grid(Material, forwardTopInd), vActs(Material, forwardTopInd), surfaceAreas(Material, forwardTopInd), timeStep);
	F111 = isolatedSpeedFunc(grid(Material, forwardTopRightInd), vActs(Material, forwardTopRightInd), surfaceAreas(Material, forwardTopRightInd), timeStep);

	//double xPos = sigmoid( ((i - iStart)/(double)ns - 0.5) * iDir, 2.5 ),
	//       yPos = sigmoid( ((j - jStart)/(double)ns - 0.5) * jDir, 2.5 ),
	//       zPos = sigmoid( ((k - kStart)/(double)ns - 0.5) * kDir, 2.5 );
	//double xPos = ((i - iStart)/(double)ns - 0.5) * iDir,
	//       yPos = ((j - jStart)/(double)ns - 0.5) * jDir,
	//       zPos = ((k - kStart)/(double)ns - 0.5) * kDir;

	double F = F000 + (F100 - F000) * xPos + // 1, x
		(F010 - F000) * yPos + (F001 - F000) * zPos + // y, z
		(F110 - F010 - F100 + F000) * xPos * yPos + // xy
		(F101 - F001 - F100 + F000) * xPos * zPos + // xz
		(F011 - F001 - F010 + F000) * yPos * zPos + // yz
		(F111 - F110 - F101 - F011 + F100 + F010 + F001 - F000) * xPos * yPos * zPos; // xyz

	return F;
}

double Interface3D::tricubicSpeedFunc(double i, double j, double k, double timeStep) {
	int I = (int)(i / ns + 0.5), J = (int)(j / ns + 0.5), K = (int)(k / ns + 0.5);
	int iStart = (I - 1)*ns, jStart = (J - 1)*ns, kStart = (K - 1)*ns;

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
		int II = xcoords[xc], JJ = ycoords[yc], KK = zcoords[zc], Index = indexFromCoarse(II, JJ, KK);
		rhs[c] = isolatedSpeedFunc(grid(Material, Index), vActs(Material, Index), surfaceAreas(Material, Index), timeStep);
	}

	return tricubicInterpolation(xPos, yPos, zPos, rhs);
}

void Interface3D::tricubicCoeff(double * coeff, double * __restrict rhs) {
	int c = 0;
	for (int i = 0; i < 64; i++) {
		double coeffTmp = 0.0;
		for (; c < tricubicRowOffsets[i]; c++) {
			coeffTmp += sparseCoeffs[c] * rhs[tricubicSparseIndices[c]];
		}
		coeff[i] = coeffTmp;
	}
}

double Interface3D::tricubicInterpolation(double x, double y, double z, double * __restrict rhs) {
	double coeff[64];
	tricubicCoeff(coeff, rhs);
	return tricubicEval(x, y, z, coeff);
}

double Interface3D::tricubicEval(double x, double y, double z, double * __restrict coeff) {
	double sum = 0.0;
	for (int i = 3; i >= 0; i--) {
		//s *= x;
		double ysum = 0.0;
		for (int j = 3; j >= 0; j--) {
			//s *= ;
			double zsum = 0.0;
			for (int k = 3; k >= 0; k--) {
				zsum = z * zsum + coeff[i + 4 * j + 16 * k];
			}
			ysum = ysum * y + zsum;
		}
		sum = sum * x + ysum;
	}

	return sum;
}

void Interface3D::tricubicDerivatives(double x, double y, double z, double & dx, double & dy, double & dz, double * __restrict coeff) {
	dx = 0;
	for (int i = 3; i >= 1; i--) {
		//s *= x;
		double ysum = 0.0;
		for (int j = 3; j >= 0; j--) {
			//s *= ;
			double zsum = 0.0;
			for (int k = 3; k >= 0; k--) {
				zsum = z * zsum + coeff[i + 4 * j + 16 * k];
			}
			ysum = ysum * y + zsum;
		}
		dx = dx * x*i + ysum;
	}

	dy = 0;
	for (int i = 3; i >= 0; i--) {
		//s *= x;
		double ysum = 0.0;
		for (int j = 3; j >= 1; j--) {
			//s *= ;
			double zsum = 0.0;
			for (int k = 3; k >= 0; k--) {
				zsum = z * zsum + coeff[i + 4 * j + 16 * k];
			}
			ysum = ysum * y*j + zsum;
		}
		dy = dy * x + ysum;
	}

	dz = 0;
	for (int i = 3; i >= 0; i--) {
		//s *= x;
		double ysum = 0.0;
		for (int j = 3; j >= 0; j--) {
			//s *= ;
			double zsum = 0.0;
			for (int k = 3; k >= 1; k--) {
				zsum = z * zsum*k + coeff[i + 4 * j + 16 * k];
			}
			ysum = ysum * y*j + zsum;
		}
		dz = dz * x + ysum;
	}
}

void Interface3D::tricubicDistanceSetup(double x0, double y0, double z0, int i, int j, int k, double coeff[64], double timeStep, EikonalHeap3D & setup) {
	const int maxNewtonIter = 20; // from Chopp
	double x = x0, y = y0, z = z0;
	double change;
	int newtonIter;
	for (newtonIter = 0; newtonIter < maxNewtonIter; newtonIter++) {
		// Compute gradient of p(x,y,z)
		double dx, dy, dz;
		tricubicDerivatives(x, y, z, dx, dy, dz, coeff);

		double eval = tricubicEval(x, y, z, coeff),
			gradMag = (dx*dx + dy * dy + dz * dz);
		if (fabs(gradMag) < 1e-10) { gradMag = 1.0; }

		double delta1Coeff = -eval / gradMag; // delta1 = delta1Coeff * grad p

		double vx = x0 - x, vy = y0 - y, vz = z0 - z;
		double delta2Coeff = (vx * dx + vy * dy + vz * dz) / gradMag; // also the coeff of gradP

		double delta2x = vx - dx * delta2Coeff, delta2y = vy - dy * delta2Coeff, delta2z = vz - dz * delta2Coeff;

		x += dx * (delta1Coeff - delta2Coeff) + vx;
		y += dy * (delta1Coeff - delta2Coeff) + vy;
		z += dz * (delta1Coeff - delta2Coeff) + vz;

		change = sqrt(delta2x*delta2x + delta2y * delta2y + delta2z * delta2z + eval * eval / gradMag);
		if (isNaN(change)) {
			cerr << "Breakdance!\n";
			change = 1.0;
			break;
		}
		if (change < 0.001) {
			break;
		}
	}

	if (change > 0.001) {// || x < -0.5 || x > 1.0 || y < 0.0 || y > 1.0 || z < 0.0 || z > 1.0) {
		// Divergence
		//cerr << "Diverged, son!\n";
	}
	else {
		double vx = x - x0, vy = y - y0, vz = z - z0,
			dist = sqrt(vx*vx + vy * vy + vz * vz),
			velocity = //tricubicSpeedFunc(i + vx, j + vy, k + vz, timeStep);
			smoothedSpeedFunc(i + vx, j + vy, k + vz, timeStep);
		setup.insert(Node3D(hx * dist, velocity, indexFromFine(i, j, k)));
		//cerr << "Tricubic found for point " << i << " " << j << " " << k << " " << dist << endl;
	}
}

// Depends on value of Material
void Interface3D::updateGradients(gradientElement & ge, EikonalHeap3D & workHeap) {
	// Algorithm: Compute volume flow at point (i,j,k) (store in workGrid) with some small tau
	// Compute tetrahedon changes in surrounding tetra
	// or maybe just surrounding cubes if we are lazy
	// then do get.setContribution(offsets, dvol);
	int i = ge.X(), j = ge.Y(), k = ge.Z();
	int I = i / ns + 1, J = j / ns + 1, K = k / ns + 1;
	int centerIndex = indexFromFine(i, j, k);

	double backupVal = interfaceGrid(Material, centerIndex);
	const double tau = 1e-3;

	double right, left, up, down, forward, back, center;

	center = interfaceGrid(Material, centerIndex);
	right = interfaceGrid(Material, centerIndex + interfaceGrid.right);
	left = interfaceGrid(Material, centerIndex + interfaceGrid.left);
	up = interfaceGrid(Material, centerIndex + interfaceGrid.top);
	down = interfaceGrid(Material, centerIndex + interfaceGrid.bottom);
	forward = interfaceGrid(Material, centerIndex + interfaceGrid.forward);
	back = interfaceGrid(Material, centerIndex + interfaceGrid.backward);

	double ux = (right - left) / (2 * hx),
		uy = (up - down) / (2 * hy),
		uz = (forward - back) / (2 * hz);

	double gradLen = ux * ux + uy * uy + uz * uz;

	double dxp = (right - center) / hx,
		dxm = (center - left) / hx,
		dyp = (up - center) / hy,
		dym = (center - down) / hy,
		dzp = (forward - center) / hz,
		dzm = (center - back) / hz;
	double gradPlus = 0.0;
	double gPlusScratch = max(dxm*fabs(dxm), 0.0) - min(dxp*fabs(dxp), 0.0)
		+ max(dym*fabs(dym), 0.0) - min(dyp*fabs(dyp), 0.0)
		+ max(dzm*fabs(dzm), 0.0) - min(dzp*fabs(dzp), 0.0);
	if (gPlusScratch >= 0) { gradPlus = sqrt(gPlusScratch); }
	//else { cerr << gPlusScratch << endl; system("PAUSE");}
	double newVal = interfaceGrid(Material, centerIndex) - tau * gradPlus;

	// Find change in volume in those boxes sharing a vertex in common with (i,j,k)
	for (Index3 iprimes(i - 1, j - 1, k - 1, i, j, k, m, n, p); iprimes.valid(); ++iprimes) {
		int ii = iprimes.i, jj = iprimes.j, kk = iprimes.k;
		int Ip = 1 + (ii) / ns - I, Jp = 1 + (jj) / ns - J, Kp = 1 + (kk) / ns - K;
		int offset = offsetFromDeltas(Ip, Jp, Kp);
		double newVolume, oldVolume, dummyArea;

		// Change to outward flow at that vertex
		interfaceGrid(Material, centerIndex) = newVal;

		newVolume = fluidAndAreaBox(interfaceGrid(Material, iprimes),
			interfaceGrid(Material, iprimes + interfaceGrid.right),
			interfaceGrid(Material, iprimes + interfaceGrid.top),
			interfaceGrid(Material, iprimes + interfaceGrid.top + interfaceGrid.right),
			interfaceGrid(Material, iprimes + interfaceGrid.forward),
			interfaceGrid(Material, iprimes + interfaceGrid.right + interfaceGrid.forward),
			interfaceGrid(Material, iprimes + interfaceGrid.top + interfaceGrid.forward),
			interfaceGrid(Material, iprimes + interfaceGrid.top + interfaceGrid.right + interfaceGrid.forward), dummyArea);

		// Check original volume in cell
		interfaceGrid(Material, centerIndex) = backupVal;
		oldVolume = fluidAndAreaBox(interfaceGrid(Material, iprimes),
			interfaceGrid(Material, iprimes + interfaceGrid.right),
			interfaceGrid(Material, iprimes + interfaceGrid.top),
			interfaceGrid(Material, iprimes + interfaceGrid.top + interfaceGrid.right),
			interfaceGrid(Material, iprimes + interfaceGrid.forward),
			interfaceGrid(Material, iprimes + interfaceGrid.right + interfaceGrid.forward),
			interfaceGrid(Material, iprimes + interfaceGrid.top + interfaceGrid.forward),
			interfaceGrid(Material, iprimes + interfaceGrid.top + interfaceGrid.right + interfaceGrid.forward), dummyArea);

		ge.addContribution(offset, fabs(newVolume - oldVolume) / (tau * ns * ns * ns));
	}

	//if (i == 22 && j == 34 && k == 35) {
	//	cerr << "Check it, yo\n";
	//}

	ge.normalize();
}

void Interface3D::computeVolumeOfGradients(EikonalHeap3D & workHeap) {
	int num = 0;

	// For each gradient we must do a FMM
	// but makes sense not to initialize each time
	// so that is why we input workHeap

	// We actually require a for-each-cell calculation here
	for (Index3 coarseInds(1, 1, 1, coarseM, coarseN, coarseP, coarseM + 1, coarseN + 1, coarseP + 1); coarseInds.valid(); ++coarseInds) {
		if (surfaceAreas(Material, coarseInds) <= 1e-14) continue;
		// Update interface in this cell by adding a constant tau = 0.01 of our volume gradients
		// Store result in workGrid
		// (f(x + tau v) - f(x)) / tau
		double tau = 0.005; //min(max(fabs(grid[I][J][K] - vActs[I][J][K]), 0.0001), 0.2);

		int I = coarseInds.i, J = coarseInds.j, K = coarseInds.k;
		int cellStartI = (I - 1)*ns, cellStartJ = (J - 1)*ns, cellStartK = (K - 1)*ns;
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
				IStar = i / ns + 1, JStar = j / ns + 1, KStar = k / ns + 1;

			int Ip = I - IStar, Jp = J - JStar, Kp = K - KStar;
			double speed = ge->getContribution(offsetFromDeltas(Ip, Jp, Kp));
			//			..if (isNaN(speed)) {
			if (speed < 0) {
				cerr << "Come on, then!\n";
			}
			workHeap.insert(Node3D(ge->U(), speed, indexFromFine(i, j, k)));
		}

		doFastMarch3D(workGrid, velocities, workHeap, false, 2.0 * hx);
		workHeap.discardAll();

		int flowStartI = cellStartI - 2, flowStartJ = cellStartJ - 2, flowStartK = cellStartK - 2,
			flowEndI = cellStartI + ns + 2, flowEndJ = cellStartJ + ns + 2, flowEndK = cellStartK + ns + 2;

		flowStartI = max(1, flowStartI);
		flowStartJ = max(1, flowStartJ);
		flowStartK = max(1, flowStartK);
		flowEndI = min(m - 1, flowEndI);
		flowEndJ = min(n - 1, flowEndJ);
		flowEndK = min(p - 1, flowEndK);

		for (Index3 flowIter(flowStartI, flowStartJ, flowStartK, flowEndI, flowEndJ, flowEndK, m, n, p); flowIter.valid(); ++flowIter) {

			if (workGrid(flowIter) == INFTY) {
				velocities(flowIter) = 0.0;
			}

			double right, left, up, down, forward, back, center;
			center = interfaceGrid(Material, flowIter);
			right = interfaceGrid(Material, flowIter + interfaceGrid.right);
			left = interfaceGrid(Material, flowIter + interfaceGrid.left);
			up = interfaceGrid(Material, flowIter + interfaceGrid.top);
			down = interfaceGrid(Material, flowIter + interfaceGrid.bottom);
			forward = interfaceGrid(Material, flowIter + interfaceGrid.forward);
			back = interfaceGrid(Material, flowIter + interfaceGrid.backward);

			double ux = (right - left) / (2 * hx),
				uy = (up - down) / (2 * hy),
				uz = (forward - back) / (2 * hz);

			double gradLen = ux * ux + uy * uy + uz * uz;
			double dxp = (right - center) / hx,
				dxm = (center - left) / hx,
				dyp = (up - center) / hy,
				dym = (center - down) / hy,
				dzp = (forward - center) / hz,
				dzm = (center - back) / hz;

			double gradPlus = 0.0;

			double gPlusScratch = max(dxm*fabs(dxm), 0.0) - min(dxp*fabs(dxp), 0.0)
				+ max(dym*fabs(dym), 0.0) - min(dyp*fabs(dyp), 0.0)
				+ max(dzm*fabs(dzm), 0.0) - min(dzp*fabs(dzp), 0.0);
			if (gPlusScratch >= 0) { gradPlus = sqrt(gPlusScratch); }

			//workGrid[i][j][k] -= tau * velocities[i][j][k] * gradPlus;
			workGrid(flowIter) = interfaceGrid(Material, flowIter) - tau * velocities(flowIter) * gradPlus;
			if (velocities(flowIter) < 0.0) {
				cerr << "Chubby cheeks\n";
			}

			if (isNaN(workGrid(flowIter))) {
				cerr << "Ok dudes\n";
			}

			//interfaceGrid[i][j][k] -= tau * velocities[i][j][k] * gradPlus;
		}

		// Update ghost layer on cells, if applicable.
		// So far, corners not taken into account
		if (flowStartI == 1 || flowEndI == m - 1) {
			int ipos = flowStartI == 1 ? 0 : m, idir = flowStartI == 1 ? 1 : -1;
			// Extrapolate in X.
			for (Index3 extrap(ipos, flowStartJ, flowStartK, ipos, flowEndJ, flowEndK, m, n, p); extrap.valid(); ++extrap) {
				workGrid(extrap) = 2.0 * workGrid(extrap + idir * interfaceGrid.right) - workGrid(extrap + 2 * idir * interfaceGrid.right);
			}
		}
		if (flowStartJ == 1 || flowEndJ == n - 1) {
			int jpos = flowStartJ == 1 ? 0 : n, jdir = flowStartJ == 1 ? 1 : -1;
			// Extrapolate in Y.
			for (Index3 extrap(flowStartI, jpos, flowStartK, flowEndI, jpos, flowEndK, m, n, p); extrap.valid(); ++extrap) {
				workGrid(extrap) = 2.0 * workGrid(extrap + jdir * interfaceGrid.top) - workGrid(extrap + 2 * jdir * interfaceGrid.top);
			}
		}
		if (flowStartK == 1 || flowEndK == p - 1) {
			int kpos = flowStartK == 1 ? 0 : m, kdir = flowStartK == 1 ? 1 : -1;
			// Extrapolate in Zed.
			for (Index3 extrap(flowStartI, flowStartJ, kpos, flowEndI, flowEndJ, kpos, m, n, p); extrap.valid(); ++extrap) {
				workGrid(extrap) = 2.0 * workGrid(extrap + kdir * interfaceGrid.forward) - workGrid(extrap + 2 * kdir * interfaceGrid.forward);
			}
		}

		num++;

		// Clear entries
		gradientElement & whichElement = volumeOfGradients(coarseInds);
		whichElement.clearContributions();

		// Then compute volume difference between interfaceGrid and workGrid for each tetrahedron
		// startin from iStart-1 to iStart + ns (+1 grid cell around coarse cell)
		// Add changed volume to appropriate volumeOfGradient

		int volStartI = cellStartI - 1, volStartJ = cellStartJ - 1, volStartK = cellStartK - 1,
			volEndI = cellStartI + ns, volEndJ = cellStartJ + ns, volEndK = cellStartK + ns;

		volStartI = max(0, volStartI);
		volStartJ = max(0, volStartJ);
		volStartK = max(0, volStartK);
		volEndI = min(m - 1, volEndI);
		volEndJ = min(n - 1, volEndJ);
		volEndK = min(p - 1, volEndK);

		for (Index3 deltaVIndices(volStartI, volStartJ, volStartK, volEndI, volEndJ, volEndK, m, n, p); deltaVIndices.valid(); ++deltaVIndices) {
			int Ip = -1 + (deltaVIndices.i - cellStartI + ns) / ns,
				Jp = -1 + (deltaVIndices.j - cellStartJ + ns) / ns,
				Kp = -1 + (deltaVIndices.k - cellStartK + ns) / ns;
			int index = offsetFromDeltas(Ip, Jp, Kp);

			double dummyarea = 0.0;
			double oldVolume = fluidAndAreaBox(interfaceGrid(Material, deltaVIndices),
				interfaceGrid(Material, deltaVIndices + interfaceGrid.right),
				interfaceGrid(Material, deltaVIndices + interfaceGrid.top),
				interfaceGrid(Material, deltaVIndices + interfaceGrid.right + interfaceGrid.top),
				interfaceGrid(Material, deltaVIndices + interfaceGrid.forward),
				interfaceGrid(Material, deltaVIndices + interfaceGrid.forward + interfaceGrid.right),
				interfaceGrid(Material, deltaVIndices + interfaceGrid.forward + interfaceGrid.top),
				interfaceGrid(Material, deltaVIndices + interfaceGrid.forward + interfaceGrid.right + interfaceGrid.top), dummyarea);

			double newVolume = fluidAndAreaBox(workGrid(deltaVIndices),
				workGrid(deltaVIndices + interfaceGrid.right),
				workGrid(deltaVIndices + interfaceGrid.top),
				workGrid(deltaVIndices + interfaceGrid.right + interfaceGrid.top),
				workGrid(deltaVIndices + interfaceGrid.forward),
				workGrid(deltaVIndices + interfaceGrid.forward + interfaceGrid.right),
				workGrid(deltaVIndices + interfaceGrid.forward + interfaceGrid.top),
				workGrid(deltaVIndices + interfaceGrid.forward + interfaceGrid.right + interfaceGrid.top), dummyarea);
			double changeInVol = (newVolume - oldVolume) / (tau * ns * ns * ns);

			//if (changeInVol < 0) {
				//cerr << "Tell me what happened: Cell " << I << " " << J << " " << K << "\n";
			//}

			whichElement.addContribution(index, changeInVol);
		} // Done with changed volume

		// not quite the right fix
		if (whichElement.getContribution(middleCell) < 0.0) {
			std::cerr << "Negative: " << I << " " << J << " " << K << " (num)\n";
			//			double fix = surfaceAreas(Material, coarseInds) / 6.0 - whichElement.getContribution(middleCell);
			//			whichElement.addContribution(middleCell, fix);
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
	int * matrixRow = new int[(coarseM + 2)*(coarseN + 2)*(coarseP + 2)];
	double * A, *rhs;

	int rows = 0;

	for (int inds = 0; inds < coarseSize; inds++) {
		lambdas(inds) = 0.0;
		residual(inds) = grid(Material, inds) - vActs(Material, inds);
		double diag = fabs(volumeOfGradients(inds).getContribution(middleCell));

		if (surfaceAreas(Material, inds) >= 1e-14 && diag > 0) {
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
			IR = indR / ((coarseN + 2) * (coarseP + 2)),
			JR = (indR / (coarseP + 2)) % (coarseN + 2),
			KR = indR % (coarseP + 2);

		for (int c = 0; c < rows; c++) {
			int indC = matrixRow[c],
				IC = indC / ((coarseN + 2) * (coarseP + 2)),
				JC = (indC / (coarseP + 2)) % (coarseN + 2),
				KC = indC % (coarseP + 2);

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
			IR = indR / ((coarseN + 2) * (coarseP + 2)),
			JR = (indR / (coarseP + 2)) % (coarseN + 2),
			KR = indR % (coarseP + 2);
		//const double maxSpeed = 0.1 * hx;//0.02 * hx;
		double maxSpeed = surfaceAreas(Material, indR) * 0.25 * hx;
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
		if (diag < 0.0 && surfaceAreas(Material, inds) > 1e-8) { std::cerr << "Diag = " << diag << "\n"; }

		lambdas(inds) = 0.0;
		double r = grid(Material, inds) - vActs(Material, inds);
		residual(inds) = r;
		L1Res += fabs(r);
	}

	std::cerr << "Residual starts at " << L1Res << endl;

	const int right = residual.right, top = residual.top, forward = residual.forward;

	for (int step = 0; step < numIters; step++) {

		for (int inds = 0; inds < coarseSize; inds++) {
			double diag = volumeOfGradients(inds).getContribution(middleCell);
			if (surfaceAreas(Material, inds) <= 1e-8 || fabs(diag) <= 1e-8) {
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

				if (surfaceAreas(Material, neighbor) <= 1e-8) {
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
		double area = surfaceAreas(Material, inds);
		if (area < 1e-8) { continue; }

		double maxSpeed = 0.1 * hx;

		if (lambdas(inds) > maxSpeed) { lambdas(inds) = maxSpeed; }
		else if (lambdas(inds) < -maxSpeed) { lambdas(inds) = -maxSpeed; }
	}
}

// Iterates over material
void Interface3D::checkResult(bool verbose, bool voronoiVol) {

	double maxError = 0.0, totalRequiredVol = 0.0, totalComputedVol = 0.0,
		sumOfErrors = 0.0;

	if (!voronoiVol) {
		computeAllVActs();
	}
	else {
		computeVActsExactMulti();
	}
	//adjustVolumes(vActs);

	int mixed = 0;

	for (Index3 inds(1, 1, 1, coarseM, coarseN, coarseP, coarseM + 1, coarseN + 1, coarseP + 1); inds.valid(); ++inds) {
		bool isMixed = false;
		for (int mat = 0; mat < numPhases; mat++) {
			if (grid(mat, inds) > 1e-8 && grid(mat, inds) < 1 - 1e-8) {
				isMixed = true;
				break;
			}
		}

		if (isMixed) {
			mixed++;
			if (verbose) { cout << "Cell (" << inds.i << "," << inds.j << "," << inds.k << "): "; }
			double maxOverMats = 0.0;
			for (int mat = 0; mat < numPhases; mat++) {
				double vAct = vActs(mat, inds);
				double vReq = grid(mat, inds);
				totalComputedVol += vAct;
				totalRequiredVol += vReq;

				double error = fabs(vReq - vAct);

				if (error > maxOverMats) { maxOverMats = error; }
				if (error > maxError) { maxError = error; }
				if (verbose) { cout << "Required " << vReq << "; actual " << vAct << endl; }
				//else if (error > 0.01) {
				//   cout << "1%+ error in cell (" << I << " " << J << " " << K << "): Required " << vReq << ", error " << error << endl;
				//}
			}

			sumOfErrors += maxOverMats;
		} // if mixed
	}

	// Count faces
	/*int faces = 0;
	for (int i = 0; i < m-1; i++) {
		for (int j = 0; j < n-1; j++) {
			for (int k = 0; k < p-1; k++) {
				// (0,0,0), (0,1,0), (1,1,0), (0,0,1)
				double dummyArea = 0.0;
				if (fluidTetrahedron(interfaceGrid[i][j][k], interfaceGrid[i][j+1][k],
									 interfaceGrid[i+1][j+1][k], interfaceGrid[i][j][k+1], dummyArea) > 1e-9) {
					faces++;
				}
				if(fluidTetrahedron(interfaceGrid[i][j+1][k+1], interfaceGrid[i][j+1][k],
									 interfaceGrid[i+1][j+1][k], interfaceGrid[i][j][k+1], dummyArea) > 1e-9) {
					faces++;
				}
				if(fluidTetrahedron(interfaceGrid[i+1][j+1][k+1], interfaceGrid[i][j+1][k+1],
									 interfaceGrid[i][j][k+1], interfaceGrid[i+1][j+1][k], dummyArea) > 1e-9) {
					faces++;
				}
				if(fluidTetrahedron(interfaceGrid[i+1][j+1][k+1], interfaceGrid[i+1][j][k+1],
									 interfaceGrid[i][j][k+1], interfaceGrid[i+1][j+1][k], dummyArea) > 1e-9) {
					faces++;
				}
				if (fluidTetrahedron(interfaceGrid[i+1][j][k], interfaceGrid[i+1][j][k+1],
									 interfaceGrid[i][j][k+1], interfaceGrid[i+1][j+1][k], dummyArea) > 1e-9) {
					faces++;
				}
				if(fluidTetrahedron(interfaceGrid[i][j][k], interfaceGrid[i+1][j][k],
									 interfaceGrid[i+1][j+1][k], interfaceGrid[i][j][k+1], dummyArea) > 1e-9) {
					faces++;
				}
			}
		}
	}*/


	cout << "Max error = " << maxError << "; L1 error = " << sumOfErrors << "\nTotal volume in interface = "
		<< totalComputedVol << "; required vol = " << totalRequiredVol << endl;

	cerr << "Mixed " << ((double)mixed / (coarseM * coarseN * coarseP)) << "\n"; //<< " | Faces " << faces << endl;
	cerr << "Average error = " << (sumOfErrors / mixed) << endl;
}

double Interface3D::calculateVolumeImmersed(
	int I, int J, int K) {
	double vol = 0.0, area = 0.0;
	int span = ns;
	double unitVol = 1.0 / (span*span*span), unitArea = 1.0 / (span*span);

	int startIndexX = (I - 1)*ns, startIndexY = (J - 1)*ns, startIndexZ = (K - 1)*ns;

	for (Index3 inds(startIndexX, startIndexY, startIndexZ, startIndexX + span - 1, startIndexY + span - 1, startIndexZ + span - 1,
		m, n, p); inds.valid(); ++inds) {

		const int center = (int)inds, right = interfaceGrid.right, top = interfaceGrid.top, forward = interfaceGrid.forward;
		// 000, 100, 010, 110, 001, 101, 011, 111
		vol += fluidAndAreaBox(interfaceGrid(Material, center), interfaceGrid(Material, center + right),
			interfaceGrid(Material, center + top), interfaceGrid(Material, center + right + top),
			interfaceGrid(Material, center + forward), interfaceGrid(Material, center + right + forward),
			interfaceGrid(Material, center + top + forward),
			interfaceGrid(Material, center + right + top + forward),
			area);
	}

	vol *= unitVol;
	area *= unitArea;

	int IJK = indexFromCoarse(I, J, K);

	double change = fabs(vol - vActs(Material, IJK));

	vActs(Material, IJK) = vol;
	// Note I think we need surface areas by material too
	surfaceAreas(Material, IJK) = area;

	totalSurfaceArea += area;

	return change;
}

Interface3D::~Interface3D() {
	// Delete all that memory we allocated. Automatic with our grid classes;
   // but delete narrow band.
	for (Material = 0; Material < numPhases; Material++) {
		delete[] narrowBandIterator[Material];
	}
	delete[] narrowBandIterator;
	delete[] narrowBandSize;
}

void Interface3D::saveInterface() {
	// For each node in the narrow band
	for (int mat = 0; mat < numPhases; mat++) {
		int nbSize = narrowBandSize[mat];
		for (int nbInd = 0; nbInd < nbSize; nbInd++) {
			int index = narrowBandIterator[mat][nbInd].Index();
			oldInterface(mat, index) = interfaceGrid(mat, index);
		}
	}
}

double Interface3D::interfaceChange() {
	double change = 0.0;

	// For each node in the narrow band
	for (int mat = 0; mat < numPhases; mat++) {
		int nbSize = narrowBandSize[mat];
		for (int nbInd = 0; nbInd < nbSize; nbInd++) {
			int index = narrowBandIterator[mat][nbInd].Index();

			double phi = interfaceGrid(mat, index),
				psi = oldInterface(mat, index);
			change += fabs(phi - psi) * exp(-0.125 * psi * psi / (hx*hx));
		}
	}

	return change * hx * hy * hz;
}
