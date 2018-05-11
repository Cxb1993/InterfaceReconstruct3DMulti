#include "stdafx.h"
#include "Interface3D.h"

using namespace std;

// Bugs to not forget (even though they disappeared)
// 1. grid had a value < 0
void Interface3D::computeVActsExactMulti() { // Also computes vActs, but computes the exact triangularization and finds its volume.
	// First, let's do this not-so efficiently, but just quickly

	// First, zero out the computed volumes.
	for (int mat = 0; mat < numPhases; mat++) {
		for (int ind = 0; ind < coarseSize; ind++) {
			if (coarseNarrowBand(mat, ind)) { 
				vActsBack(mat, ind) = vActs(mat, ind);
				vActs(mat, ind) = 0.0; 
			}
		}
	}

	double (* verticesOut)[3], (* verticesIn)[3];
	int ** facesOut, ** facesIn, nVertsIn, nFacesIn, nVertsOut, nFacesOut;
	// Please init all arrays, foo'.
	facesOut = new int*[50]; verticesOut = new double[50][3];
	facesIn  = new int*[50]; verticesIn  = new double[50][3];
	for (int i = 0; i < 50; i++) {
	    facesOut[i] = new int[50];
		facesIn[i]  = new int[50];
	}
	
	const int right = interfaceGrid.right, top = interfaceGrid.top, forward = interfaceGrid.forward;

	int tetIndices[] = {
		0, top, right + top, forward, // (0,0,0), (0,1,0), (1,1,0), (0,0,1)
		top + forward, top, right + top, forward, // (0,1,1), (0,1,0), (1,1,0), (0,0,1)
		right + top + forward, top + forward, forward, right + top, // (1,1,1), (0,1,1), (0,0,1), (1,1,0)
		right + top + forward, right + forward, forward, right + top, // (1,1,1), (1,0,1), (0,0,1), (1,1,0)
		right, right + forward, forward, right + top, // (1,0,0), (1,0,1), (0,0,1), (1,1,0)
		0, right, right + top, forward // (0,0,0), (1,0,0), (1,1,0), (0,0,1)
	};

	for (Index3 inds(0,0,0, m-1,n-1,p-1, m,n,p); inds.valid(); ++inds) {
		for (int MatA = 0; MatA < numPhases; MatA++) {
			//if (!narrowBand(MatA, inds)) { continue; }
			int coarseInd = coarseIndexFromFine(inds.i, inds.j, inds.k);
			if (!coarseNarrowBand(MatA, coarseInd)) { continue; }
			// For each tetrahedron
			for (int T = 0; T < 6; T++) {
			    // Start with the whole tetrahedron.
				//if (inds.ind == 194445 && T == 0 && MatA == 0) {
				//	cerr << "Also bp\n";
				//}
			    for (int i = 0; i < 4; i++) {
				    verticesIn[i][0] = vertices[i][0]; // vertices is the unit tetrahedron. Should refactor it more important-looking.
				    verticesIn[i][1] = vertices[i][1];
				    verticesIn[i][2] = vertices[i][2];
				    for (int j = 0; j < 5; j++) {
					    facesIn[i][j] = unitTetFaces[i][j];
				    }
			    }
			    nVertsIn = nFacesIn = 4;

			    for (int MatB = 0; MatB < numPhases; MatB++) {
					if (MatB == MatA || !coarseNarrowBand(MatB, coarseInd)) { continue; }

					double a = interfaceGrid(MatA, inds + tetIndices[4*T]) - interfaceGrid(MatB, inds + tetIndices[4*T]),
						   b = interfaceGrid(MatA, inds + tetIndices[4*T+1]) - interfaceGrid(MatB, inds + tetIndices[4*T+1]),
						   c = interfaceGrid(MatA, inds + tetIndices[4*T+2]) - interfaceGrid(MatB, inds + tetIndices[4*T+2]),
						   d = interfaceGrid(MatA, inds + tetIndices[4*T+3]) - interfaceGrid(MatB, inds + tetIndices[4*T+3]);

					int nMinus = 0;
					if (a <= 0.0) { nMinus++; }
					if (b <= 0.0) { nMinus++; }
					if (c <= 0.0) { nMinus++; }
					if (d <= 0.0) { nMinus++; }
					if (nMinus == 0) { 
						nFacesIn = 0; break; 
					} else if (nMinus == 4) {
						continue;
					}

				    chopPolytope(verticesIn, nVertsIn, facesIn, nFacesIn, verticesOut, nVertsOut, facesOut, nFacesOut, c-b, b-a, d-a, a);

					// Move data back over to input.
					for (int i = 0; i < nVertsOut; i++) {
						verticesIn[i][0] = verticesOut[i][0];
						verticesIn[i][1] = verticesOut[i][1];
						verticesIn[i][2] = verticesOut[i][2];
					}
					for (int i = 0; i < nFacesOut; i++) {
						for (int j = 0; j < 50; j++) {
						    facesIn[i][j] = facesOut[i][j];
						}
					}
					nVertsIn = nVertsOut;
					nFacesIn = nFacesOut;
			    } // For each other material

				if (nFacesIn > 0) {
					double volume = volumeOfPolytope(verticesIn, nVertsIn, facesIn, nFacesIn);
					if (volume < 1.0/6.0 - 1e-12) {
						//cerr << "BP\n";
					}
					vActs(MatA, coarseInd) += volume / (ns*ns*ns);
				}
			} // For each Tetrahedron
		} // For each material
	} // For each fine cube

	// Delete
	for (int i = 0; i < 50; i++) {
	    delete[] facesOut[i];
		delete[] facesIn[i];
	}

	delete[] facesOut;
	delete[] facesIn;
	delete[] verticesOut;
	delete[] verticesIn;
}

void Interface3D::computeSingleExactVAct(int MatA) {

	// First, zero out the computed volumes.
		for (int ind = 0; ind < coarseSize; ind++) {
			if (coarseNarrowBand(MatA, ind)) { 
				vActsBack(MatA, ind) = vActs(MatA, ind);
				vActs(MatA, ind) = 0.0; 
			}
		}

	double (* verticesOut)[3], (* verticesIn)[3];
	int ** facesOut, ** facesIn, nVertsIn, nFacesIn, nVertsOut, nFacesOut;
	// Please init all arrays, foo'.
	facesOut = new int*[50]; verticesOut = new double[50][3];
	facesIn  = new int*[50]; verticesIn  = new double[50][3];
	for (int i = 0; i < 50; i++) {
	    facesOut[i] = new int[50];
		facesIn[i]  = new int[50];
	}
	
	const int right = interfaceGrid.right, top = interfaceGrid.top, forward = interfaceGrid.forward;

	int tetIndices[] = {
		0, top, right + top, forward, // (0,0,0), (0,1,0), (1,1,0), (0,0,1)
		top + forward, top, right + top, forward, // (0,1,1), (0,1,0), (1,1,0), (0,0,1)
		right + top + forward, top + forward, forward, right + top, // (1,1,1), (0,1,1), (0,0,1), (1,1,0)
		right + top + forward, right + forward, forward, right + top, // (1,1,1), (1,0,1), (0,0,1), (1,1,0)
		right, right + forward, forward, right + top, // (1,0,0), (1,0,1), (0,0,1), (1,1,0)
		0, right, right + top, forward // (0,0,0), (1,0,0), (1,1,0), (0,0,1)
	};

	for (Index3 coarseInds(1,1,1, coarseM,coarseN,coarseP, coarseM+1,coarseN+1,coarseP+1); coarseInds.valid(); ++coarseInds) {
		if (!coarseNarrowBand(MatA, coarseInds)) { continue; }
		int iStart = (coarseInds.i-1)*ns, jStart = (coarseInds.j-1)*ns, kStart = (coarseInds.k-1)*ns;
		for (Index3 inds(iStart, jStart, kStart, iStart + ns-1, jStart + ns-1, kStart + ns-1, m, n, p); inds.valid(); ++inds) {
			if(interfaceGrid(MatA, inds) > ns * hx) {
				continue;
			} else if (interfaceGrid(MatA, inds) < -ns * hx) {
				vActs(MatA, coarseInds) += 1.0 / (ns * ns * ns);
				continue;
			}
		    // For each tetrahedron
			for (int T = 0; T < 6; T++) {
			    // Start with the whole tetrahedron.
				//if (inds.ind == 194445 && T == 0 && MatA == 0) {
				//	cerr << "Also bp\n";
				//}
			    for (int i = 0; i < 4; i++) {
				    verticesIn[i][0] = vertices[i][0]; // vertices is the unit tetrahedron. Should refactor it more important-looking.
				    verticesIn[i][1] = vertices[i][1];
				    verticesIn[i][2] = vertices[i][2];
				    for (int j = 0; j < 5; j++) {
					    facesIn[i][j] = unitTetFaces[i][j];
				    }
			    }
			    nVertsIn = nFacesIn = 4;

			    for (int MatB = 0; MatB < numPhases; MatB++) {
					if (MatB == MatA || !coarseNarrowBand(MatB, coarseInds)) { continue; }

					double a = interfaceGrid(MatA, inds + tetIndices[4*T]) - interfaceGrid(MatB, inds + tetIndices[4*T]),
						   b = interfaceGrid(MatA, inds + tetIndices[4*T+1]) - interfaceGrid(MatB, inds + tetIndices[4*T+1]),
						   c = interfaceGrid(MatA, inds + tetIndices[4*T+2]) - interfaceGrid(MatB, inds + tetIndices[4*T+2]),
						   d = interfaceGrid(MatA, inds + tetIndices[4*T+3]) - interfaceGrid(MatB, inds + tetIndices[4*T+3]);

					int nMinus = 0;
					if (a <= 0.0) { nMinus++; }
					if (b <= 0.0) { nMinus++; }
					if (c <= 0.0) { nMinus++; }
					if (d <= 0.0) { nMinus++; }
					if (nMinus == 0) { 
						nFacesIn = 0; break; 
					} else if (nMinus == 4) {
						continue;
					}

				    chopPolytope(verticesIn, nVertsIn, facesIn, nFacesIn, verticesOut, nVertsOut, facesOut, nFacesOut, c-b, b-a, d-a, a);

					// Move data back over to input.
					for (int i = 0; i < nVertsOut; i++) {
						verticesIn[i][0] = verticesOut[i][0];
						verticesIn[i][1] = verticesOut[i][1];
						verticesIn[i][2] = verticesOut[i][2];
					}
					for (int i = 0; i < nFacesOut; i++) {
						for (int j = 0; j < 50; j++) {
						    facesIn[i][j] = facesOut[i][j];
						}
					}
					nVertsIn = nVertsOut;
					nFacesIn = nFacesOut;
			    } // For each other material

				if (nFacesIn > 0) {
					double volume = volumeOfPolytope(verticesIn, nVertsIn, facesIn, nFacesIn);
					vActs(MatA, coarseInds) += volume / (ns*ns*ns);
				}
			} // For each Tetrahedron
		} // For each material
	} // For each fine cube

	// Delete
	for (int i = 0; i < 50; i++) {
	    delete[] facesOut[i];
		delete[] facesIn[i];
	}

	delete[] facesOut;
	delete[] facesIn;
	delete[] verticesOut;
	delete[] verticesIn;
}

void Interface3D::addPolygonsForTetrahedron(std::vector<double> & points, std::vector<int> & polySizes, std::vector<int> & faceTypes, 
	                                        std::vector<double> & errorVals, int indices[4][3], int transform[3], int coarseInd) {
	// For each pair of Materials
	for (int MatI = 0; MatI < numPhases - 1; MatI++) {
		for (int MatJ = MatI + 1; MatJ < numPhases; MatJ++) {
			double phiIa = interfaceGrid(MatI, indices[0][0], indices[0][1], indices[0][2]),
				   phiJa = interfaceGrid(MatJ, indices[0][0], indices[0][1], indices[0][2]),
				   phiIb = interfaceGrid(MatI, indices[1][0], indices[1][1], indices[1][2]),
				   phiJb = interfaceGrid(MatJ, indices[1][0], indices[1][1], indices[1][2]),
				   phiIc = interfaceGrid(MatI, indices[2][0], indices[2][1], indices[2][2]),
				   phiJc = interfaceGrid(MatJ, indices[2][0], indices[2][1], indices[2][2]),
				   phiId = interfaceGrid(MatI, indices[3][0], indices[3][1], indices[3][2]),
				   phiJd = interfaceGrid(MatJ, indices[3][0], indices[3][1], indices[3][2]);

			if (fabs(phiIa) > ns * hx * 1.5 || fabs(phiJa) > ns * hx * 1.5) {
				continue; 
			}

			double a = phiIa - phiJa,
				   b = phiIb - phiJb,
				   c = phiIc - phiJc,
				   d = phiId - phiJd;

			int nBase;
			double base[10][3], baseTmp[10][3];
			getInterfacePolygon(a, b, c, d, nBase, base);
			bool hasAnything = true;

			// Hackjob
			for (int v = 0; v < nBase; v++) {
				for (int x = 0; x < 3; x++) {
					if (isNaN(base[v][x])) {
						hasAnything = false;
					}
				}
			}
			if (!hasAnything) { continue; }

			for (int MatK = 0; MatK < numPhases; MatK++) {
				if (MatK == MatI || MatK == MatJ) { continue; }
				// If MatC at any of those base points is smaller than MatA at any of them...
				// Then trim!
                double phiKa = interfaceGrid(MatK, indices[0][0], indices[0][1], indices[0][2]), 
				       phiKb = interfaceGrid(MatK, indices[1][0], indices[1][1], indices[1][2]),
				       phiKc = interfaceGrid(MatK, indices[2][0], indices[2][1], indices[2][2]),
				       phiKd = interfaceGrid(MatK, indices[3][0], indices[3][1], indices[3][2]);
				// Interpolate onto base points
				// phiK(x,y,z) = phiKa + (phiKb - phiKa) y + (phiKc - phiKb) x + (phiKd - phiKa) z
				double phiDiffAtBase[10];
				bool trim;
				int positivePt = -1;
				for (int c = 0; c < nBase; c++) {
					phiDiffAtBase[c] = (phiKa + (phiKb - phiKa) * base[c][1] + (phiKc - phiKb) * base[c][0] + (phiKd - phiKa) * base[c][2]) -
						               (phiIa + (phiIb - phiIa) * base[c][1] + (phiIc - phiIb) * base[c][0] + (phiId - phiIa) * base[c][2]);
					if (phiDiffAtBase[c] < 0.0) { trim = true; } else if (positivePt < 0) { positivePt = c; }
				}

				if (trim) {
					if (positivePt < 0) {
						// Nothing, yo.
						hasAnything = false;
						break;
					}
					// Do linear interpolation along edges to find points when phiK = phiI.
					// Algorithm: Start with a known point where phiI,J < phiK. Walk along polygon until an edge has an intersection point.
					// Add that to new base. Then...
					baseTmp[0][0] = base[positivePt][0];
					baseTmp[0][1] = base[positivePt][1];
					baseTmp[0][2] = base[positivePt][2];
					int numNewPoly = 1;
					int c = positivePt;
					do {
						int next = c + 1;
						if (next >= nBase) { next = 0; }

						double t = -phiDiffAtBase[c] / (phiDiffAtBase[next] - phiDiffAtBase[c]);
						if (t > 0.0 && t < 1.0) {
							baseTmp[numNewPoly][0] = base[c][0] + t * (base[next][0] - base[c][0]);
							baseTmp[numNewPoly][1] = base[c][1] + t * (base[next][1] - base[c][1]);
							baseTmp[numNewPoly][2] = base[c][2] + t * (base[next][2] - base[c][2]);
							numNewPoly++;
						} 

						if (phiDiffAtBase[next] >= 0.0 && next != positivePt) {
							baseTmp[numNewPoly][0] = base[next][0];
							baseTmp[numNewPoly][1] = base[next][1];
							baseTmp[numNewPoly][2] = base[next][2];
							numNewPoly++;
						}

						c = next;
					} while (c != positivePt);

					// need to swap base and baseTmp
					// wanted a shallow copy but a deep one will do for now
					for (int c = 0; c < numNewPoly; c++) {
						for (int c2 = 0; c2 < 3; c2++) {
							base[c][c2] = baseTmp[c][c2];
						}
					}
					nBase = numNewPoly;
				}
				// Now, onto next Material to see if any more intersections.

			} // Iterate over MatK

			if (!hasAnything || nBase == 0) {
				continue;
			}

			// Fix front/back face. Important for rendering.
			double orientx1 = base[1][0] - base[0][0], orienty1 = base[1][1] - base[0][1], orientz1 = base[1][2] - base[0][2],
				   orientx2 = base[2][0] - base[0][0], orienty2 = base[2][1] - base[0][1], orientz2 = base[2][2] - base[0][2];
			double crossX = orienty1 * orientz2 - orientz1 * orienty2, 
				   crossY = orientz1 * orientx2 - orientx1 * orientz2,
				   crossZ = orientx1 * orienty2 - orienty1 * orientx2;
			// Now define the front face as the one where MatA < MatB.
			// Normal vector for phiA - phiB is (c-b, b-a, d-a).
			// No need to adjust for transform because
			// T(n) dot T(xprod) = T(n dot xprod)
			double orientation = crossX * (c-b) + crossY * (b-a) + crossZ * (d-a);

			int cStart = 0, cStep = 1, cEnd = nBase;
			if (orientation <= 0.0) {
				cStart = nBase - 1;
				cStep = -1;
				cEnd = -1;
			}			

			// Resulting polygon passes test. Add it to list.
			for (int c = cStart; c != cEnd; c += cStep) {
				// Oh no, you need to transform.
				double transformed[3];
				computeCoord(indices[0], transform, base[c], transformed);
				points.push_back(transformed[0]);
				points.push_back(transformed[1]);
				points.push_back(transformed[2]);

				if (isNaN(transformed[0]) || isNaN(transformed[1]) || isNaN(transformed[2])) {
					cerr << "Well\n";
				}
			}
			polySizes.push_back(nBase);
			faceTypes.push_back(MatI * numPhases + MatJ);
			errorVals.push_back((grid(MatI, coarseInd) - vActs(MatI, coarseInd)));
		} // MatJ
	} // MatI
}

// For the tetrahedron defined by (0,0,0), (0,1,0), (1,1,0), (0,0,1), assuming f(0,0,0) = a, f(0,1,0) = b, etc,
// and that f = c_0 + c_1 x + c_2 y + c_3 z, finds the plane surface f(x,y,z) = 0.
void Interface3D::getInterfacePolygon(double a, double b, double c, double d, int & nBase, double base[4][3]) {

    double nMinus = 0;
    if (a <= 0) { nMinus++; }
    if (b <= 0) { nMinus++; }
    if (c <= 0) { nMinus++; }
    if (d <= 0) { nMinus++; }
           
    if (nMinus == 0 || nMinus == 4) { nBase = 0; return; }

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
    int e0;
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
    nBase = 3;
    if (edges[oppositeEdge[e0]]) {                 
        nBase = 4;
        base[3][0] = v[oppositeEdge[e0]][0];
        base[3][1] = v[oppositeEdge[e0]][1];
        base[3][2] = v[oppositeEdge[e0]][2];    
    }
}

void Interface3D::chopPolytope(double Verts[][3], int numVerts, int ** faces, int nFaceIn, 
	                           double newVerts[][3], int & numNewVerts, int ** newFaces, int & nFaceOut,
	                           double a, double b, double c, double d) {
	double * funcValues = new double[numVerts];
	for (int i = 0; i < numVerts; i++) {
		funcValues[i] = a * Verts[i][0] + b * Verts[i][1] + c * Verts[i][2] + d;
	}

	// Edges identified by which two vertices they connect.
	// Technically only the upper triangle of this matrix is used.
	int * edgeToNewVertIndex = new int[numVerts * numVerts];
	for (int i = 0; i < numVerts*numVerts; i++) {
		edgeToNewVertIndex[i] = 0;
	}

	int * newVerticesList = new int[2 * nFaceIn];
	for (int i = 0; i < 2 * nFaceIn; i++) {
		newVerticesList[i] = 0;
	}

	// Assume newVerts. newFaces already allocated.
	// Copy all old vertices, even if we don't need them. What's the harm?
	for (int v = 0; v < numVerts; v++) {
		newVerts[v][0] = Verts[v][0];
		newVerts[v][1] = Verts[v][1];
		newVerts[v][2] = Verts[v][2];
	}
	numNewVerts = numVerts;

	// faces is specified as follows: faces[f]: list of vertices from Verts, in order, and having first and last as the same index.
	// Also, end the list with -1. (Can probably work around that by checking if same as start, but whatever)

	int numChangedFaces = 0;
	int numNewFaces = 0;

	// For each face
	for (int f = 0; f < nFaceIn; f++) {
		// Make new face out of chopped one, start with 0 vertices
		int numVInNewFace = 0;
		int numNewVUsed = 0;
		// For each vertex pair
		for (int e = 1; faces[f][e] >= 0; e++) {
			// The vertex pair (edge)
			int v0 = faces[f][e-1], v1 = faces[f][e];
			int whichEdge = min(v0,v1) * numVerts + max(v0,v1);
			double valA = funcValues[v0], valB = funcValues[v1];
			// Survived chop
			if (valA <= 0.0) {
				newFaces[numNewFaces][numVInNewFace++] = v0;
			}
			double tau = -valA / (valB - valA);
			// New chop point
			if (tau > 0.0 && tau < 1.0) {
				// Note that there is a vertex between two other vertices now because of the chopping
				int newVertIndex = edgeToNewVertIndex[whichEdge];
				if (newVertIndex == 0) {
					edgeToNewVertIndex[whichEdge] = newVertIndex = numNewVerts;
					newVerts[newVertIndex][0] = (1.0 - tau) * Verts[v0][0] + tau * Verts[v1][0];
					newVerts[newVertIndex][1] = (1.0 - tau) * Verts[v0][1] + tau * Verts[v1][1];
					newVerts[newVertIndex][2] = (1.0 - tau) * Verts[v0][2] + tau * Verts[v1][2];
					++numNewVerts;
				}
				// Add new vertex, since its value is zero it is part of the new polytope
				newFaces[numNewFaces][numVInNewFace++] = newVertIndex;

				newVerticesList[2 * f + numNewVUsed++] = newVertIndex;
			}
		} // New face has been processed
		
		if (numNewVUsed > 0) { numChangedFaces++; }

		if (numVInNewFace > 0) {
			newFaces[numNewFaces][numVInNewFace] = newFaces[numNewFaces][0];
			newFaces[numNewFaces][numVInNewFace+1] = -1;
			numNewFaces++;
		} // This condition is basically only so, if nothing in new face, don't add a new face.
	} // Done processing new faces.

	// Now that all old faces were chopped, we must add the final face, which is on the plane ax + by + cz + d = 0.
	// Any new edges that were created are part of this new face. Connect them in order.
	// Then check the normal at the end to make sure you return it counter-clockwise.

	// Ok, do the first two vertices in new face

	if (numChangedFaces > 0) {
    	int prevV;
    	for (int nv = 0; nv < nFaceIn; nv++) {
	    	int whatV = newVerticesList[2*nv];
		    if (newVerticesList[2*nv] > 0) {
				prevV = newVerticesList[2*nv+1];
	    		newFaces[numNewFaces][0] = whatV;
		    	newFaces[numNewFaces][1] = newVerticesList[2*nv+1];

				// Zero them out
				newVerticesList[2*nv] = newVerticesList[2*nv+1] = 0;
				break;
		    }
	    }

	    for (int currentV = 2; currentV < numChangedFaces; currentV++) {
			for (int nv = 0; nv < nFaceIn; nv++) {
				// Find the matching vertex
				int vA = newVerticesList[2*nv], vB = newVerticesList[2*nv+1];
				if (vA == prevV) {
					newFaces[numNewFaces][currentV] = vB;
					prevV = vB;
					newVerticesList[2*nv] = newVerticesList[2*nv+1] = 0;
					break;
				} else if (vB == prevV) {
					newFaces[numNewFaces][currentV] = vA;
					prevV = vA;
					newVerticesList[2*nv] = newVerticesList[2*nv+1] = 0;
					break;
				}
			}
	    }

		newFaces[numNewFaces][numChangedFaces] = newFaces[numNewFaces][0];
		newFaces[numNewFaces][numChangedFaces+1] = -1;

		// Ok, now find out if the orientation is correct or not. Is easy; the outer normal is (a,b,c)
		// and the normal to the face is (v1-v0)x(v2-v0).
		// Fix front/back face. Important for rendering.
		double orientx1 = newVerts[newFaces[numNewFaces][1]][0] - newVerts[newFaces[numNewFaces][0]][0], 
			   orienty1 = newVerts[newFaces[numNewFaces][1]][1] - newVerts[newFaces[numNewFaces][0]][1], 
			   orientz1 = newVerts[newFaces[numNewFaces][1]][2] - newVerts[newFaces[numNewFaces][0]][2], 
			   orientx2 = newVerts[newFaces[numNewFaces][2]][0] - newVerts[newFaces[numNewFaces][0]][0], 
			   orienty2 = newVerts[newFaces[numNewFaces][2]][1] - newVerts[newFaces[numNewFaces][0]][1], 
			   orientz2 = newVerts[newFaces[numNewFaces][2]][2] - newVerts[newFaces[numNewFaces][0]][2];
		double crossX = orienty1 * orientz2 - orientz1 * orienty2, 
			   crossY = orientz1 * orientx2 - orientx1 * orientz2,
			   crossZ = orientx1 * orienty2 - orienty1 * orientx2;
		if (crossX * a + crossY * b + crossZ * c < 0) {
			// Reverse it
			for (int i = 1; i <= numChangedFaces / 2; i++) {
				int tmp = newFaces[numNewFaces][i];
				newFaces[numNewFaces][i] = newFaces[numNewFaces][numChangedFaces - i];
				newFaces[numNewFaces][numChangedFaces - i] = tmp;
			}
		}

		numNewFaces++;
	}

	nFaceOut = numNewFaces;
	// Ok... I think that's it.
	delete[] funcValues;
	delete[] edgeToNewVertIndex;
	delete[] newVerticesList;
}

double Interface3D::volumeOfPolytope(const double Vertices[][3], int numVerts, int ** faces, int numFaces) {
	// Volume is calculated as follows:
	// Integral over P of 1 dx dy dz =
	// Integral over P of 1/3 div(x, y, z) dx dy dz =
	// Integral over dP of 1/3 (x,y,z) dot n dS =
	// Sum_i of 1/3 area(face_i) centroid(face_i) dot normal(face_i)
	double volume = 0.0;
	for (int i = 0; i < numFaces; i++) {
		double area[3], centroid[3];
		
		// Any uninitialized coordinates will be the same as the first vertex; this is
		// for those exceptional cases where a face is exactly on x = 0, y = 0, or z = 0.
		int firstV = faces[i][0];
		centroid[0] = Vertices[firstV][0];
		centroid[1] = Vertices[firstV][1];
		centroid[2] = Vertices[firstV][2];
		// Pick two coordinates (two of x,y,z). For example, x,z.
		// Compute centroid of projection (x,y,z) -> (x,z)
		// That gives xbar and zbar. Then do (x,y,z) -> (x,y).
		// Then find ybar.
		// If Area = 0, then try another combo.
		for (int coordPair = 0; coordPair < 3; coordPair++) {
			int first = (coordPair+1)%3, second = (coordPair+2)%3;
			// Perform calculation.
			double areaTmp = 0.0; 
			double centroid_first = 0.0, centroid_second = 0.0;
			for (int j = 1; faces[i][j] > 0; j++) {
				int v0 = faces[i][j-1], v1 = faces[i][j];
				double u_cur = Vertices[v0][first], u_next = Vertices[v1][first],
					   v_cur = Vertices[v0][second], v_next = Vertices[v1][second];
				double areaDelta = u_cur * v_next - u_next * v_cur;
				areaTmp += areaDelta;
				centroid_first += (u_cur + u_next) * areaDelta;
				centroid_second += (v_cur + v_next) * areaDelta;
			}

			areaTmp *= 0.5;
			area[coordPair] = areaTmp;
			if (fabs(areaTmp) > 1e-10) { 	
				centroid_first /= (6.0 * areaTmp);
				centroid_second /= (6.0 * areaTmp);
				centroid[first]  = centroid_first;
				centroid[second] = centroid_second;
			}
		}

		//double Area = sqrt(area[0]*area[0] + area[1]*area[1] + area[2]*area[2]);
		double volContribution = area[0] * centroid[0] + area[1] * centroid[1] + area[2] * centroid[2];
		volume += volContribution;
	}

	return volume / 3.0;
}


void Interface3D::normalToInterface(int mat, int i, int j, int k,
                                        double & n1, double & n2, double & n3) {
     int numVectors = 0;
     n1 = 0.0; n2 = 0.0; n3 = 0.0;
     double n1t, n2t, n3t;

	 int ind000 = indexFromFine(i, j, k), 
		 ind100 = ind000 + interfaceGrid.right,
		 ind010 = ind000 + interfaceGrid.top,
		 ind110 = ind100 + interfaceGrid.top,
		 ind001 = ind000 + interfaceGrid.forward,
		 ind101 = ind100 + interfaceGrid.forward,
		 ind011 = ind010 + interfaceGrid.forward,
		 ind111 = ind110 + interfaceGrid.forward;
           
     double A, B, C, D;
     // (0,0,0), (0,1,0), (1,1,0), (0,0,1)
     A = interfaceGrid(mat, ind000); B = interfaceGrid(mat, ind010); 
     C = interfaceGrid(mat, ind110); D = interfaceGrid(mat, ind001);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = C-B;
        n2t = B-A;
        n3t = D-A;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }           
     // (1,1,1), (0,1,1), (0,0,1), (1,1,0)
     A = interfaceGrid(mat, ind111); B = interfaceGrid(mat, ind011); 
     C = interfaceGrid(mat, ind001); D = interfaceGrid(mat, ind110);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = -B+A;//C-B;
        n2t = B-C;//B;
        n3t = -D+A;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";        
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }   
     // (1,1,1), (1,0,1), (0,0,1), (1,1,0)
     A = interfaceGrid(mat, ind111); B = interfaceGrid(mat, ind101); 
     C = interfaceGrid(mat, ind001); D = interfaceGrid(mat, ind110);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = B-C;
        n2t = -B+A;
        n3t = -D+A;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";        
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }   
     // (1,0,0), (1,0,1), (0,0,1), (1,1,0)
     A = interfaceGrid(mat, ind100); B = interfaceGrid(mat, ind101); 
     C = interfaceGrid(mat, ind001); D = interfaceGrid(mat, ind110);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = B-C;
        n2t = D-A;
        n3t = B-A;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";        
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }   
     // (0,1,1), (0,1,0), (1,1,0), (0,0,1)
     A = interfaceGrid(mat, ind011); B = interfaceGrid(mat, ind010); 
     C = interfaceGrid(mat, ind110); D = interfaceGrid(mat, ind001);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = C-B;
        n2t = -D+A; //B;
        n3t = -B+A; //D;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";        
        n1 += n1t/L;
        n2 += n2t/L;
        n3 += n3t/L;
        
        numVectors++;
     }   
     // (0,0,0), (1,0,0), (1,1,0), (0,0,1)
     A = interfaceGrid(mat, ind000); B = interfaceGrid(mat, ind100); 
     C = interfaceGrid(mat, ind110); D = interfaceGrid(mat, ind001);
     if (A*B <= 0 || B*C <= 0 || A * C <= 0 || A * D <= 0 || B * D <= 0 || C * D <= 0) {
        n1t = B-A;//C-B;
        n2t = C-B;//B;
        n3t = D-A;
        double L = sqrt(n1t*n1t + n2t*n2t + n3t*n3t);
        cerr << "n = (" << (n1t/L) << "," << (n2t/L) << "," << (n3t/L) << ")\n";        
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
        cerr << "No normal in this cell...\n";
        //for (int ii = i-2; ii <= i+2; ii++) {
        //    for (int jj = j-2; jj <= j+2; jj++) {
        //        for (int kk = k-2; kk <= k+2; kk++) {
        //            cerr << interfaceGrid[ii][jj][kk] << " ";
        //        }
        //        cerr << endl;
        //    }
        //    cerr << endl;
        //}
     }
}
