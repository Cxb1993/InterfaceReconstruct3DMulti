#include "stdafx.h"
#include "Interface3D.h"
#include <sstream>

using namespace std;

char base64table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

void Interface3D::outputInterface(const char * filename) {
	ofstream outfile(filename);

	outfile << "<VTKFile type=\"ImageData\">" << endl;
	outfile << " <ImageData WholeExtent=\"0 " << p << " 0 " << m << " 0 " << n << "\"" << endl;
	outfile << "  Origin=\"0 0 0\" Spacing=\"1 1 1\">" << endl;
	outfile << "   <Piece Extent=\"0 " << p << " 0 " << m << " 0 " << n << "\">" << endl;
	outfile << "    <PointData>" << endl;

	for (Material = 0; Material < numPhases; Material++) {
		outfile << "     <DataArray type=\"Float32\" Name=\"Phi" << Material << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		//char base64[ceil((m+1)*(n+1)*(p+1)*sizeof(float)*4.0/3.0)];
		//int nleftoverbits=0, leftvogmaerval, current = 0;
		int inds = 0;
		for (int i = 0; i <= m; i++) {
			for (int j = 0; j <= n; j++) {
				for (int k = 0; k <= p; k++) {
					outfile << interfaceGrid(Material, inds) << " ";
					inds++;
				}
				outfile << endl;
			}
		}
		outfile << "    </DataArray>" << endl;
	}
	/*outfile << "     <DataArray type=\"Float32\" Name=\"Error\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				int I = i/ns, J = j/ns, K = k/ns;
				double err = (grid[I+1][J+1][K+1] - vActs[I+1][J+1][K+1]);
				if (i % ns == 0 || j % ns == 0 || k % ns == 0) {
					err = 0.0;
				}
				outfile << err << " ";
			}
			outfile << endl;
		}
	}
	outfile << "    </DataArray>" << endl;*/


	/*outfile << "     <DataArray type=\"Float32\" Name=\"Curvature\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= p; k++) {
				double right, left, up, down, forward, back;
				double topright, topleft, botright, botleft, center;
				double forwardright, forwardleft, backright, backleft;
				double forwardtop, forwardbot, backtop, backbot;

				// periodic... why is it periodic again?
				center = interfaceGrid[i][j][k];
				right = interfaceGrid[(i+1)%m][j][k];
				left = interfaceGrid[(i+m-1)%m][j][k];
				up = interfaceGrid[i][(j+1)%n][k];
				down = interfaceGrid[i][(j+n-1)%n][k];
				forward = interfaceGrid[i][j][(k+1)%p];
				back = interfaceGrid[i][j][(k+p-1)%p];

				topright = interfaceGrid[(i+1)%m][(j+1)%n][k];
				topleft = interfaceGrid[(i+m-1)%m][(j+1)%n][k];
				botright = interfaceGrid[(i+1)%m][(j+n-1)%n][k];
				botleft = interfaceGrid[(i+m-1)%m][(j+n-1)%n][k];

				forwardright = interfaceGrid[(i+1)%m][j][(k+1)%p];
				forwardleft = interfaceGrid[(i+m-1)%m][j][(k+1)%p];
				backright = interfaceGrid[(i+1)%m][j][(k+p-1)%p];
				backleft = interfaceGrid[(i+m-1)%m][j][(k+p-1)%p];

				forwardtop = interfaceGrid[i][(j+1)%n][(k+1)%p];
				forwardbot = interfaceGrid[i][(j+n-1)%n][(k+1)%p];
				backtop = interfaceGrid[i][(j+1)%n][(k+p-1)%p];
				backbot = interfaceGrid[i][(j+n-1)%n][(k+p-1)%p];

				double uxx = (right - 2 * center + left)/(hx*hx),
					   uyy = (up - 2 * center + down)/(hy*hy),
					   uzz = (forward - 2 * center + back)/(hz*hz);
				double uxy = ((topright - topleft)/(2*hx) - (botright - botleft)/(2*hx)) / (2*hy),
					   uxz = ((forwardright - forwardleft)/(2*hx) - (backright - backleft)/(2*hx)) / (2*hz),
					   uyz = ((forwardtop - forwardbot)/(2*hy) - (backtop - backbot)/(2*hy)) / (2*hz);
				double ux  = (right - left)/(2*hx),
					   uy  = (up - down)/(2*hy),
					   uz  = (forward - back)/(2*hz);

				double gradLen = ux*ux + uy*uy + uz*uz;
				double curvature = ( (uyy + uzz) * ux * ux + (uxx + uzz) * uy * uy + (uxx + uyy) * uz * uz
							  - 2 * ux * uy * uxy - 2 * ux * uz * uxz - 2 * uy * uz * uyz)
							  / (gradLen * sqrt(gradLen) + hx*hz);

				outfile << curvature << " ";
			}
			outfile << endl;
		}
	}
	outfile << "    </DataArray>" << endl;    */


	outfile << "   </PointData>" << endl;
	outfile << "   <CellData> </CellData>" << endl;
	outfile << "  </Piece>" << endl;
	outfile << " </ImageData>" << endl;
	outfile << "</VTKFile>" << endl;

	outfile.close();
}

void Interface3D::outputInterfaceCoarse(const char * filename) {
	int gridM = m / ns, gridN = n / ns, gridP = p / ns;

	ofstream outfile(filename);

	outfile << "<VTKFile type=\"ImageData\">" << endl;
	outfile << " <ImageData WholeExtent=\"0 " << (gridP + 1) << " 0 " << (gridM + 1) << " 0 " << (gridN + 1) << "\"" << endl;
	outfile << "  Origin=\"0 0 0\" Spacing=\"1 1 1\">" << endl;
	outfile << "   <Piece Extent=\"0 " << (gridP + 1) << " 0 " << (gridM + 1) << " 0 " << (gridN + 1) << "\">" << endl;
	outfile << "    <PointData>" << endl;

	for (Material = 0; Material < numPhases; Material++) {
		outfile << "     <DataArray type=\"Float32\" Name=\"Phi" << Material << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;

		for (int I = 0; I < gridM + 2; I++) {
			for (int J = 0; J < gridN + 2; J++) {
				for (int K = 0; K < gridP + 2; K++) {
					if (I == 0 || I == gridM + 1 || J == 0 || J == gridN + 1 || K == 0 || K == gridP + 1) {
						outfile << 50.0 << " ";
						continue;
					}
					double sum = 0.0;
					for (int i = 0; i < ns; i++) {
						for (int j = 0; j < ns; j++) {
							for (int k = 0; k < ns; k++) {
								sum += interfaceGrid(Material, (I - 1)*ns + i, (J - 1)*ns + j, (K - 1)*ns + k);
							}
						}
					}
					sum /= ns * ns * ns;
					outfile << sum << " ";
				}
				outfile << endl;
			}
		}
		outfile << "    </DataArray>" << endl;

		outfile << "     <DataArray type=\"Float32\" Name=\"Error" << Material << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;

		int coarseInd = 0;
		for (int I = 0; I < gridM + 2; I++) {
			for (int J = 0; J < gridN + 2; J++) {
				for (int K = 0; K < gridP + 2; K++) {
					outfile << (grid(Material, coarseInd) - vActs(Material, coarseInd)) << " ";
					coarseInd++;
				}
				outfile << endl;
			}
		}
		outfile << "    </DataArray>" << endl;
	}

	outfile << "   </PointData>" << endl;
	outfile << "   <CellData> </CellData>" << endl;
	outfile << "  </Piece>" << endl;
	outfile << " </ImageData>" << endl;
	outfile << "</VTKFile>" << endl;

	outfile.close();
}

void Interface3D::outputNarrowBand(const char * filename) {
	ofstream outfile(filename);

	outfile << "<VTKFile type=\"ImageData\">" << endl;
	outfile << " <ImageData WholeExtent=\"0 " << p << " 0 " << m << " 0 " << n << "\"" << endl;
	outfile << "  Origin=\"0 0 0\" Spacing=\"1 1 1\">" << endl;
	outfile << "   <Piece Extent=\"0 " << p << " 0 " << m << " 0 " << n << "\">" << endl;
	outfile << "    <PointData>" << endl;

	for (int mat = 0; mat < numPhases; mat++) {
		outfile << "     <DataArray type=\"Float32\" Name=\"NB" << mat << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		//char base64[ceil((m+1)*(n+1)*(p+1)*sizeof(float)*4.0/3.0)];
		//int nleftoverbits=0, leftvogmaerval, current = 0;
		int ind = 0;
		for (int i = 0; i <= m; i++) {
			for (int j = 0; j <= n; j++) {
				for (int k = 0; k <= p; k++) {
					double val = narrowBand(mat, ind) ? -1.0 : 1.0;
					outfile << val << " ";
					ind++;
				}
				outfile << endl;
			}
		}
		outfile << "    </DataArray>" << endl;
	}

	outfile << "   </PointData>" << endl;
	outfile << "   <CellData> </CellData>" << endl;
	outfile << "  </Piece>" << endl;
	outfile << " </ImageData>" << endl;
	outfile << "</VTKFile>" << endl;

	outfile.close();
}

void Interface3D::outputVelocities(const char * filename) {
	ofstream outfile(filename);

	outfile << "<VTKFile type=\"ImageData\">" << endl;
	outfile << " <ImageData WholeExtent=\"0 " << p << " 0 " << m << " 0 " << n << "\"" << endl;
	outfile << "  Origin=\"0 0 0\" Spacing=\"1 1 1\">" << endl;
	outfile << "   <Piece Extent=\"0 " << p << " 0 " << m << " 0 " << n << "\">" << endl;
	outfile << "    <PointData>" << endl;

	outfile << "     <DataArray type=\"Float32\" Name=\"Phi\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	//char base64[ceil((m+1)*(n+1)*(p+1)*sizeof(float)*4.0/3.0)];
	//int nleftoverbits=0, leftvogmaerval, current = 0;
	for (Index3 inds(m, n, p); inds.valid(); ++inds) {
		outfile << interfaceGrid(Material, inds) << " ";
		if (inds % (p + 1) == 0) {
			outfile << endl;
		}
	}
	outfile << "    </DataArray>" << endl;

	outfile << "     <DataArray type=\"Float32\" Name=\"Psi\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	//char base64[ceil((m+1)*(n+1)*(p+1)*sizeof(float)*4.0/3.0)];
	//int nleftoverbits=0, leftvogmaerval, current = 0;
	for (Index3 inds(m, n, p); inds.valid(); ++inds) {
		outfile << velocities(inds) << " ";
		if (inds % (p + 1) == 0) {
			outfile << endl;
		}
	}
	outfile << "    </DataArray>" << endl;

	outfile << "   </PointData>" << endl;
	outfile << "   <CellData> </CellData>" << endl;
	outfile << "  </Piece>" << endl;
	outfile << " </ImageData>" << endl;
	outfile << "</VTKFile>" << endl;

	outfile.close();
}

void Interface3D::outputInterfaceMesh(const char * filename) {
	// Will be really slow for now. Can do lots of optimizations for not looking over the whole thing later.
	std::vector<double> points;
	std::vector<int> polySizes;
	std::vector<int> faceTypes;
	std::vector<double> errorVals;

	for (Index3 inds(0, 0, 0, m - 1, n - 1, p - 1, m, n, p); inds.valid(); ++inds) {
		int i = inds.i, j = inds.j, k = inds.k;
		int coarseInd = coarseIndexFromFine(i, j, k);
		// Do for all six tetrahedra
		// addPolygonsForTetrahedron(points, polySizes, ...)
		int indices[4][3], transform[3];
		// Each cube is divided into six tetrahedra
		// (1) (0,0,0), (0,1,0), (1,1,0), (0,0,1)
		fillIndices(indices, i, j, k, i, j + 1, k, i + 1, j + 1, k, i, j, k + 1);
		// This tetrahedron is not rotated or flipped. (Identity transformation)
		transform[0] = 1;  transform[1] = 2; transform[2] = 3;
		addPolygonsForTetrahedron(points, polySizes, faceTypes, errorVals, indices, transform, coarseInd);

		// (2) (0,1,1), (0,1,0), (1,1,0), (0,0,1)
		fillIndices(indices, i, j + 1, k + 1, i, j + 1, k, i + 1, j + 1, k, i, j, k + 1);
		transform[0] = 1;  transform[1] = -3; transform[2] = -2;
		addPolygonsForTetrahedron(points, polySizes, faceTypes, errorVals, indices, transform, coarseInd);

		// (3) (1,1,1), (0,1,1), (0,0,1), (1,1,0)
		fillIndices(indices, i + 1, j + 1, k + 1, i, j + 1, k + 1, i, j, k + 1, i + 1, j + 1, k);
		transform[0] = -2;  transform[1] = -1; transform[2] = -3;
		addPolygonsForTetrahedron(points, polySizes, faceTypes, errorVals, indices, transform, coarseInd);

		// (4) (1,1,1), (1,0,1), (0,0,1), (1,1,0)
		fillIndices(indices, i + 1, j + 1, k + 1, i + 1, j, k + 1, i, j, k + 1, i + 1, j + 1, k);
		transform[0] = -1;  transform[1] = -2;  transform[2] = -3;
		addPolygonsForTetrahedron(points, polySizes, faceTypes, errorVals, indices, transform, coarseInd);

		// (5) (1,0,0), (1,0,1), (0,0,1), (1,1,0)
		fillIndices(indices, i + 1, j, k, i + 1, j, k + 1, i, j, k + 1, i + 1, j + 1, k);
		transform[0] = -1;  transform[1] = 3;  transform[2] = 2;
		addPolygonsForTetrahedron(points, polySizes, faceTypes, errorVals, indices, transform, coarseInd);
		// (6) (0,0,0), (1,0,0), (1,1,0), (0,0,1)
		fillIndices(indices, i, j, k, i + 1, j, k, i + 1, j + 1, k, i, j, k + 1);
		transform[0] = 2;  transform[1] = 1;  transform[2] = 3;
		addPolygonsForTetrahedron(points, polySizes, faceTypes, errorVals, indices, transform, coarseInd);
	}

	ofstream outfile(filename);
	outfile << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">";
	outfile << " <PolyData>\n";
	outfile << "  <Piece NumberOfPoints=\"" << (points.size() / 3) << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\""
		<< " NumberOfStrips=\"0\" NumberOfPolys=\"" << polySizes.size() << "\">\n";
	outfile << "   <Points>\n";
	outfile << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

	for (std::vector<double>::iterator it = points.begin(); it != points.end(); it++) {
		double val = *it;
		outfile << val << " ";
	}
	outfile << endl;

	outfile << "    </DataArray>\n";
	outfile << "   </Points>\n";
	outfile << "   <CellData>\n";
	outfile << "    <DataArray type=\"Int32\" Name=\"cell_scalars\" format=\"ascii\">\n";
	for (std::vector<int>::iterator it = faceTypes.begin(); it != faceTypes.end(); it++) {
		outfile << (*it) << " ";
	}
	outfile << "\n    </DataArray>\n";
	outfile << "    <DataArray type=\"Float32\" Name=\"error\" format=\"ascii\">\n";
	for (std::vector<double>::iterator it = errorVals.begin(); it != errorVals.end(); it++) {
		outfile << (*it) << " ";
	}
	outfile << "\n    </DataArray>\n";
	outfile << "   </CellData>\n";
	outfile << "   <Polys>\n";
	outfile << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

	int start = 0;
	for (std::vector<int>::iterator it = polySizes.begin(); it != polySizes.end(); it++) {
		int val = *it;
		for (int c = 0; c < val; c++) {
			outfile << start << " ";
			start++;
		}
	}
	outfile << "\n";

	outfile << "    </DataArray>\n";
	outfile << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

	start = 0;
	for (std::vector<int>::iterator it = polySizes.begin(); it != polySizes.end(); it++) {
		int val = *it;
		start += val;
		outfile << start << " ";
	}
	outfile << "\n";

	outfile << "    </DataArray>\n";
	outfile << "   </Polys>\n";
	outfile << "  </Piece>\n";
	outfile << " </PolyData>\n";
	outfile << "</VTKFile>\n";
}