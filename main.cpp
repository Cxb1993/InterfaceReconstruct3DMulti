#include "stdafx.h"
#include <cstdlib>
#include <sstream>
#include <string>
#include <ctime>
#include "Interface3D.h"

#include "MultiGrid.h"

using namespace std;

int main(int argc, char *argv[])
{
	//    if (argc < 3) {
	//        cerr << "Specify input file name and number of subdivisions.\n";
	//        return 1;
	//    }

	/*    ifstream infile(argv[1]);
		int gridM, gridN, gridP;
		double H;

		infile >> gridM >> gridN >> gridP >> H;

		int ns = atoi(argv[2]);
		double h = H/(double)ns;

		cerr << "Ok...\n";

		Interface3D Ivf3(gridM, gridN, gridP, ns, h, infile);*/
	int ns = 3;

	if (argc > 1) {
		ns = atoi(argv[1]);
	}

	Interface3D Ivf3(4, 24, 24, 24, ns);
	//Ivf3.initTripleLine(0.11, 0.0, 0.0, 1.0, 0.5, 0.0);
	//Ivf3.initTripleLineFlat(0.11, 0.0, 0.0, 1.0, 0.5, 0.0);
	Ivf3.initQuadruplePoint();
	//Ivf3.initRandomVoronoi(11);
	//Ivf3.initHemispheres(2.0, 0.0, 0.0, 0.0, M_PI / 5.0, M_PI / 5.0, 0.0, 6.0);
	//Ivf3.initCubeSphereIntersect();
	//Ivf3.initConcentricCircles();
	//Ivf3.initCube(3.0, 0.0, 0.0, 0.0, M_PI / 5.0, M_PI / 5.0, 0.0, 6.0);
	//Ivf3.initSphere(2.0, 0.0, 0.0, 0.0);
	//Ivf3.initDeltoidish(3.0, 0.0, 0.0, 0.0, M_PI / 5.0, M_PI / 5.0, 0.0, 6.0);
	//Ivf3.initCone(1.5, 4.0, 0.0, 0.0, 0.0, M_PI / 5.0, M_PI / 5.0);
	//Ivf3.initTetra();
	//Interface3D * theBay = Interface3D::SFBay(ns);
	//Interface3D & Ivf3 = *theBay;
	double h = Ivf3.hz;
	//Ivf3.advancedSeeding();

	Ivf3.seedGrid(0.5);

	Ivf3.outputNarrowBand("HelloWorld.vti");

	Ivf3.reinitialize();

	Ivf3.outputInterface("seed.vti");

	// seed with sphere, then call reinitialize.    

	srand((int)time(NULL));
	// epsilon = 0.04 = good
	//double epsilon = 0.04, timeStep = h*h/(4*epsilon)/1.3; 
	double epsilon = 0.04, timeStep = h * h / (4 * epsilon) / 3.0;

	cerr << "Initialization complete\n";

	time_t tstart = clock();

	const int stageOneIters = 20, stageTwoIters = 15;

	for (int iter = 0; iter < stageOneIters + stageTwoIters; iter++) {
		Ivf3.saveInterface();

		if (iter == stageOneIters / 2 || iter == (stageOneIters + stageTwoIters / 2)) {
			Ivf3.reinitialize();
		}
		else if (iter == stageOneIters) {
			Ivf3.checkResult(false);
			Ivf3.reinitialize();
			Ivf3.outputInterface("midterm.vti");
			Ivf3.outputInterfaceMesh("MidtermMesh.vtp");
		}
		else if (iter == stageOneIters + stageTwoIters - 3) {
			//timeStep /= 3.0;
		}

		//if (iter % 5 == 0 || iter > 30) {
		//    stringstream filename;
		//    filename << "result" << iter << ".vti";
		//    Ivf3.outputInterfaceCoarse(filename.str().c_str());
		//  
		//    //Ivf3.checkResult(true);
		//    int tend = clock();
		//    cerr << "Clock ticks = " << (tend - tstart) << endl;
		//} 

		//Ivf3.checkResult(false);

		//if (iter == 200) {
		//    Ivf3.outputInterface("midterm.vti");
		//}
		if (iter < stageOneIters) {
			Ivf3.iterateLevelSetMethod(epsilon, timeStep);
		}
		else {
			Ivf3.iterateLevelSetMethod(epsilon, timeStep, true); // with volume projections
		}
		time_t tend = clock() - tstart;
		cerr << tend << ": Completed step " << iter << " with change " << Ivf3.interfaceChange() << endl;
		//Ivf3.checkResult(false);

		/*stringstream resultFile;
		resultFile << "result" << iter << ".vti";
		Ivf3.outputInterface(resultFile.str().c_str());*/
	}

	// Hmm...
	//Ivf3.computeCurvatureFlow(epsilon);
	//Ivf3.updateInterface(timeStep);
	//Ivf3.outputInterface("CurvCheck.vti");
	//Ivf3.computeVolumeFlow(0.0, timeStep, false);
	//Ivf3.updateInterface(timeStep);

	cerr << "Checking result: " << endl;

	Ivf3.outputNarrowBand("NarrowBandFinal.vti");
	Ivf3.outputInterface("final.vti");

	Ivf3.outputInterfaceCoarse("coarse.vti");

	time_t tend = clock();
	cerr << "Clock ticks = " << (tend - tstart) << endl;

	Ivf3.outputInterfaceMesh("FinalMesh.vtp");
	system("pause");
	Ivf3.checkResult(false, true);

	// Let's do a bunch of random test cases.
/*    for (int i = 0; i < 50; i++) {
		double a = rand() / 32767.0 - 0.5;
		double b = rand() / 32767.0 - 0.5;
		double c = rand() / 32767.0 - 0.5;
		double d = rand() / 32767.0 - 0.5;

		int numMinus = 0;
		if (a <= 0) numMinus++;
		if (b <= 0) numMinus++;
		if (c <= 0) numMinus++;
		if (d <= 0) numMinus++;

		cerr << "Comparison: Exact vs monte carlo: " << Ivf3.fluidTetrahedron(a,b,c,d) << " " <<
														Ivf3.volApprox(a,b,c,d) << "; case " << numMinus << endl;
	}*/

	//if (strcmp(argv[1], "circ3DR2.txt") == 0) {
	Ivf3.checkSphereR2Normal();
	//}

	system("PAUSE");
	return EXIT_SUCCESS;
}
