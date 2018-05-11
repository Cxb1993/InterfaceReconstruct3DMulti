#pragma once
class gradientElement
{
	int i, j, k;
	double dist;
	double * volumeContribution;
public:

	gradientElement(void);
	gradientElement(double d): dist(d), volumeContribution(nullptr) { clearContributions(); }
	~gradientElement(void);

	void normalize();
	void clearContributions();

	// Simple for now -- in the future may not store all 27 (since most will have only 1 defined, and at most 8 are possible)
	// but if it's never a memory or otherwise bottleneck, then it will remain this way
	double getContribution(int whichCell) { 
		if (volumeContribution == nullptr) { 
			return 0.0; 
		} else return volumeContribution[whichCell]; 
	}

	void   addContribution(int whichCell, double val) { 
		if (volumeContribution == nullptr) {
			clearContributions();
		}
		volumeContribution[whichCell] += val; 
	}

	int X() { return i; }
	int Y() { return j; }
	int Z() { return k; }
	double U() { return dist; }
	void setU(double u) { dist = u; }

	void setXYZ(int x, int y, int z) { i = x; j = y; k = z; }
};

