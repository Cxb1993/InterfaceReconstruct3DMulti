#pragma once

template <typename T>
class MultiGrid3
{
	T ** grids;
	int M, N, P;
public:
	const int numGrids, dataSize;
	const int left, right, top, bottom, forward, backward;

	MultiGrid3(int nG, int m, int n, int p): numGrids(nG), dataSize((m+1)*(n+1)*(p+1)), forward(1), backward(-1), top(p+1), 
		                                    bottom(-p-1), right((n+1)*(p+1)), left(-(n+1)*(p+1)), M(m), N(n), P(p) {
	   grids = new T*[nG];
	   for (int G = 0; G < nG; G++) {
		   grids[G] = new T[dataSize];
	   }
	}

	~MultiGrid3(void) {
		for (int G = 0; G < numGrids; G++) {
		    delete[] grids[G];
		}
		delete[] grids;
	}

	T & operator()(int mat, int i, int j, int k) {
		return grids[mat][k + (P+1) * (j + (N+1) * i)];
	}

	T & operator()(int mat, int index) { return grids[mat][index]; }
};

