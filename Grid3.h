#pragma once

// Three dimensional grid that indexes from [0..M][0..N][0..P] (Yes inclusive.)
template <typename T>
class Grid3
{
	T * data;

	int M, N, P;
public:
	const int left, right, top, bottom, forward, backward;

	Grid3(int m, int n, int p): M(m), N(n), P(p), forward(1), backward(-1), top(p+1), bottom(-p-1), right((n+1)*(p+1)), left(-(n+1)*(p+1))
	{
		data = new T[(m+1)*(n+1)*(p+1)];
	}
	
	~Grid3(void)
	{
		delete[] data;
	}

	T & operator()(int i, int j, int k) {
		return data[k + (P+1) * (j + (N+1) * i)];
	}

	T & operator()(int index) { return data[index]; }
};

