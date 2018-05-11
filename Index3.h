#pragma once
class Index3
{
	int m, n, p;
	int imin, jmin, kmin, imax, jmax, kmax;
public:
	int i, j, k, ind;

	Index3(int il, int jl, int kl, int ir, int jr, int kr, int M, int N, int P) :
		imin(il), jmin(jl), kmin(kl), imax(ir), jmax(jr), kmax(kr), m(M), n(N), p(P)
	{
		i = imin; j = jmin; k = kmin;
		ind = kmin + (p + 1) * (jmin + (n + 1) * imin);
	}

	Index3(int M, int N, int P) :
		imin(0), jmin(0), kmin(0), imax(M), jmax(N), kmax(P), m(M), n(N), p(P)
	{
		i = imin; j = jmin; k = kmin;
		ind = kmin + (p + 1) * (jmin + (n + 1) * imin);
	}

	Index3 & operator++() {
		if (++k > kmax) {
			k = kmin;
			if (++j > jmax) {
				j = jmin;
				++i;
			}
			ind = k + (p + 1)*(j + (n + 1)*i);
		}
		else {
			ind++;
		}
		return *this;
	}

	bool valid() { return i <= imax && j <= jmax && k <= kmax; }

	operator int() { return ind; }

	~Index3(void)
	{
	}
};

