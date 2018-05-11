#ifndef __NODE_H
#define __NODE_H

// Represents a grid cell.

class Node {
private:
	double U;  // arrival time
	int x, y, z;  // Space location

public:

	Node() : U(0.0), x(0), y(0), z(0) {}
	Node(double _U, int _x, int _y, int _z) :
		U(_U), x(_x), y(_y), z(_z) {}

	double getU() const { return U; }
	int X() const { return x; }
	int Y() const { return y; }
	int Z() const { return z; }

	void updateU(double newval) { U = newval; }

	// Arrives before or after?
	bool operator< (const Node & other) {
		return this->U < other.U;
	}

	bool operator> (const Node & other) {
		return this->U > other.U;
	}

	// code is losing it's elegance
	int heapIndex; // for back pointer purposes
};

#endif
