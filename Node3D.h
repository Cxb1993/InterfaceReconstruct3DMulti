#ifndef __NODE_H
#define __NODE_H

// Represents a grid cell.

class Node3D {
private:
	double U, V;  // arrival time
	int index;  // Space location

	// temporary: grid cell of corresponding interface
	//int gridcell;

public:

	Node3D() : U(0.0), index(0) {}
	Node3D(double _U, double _V, int _index) :
		U(_U), V(_V), index(_index) {}

	double getU() const { return U; }
	double getV() const { return V; }

	int Index() const { return index; }

	void setV(double v) { V = v; }

	void update(double newU, double newV) { U = newU; V = newV; }

	// Arrives before or after?
	bool operator< (const Node3D & other) {
		return this->U < other.U;
	}

	bool operator> (const Node3D & other) {
		return this->U > other.U;
	}

	int heapIndex; // for back pointer purposes

	// Temporary
	//void setcell(int c) { gridcell = c; }
	//Node3D(double _U, double _V, int _x, int _y, int _z, int cell):
	//            U(_U), V(_V), x(_x), y(_y), z(_z), gridcell(cell) {}
	//int getcell() { return gridcell; }
};

#endif
