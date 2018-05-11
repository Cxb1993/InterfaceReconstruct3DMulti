/* Represents the heap used for this problem. */
#include "Node3D.h"
#include "Grid3.h"

class EikonalHeap3D {
private:
    Node3D *heap;
    int capacity, size;
	
    // Back pointers.
	Grid3<int> indices;
    int nx, ny, nz; // Internal rep of grid for backpointer access
    
    // Updates back pointers too!
    void swapElements(int i, int j) {
         Node3D jth = heap[i];
         heap[i] = heap[j];
         heap[j] = jth;
         
         Node3D ith = heap[i];
		 indices(ith.Index()) = i;
         indices(jth.Index()) = j;
    }
    
    void heapify(int index);

public:
    EikonalHeap3D(int initialCapacity, int _nx, int _ny, int _nz):
	   indices(_nx, _ny, _nz) {
       heap = new Node3D[initialCapacity + 1];
       size = 0;
       capacity = initialCapacity;
	   for (int inds = 0; inds < (_nx+1) * (_ny+1) * (_nz+1); ++inds) {
		   indices(inds) = 0;
	   }
    }

    ~EikonalHeap3D() { 
       delete[] heap; 
    }
    
    Grid3<int> & getBackpointerArray() { return indices; }

    void insert(Node3D next);
    Node3D remove();

    bool isEmpty() { return size == 0; }
    
    double getElementTime(int index) {
         return heap[index].getU();
    }
    
    void updateElementTime(int index, double newtime, double newvelocity) {
         heap[index].update(newtime, newvelocity);
         heapify(index);
    }

	void discardAll();

	int currentSize() { return size; }
};
