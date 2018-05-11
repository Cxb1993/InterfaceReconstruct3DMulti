#include "stdafx.h"
#include "Heap3D.h"
/* Implementation of Heap. */

// Add an element to a min heap.
// Indexing starts from 1, FYI.
void EikonalHeap3D::insert(Node3D next) {

	int ind = next.Index();
	double U = next.getU(), V = next.getV();

	int nodeIndex = indices(ind);
	if (nodeIndex != 0) {
		if (U < getElementTime(nodeIndex)) {
			updateElementTime(nodeIndex, U, V);
		}
		return;
	}

	if (size == capacity) {
		capacity *= 2;
		Node3D * tmp = new Node3D[capacity + 1];
		for (int i = 1; i <= size; i++) {
			tmp[i] = heap[i];
		}
		delete[] heap;
		heap = tmp;
	}

	// Place new item at the bottom of the heap.
	// We'll see how this performs: Making an explicit copy.
	heap[++size] = next;
	indices(ind) = size;

	heapify(size);
}

// Remove the largest item.
Node3D EikonalHeap3D::remove() {
	// Grab necessary item - copy it
	Node3D returnVal = heap[1];
	// Put bottom heap item at root
	swapElements(1, size);
	--size;
	// Clean back pointer to removed node
	indices(returnVal.Index()) = 0;
	int i = 1;
	// Bubble the bigger of the children to this spot.
	// Repeat until the originally swapped item is bigger than its children.
	while (2 * i <= size) {
		int left = 2 * i, right = left + 1;
		// All these conditions are for things like "Right child does not exist"
		// and etc.
		// Idea: Making it so heap[size+1] val is always really big will eliminate one if statement... if that helps at all...
		if (right <= size) {
			if (heap[left] < heap[i] && heap[left] < heap[right]) {
				swapElements(i, left);
				i = left;
			}
			else if (heap[right] < heap[i]) {
				swapElements(i, right);
				i = right;
			}
			else {
				break;
			}
		}
		else if (heap[left] < heap[i]) {
			swapElements(i, left);
			i = left;
		}
		else {
			break;
		}
	}

	return returnVal;
}

void EikonalHeap3D::heapify(int i) {
	while (i >= 2 && heap[i / 2] > heap[i]) {
		swapElements(i, i / 2); // Updates back pointers automatically.
		i = i / 2;
	}
}

void EikonalHeap3D::discardAll() {
	for (int i = 1; i <= size; i++) {
		Node3D node = heap[i];
		indices(node.Index()) = 0;
	}
	size = 0;
}
