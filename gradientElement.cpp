// Stub
#include "stdafx.h"
#include "gradientElement.h"
#include <math.h>


gradientElement::gradientElement(void) : dist(-1.0)
{
	volumeContribution = nullptr;
}


gradientElement::~gradientElement(void)
{
	if (volumeContribution != nullptr) {
		delete[] volumeContribution;
	}
	volumeContribution = nullptr;
}

void gradientElement::normalize() {
	double L = 0.0;
	for (int i = 0; i < 27; i++) {
		L += fabs(volumeContribution[i]);
	}

	// Ignores larger issue: why is zero delta volume being computed for a thing supposedly close to an interface
	if (L < 1e-8) { return; }

	for (int i = 0; i < 27; i++) {
		volumeContribution[i] /= L;
	}
}

void gradientElement::clearContributions() {
	if (volumeContribution == nullptr) {
		volumeContribution = new double[27];
	}
	for (int c = 0; c < 27; c++) { volumeContribution[c] = 0.0; }
}