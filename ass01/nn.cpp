/*****************************
 * File:	warmup.cpp
 * Project:	Assignment 1
 * Class:	cs475
 * Asn. Pg:	http://marvin.cs.uidaho.edu/Teaching/CS475/pas01.pdf
 *
 *****************************/
#include <iostream>
#include "mat.h"
#include "rand.h"

// ---- constants ---- //
const double eta = 0.001;
const double slope = 0.001;
const double spread = 0.0;
const double bias = 1.0;
const int numPasses = 100;

const bool debug = false;

int main (int argc, char *argv[]) {

	initRand();

	// ---- read & process the input ---- //

	// read in the number of features to be used for each data point
	int N; cin >> N;

	// initialize and name arrays
	Matrix trI = Matrix("in");
	Matrix trO = Matrix("out");
	Matrix tsI = Matrix("in");
	Matrix tsO = Matrix("out");

	// read in the input data
	trI.read();
	tsI.read();
	// extract the expected outputs
	trO = trI.extract(0, N, 0, 0);
	tsO = tsI.extract(0, N, 0, 0);
	// truncate the data (cut out expected outputs)
	trI = Matrix(trI.extract(0, 0, 0, N));
	tsI = Matrix(tsI.extract(0, 0, 0, N));
	// normalize the data
	trI.normalize();
	tsI.normalize();
	// add bias column vector
	Matrix bV = Matrix(trI.numRows(), 1, 1.0);
	trI = trI.joinRight(bV);
		   bV = Matrix(tsI.numRows(), 1, 1.0);
	tsI = tsI.joinRight(bV);

	// --- print input that was read ---- //
	trI.print("Inputs + bias:");
	trO.print("Outputs:");
	tsI.print("Test Inputs + bias:");
	tsO.print("Expected Outputs:");


	// shuffle the training data
	trI.shuffle();
	
	// ---- main training loop ---- //
	for (int i = 0; i < numPasses; i++) {
		// do the stuff

	}	


	//m.print("Tests + bias + results + diff:");
}
