/*****************************
 * File:	warmup.cpp
 * Project:	Assignment 1
 * Class:	cs475
 * Asn. Pg:	http://marvin.cs.uidaho.edu/Teaching/CS475/pas01.pdf
 *
 *****************************/
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "mat.h"
#include "rand.h"

// ---- constants ---- //
const double eta = 0.0001;
const double slope = 0.0001;
const double bias = 1.0;
const int numPasses = 1000;
//const double spread = 0.0;

const bool debug = false;

double transfer (int, double*);

int main (int argc, char *argv[]) {

	srand(time(NULL));
	initRand();

	// ---- read & process the input ---- //

	// read in the number of features to be used for each data point
	int N; cin >> N;

	// initialize and name arrays
	Matrix trI = Matrix("in");  // training input
	Matrix trO = Matrix("out"); // training output
	Matrix tsI = Matrix("in");  // test input
	Matrix tsO = Matrix("out"); // test output
	Matrix w = Matrix("w"); // weights
	Matrix d = Matrix("d"); // delta

	Matrix x = Matrix("x"); // 
	Matrix y = Matrix("y"); // 
	Matrix t = Matrix("t"); // 
	Matrix df = Matrix("df"); // diff (t - y)

	// read in the input data
	trI.read();
	tsI.read();
	// extract the expected outputs
	trO = trI.extract(0, N, 0, 0);
	tsO = tsI.extract(0, N, 0, 0);
	// normalize the data
	trI.normalize();
	tsI.normalize();
	// truncate the data (cut out expected outputs)
	trI = Matrix(trI.extract(0, 0, 0, N));
	tsI = Matrix(tsI.extract(0, 0, 0, N));
	// add bias column vector
	Matrix bV = Matrix(trI.numRows(), 1, bias);
	trI = trI.joinRight(bV);
		   bV = Matrix(tsI.numRows(), 1, bias);
	tsI = tsI.joinRight(bV);

	// settup weights vector
	w = Matrix(trI.numCols(), trO.numCols(), 0.0).rand(-1.0, 1.0);
	
	// settup delta vector
	d = Matrix(trI.numCols(), 1, 0.0);

	// --- print input that was read ---- //
	/*trI.print("Inputs + bias:");
	trO.print("Outputs:");
	tsI.print("Test Inputs + bias:");
	tsO.print("Expected Outputs:");*/

	// shuffle the training data
	//trI.shuffle();

	// ---- main training loop ---- //
	for (int itr = 0; itr < numPasses; itr++) {	
		int i = rand() % trI.numRows();

		x = trI.extract(i, 0, 1, 0);
		t = trO.extract(i, 0, 1, 0);
		y = x.dot(w);//.mapCol(&transfer);
		df = t.sub(y); //df.print();
		
		d = x.Tdot(df).scalarMul(eta);
		w = w.add(d);

		if (/*debug*/ itr < 10 || itr > numPasses-10) {
			printf("\nitr: %d\n", itr);
			w.print();
			x.print();
			t.print();
			y.print();
			d.print();
			printf("diff:\n");
			df.print();
		}	
	}
	//m.print("Tests + bias + results + diff:");
}

double transfer (int s, double *v) {
	return 1.0 / (1.0 + exp(-slope * *v));
}
