/*****************************
 * By:			Ethan Corgatelli
 * File:		nn.cpp
 * Project:		Assignment 1 (improved!)
 * Class:		cs475
 * Asn. Pg:		http://marvin.cs.uidaho.edu/Teaching/CS475/pas01.pdf
 *
 *****************************/
#include "nn.h"

// ---- constants ---- //
const double eta = 2.0;
const double slope = 0.1;
const double bias = -1.0;
const int numPasses = 10000;
const double spread = 0.5;

const bool debug = false;
const int peek = 5;

int main (int argc, char *argv[]) {

	srand(time(NULL));
	initRand();

	// ---- read & process the input ---- //

	// read in the number of features to be used for each data point
	int N; cin >> N;

	// initialize and name arrays
		// setup matrices
	Matrix trI = Matrix("in");    // training input
	Matrix trO = Matrix("out");   // training output
	
		// for the training
	Matrix w = Matrix("w"); // weights
	Matrix x = Matrix("x"); // in
	Matrix d = Matrix("d"); // delta
	Matrix y = Matrix("y"); // out
	Matrix t = Matrix("t"); // expected
	Matrix df = Matrix("df"); // difference (t - y)

		// for the output
	Matrix r = Matrix("r"); // recursivly built row by row
	Matrix m = Matrix("m"); // final output matrix


	// read in the input data
	trI.read();

	// extract the expected outputs
	trO = trI.extract(0, N, 0, 0);

	// truncate the data (cut out expected outputs)
	trI = Matrix(trI.extract(0, 0, 0, N));

	// normalize the data
	/*Matrix norm = */trI.normalizeCols();

	// add bias column vector
	Matrix bV = Matrix(trI.numRows(), 1, bias, "bias");
	trI = trI.joinRight(bV);
	
	// settup weights vector
	w = Matrix(trI.numCols(), trO.numCols(), 0.0)
							.rand(-spread, spread);
	
	// settup delta vector
	d = Matrix(trI.numCols(), 1, 0.0);

	// --- print input that was read ---- //
	if (!debug) {
		trI.print("Inputs + bias:");
		trO.print("Outputs:");
	}


	

	// ---- main training loop ---- //
	for (int itr = 0; itr < numPasses; itr++) {
		x = Matrix(trI);
		t = Matrix(trO);

		y = x.dot(w).map(&transfer);
		df = t.sub(y);
		
		d = x.Tdot(df).scalarMul(eta);
		w = w.add(d);

		if (debug && (itr < peek || itr > numPasses - peek)) {
			printf("\nitr: %d\n", itr);
			//w.print(); x.print(); t.print(); y.print(); d.print();
			printf("diff:\n");
			df.print();
		}
	}

	// ---- get results based on test data ---- //
	for (int i = 0; i < tsI.numRows(); i++) {
		x = tsI.extract(i, 0, 1, 0);

		y = x.dot(w).map(&transfer);
		df = t.sub(y);

		y = y.joinRight(df);
		r = r.joinBottom(y);
	}

	//m = tsIU.joinRight(r);
	//m.print("Tests + bias + results + diff:");
}

double transfer (double v) {
	return 1.0 / (1.0 + exp(-slope * v)); //sigmoid function
	//return *v < 0.5 ? 0.0 : 1.0; //square cutoff
}
