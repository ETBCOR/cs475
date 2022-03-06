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
double NeuralNet::eta = 2.0;
static double slope = 0.1;
double NeuralNet::bias = -1.0;
double NeuralNet::spread = 0.5;

const int numPasses = 1;

bool NeuralNet::debug = false;

int main (int argc, char *argv[]) {

	// check that args are good
	if (argc-1 > 1) {
		printf("%d is not a valid number of inputs. Please provide a single integer (0-9) for number of nodes in the hidden layer.\n", argc-1);
		return -1;
	}

	int nodes = 1;
	if (argc-1 == 1) {
		nodes = atoi(argv[1]);
		
		if (nodes < 1 || nodes > 10) {
			printf("%d is not a valid number of nodes. Please provide a single integer (1-10) for the number of nodes in the hidden layers.\n", nodes);
			return -1;
		}
	}

	srand(time(NULL));
	initRand();

	NeuralNet nn = NeuralNet(nodes);
	nn.read();
	nn.train(numPasses);
	//nn.guess();
}

// constructor
NeuralNet::NeuralNet(int _nodes) {

	nodes = _nodes;

	trI = new Matrix("in");  // training input
	trO = new Matrix("out"); // training output

	X = new Matrix("X");
	V = new Matrix("V");
	H = new Matrix("H");
	W = new Matrix("W");
	Y = new Matrix("Y");
	T = new Matrix("T");
	
	r = new Matrix("r"); // recursivly built row by row
	m = new Matrix("m"); // final output matrix
}

// ---- read & process the input ---- //
void NeuralNet::read() {

	// read in the number of features to be used for each data point
	cin >> N;

	// read in the input data
	trI->read();

	// extract the expected outputs
	*trO = trI->extract(0, 0, 0, N+1);
	M = trO->numCols();
	//trO->print();

	// truncate the data (cut out expected outputs)
	*trI = Matrix(trI->extract(0, 0, 0, N));

	// normalize the data
	trI->normalize();

	*V = Matrix(N+1, nodes, 0.0).rand(-spread, spread);
	*W = Matrix(nodes+1, M, 0.0).rand(-spread, spread);

	//trI->print();
	trO->print("Target:");

}

void NeuralNet::train(int passes) {

	for (int i = 0; i < passes; i++) {
		// Add the bias to the inputs X by adding a -1 column. X -> Xplus.
		Matrix bv = Matrix(trI->numRows(), 1, bias, "bias");
		*X = trI->joinRight(bv);

		// Retrieve expected matrix
		*T = Matrix(*trO);

		// H = f(Xplus V)
		*H = X->dot(*V).map(&transfer);

		// Add the bias to the hidden layer by adding a -1 column. H -> Hplus.
		bv = Matrix(H->numRows(), 1, bias, "bias");
		*H = H->joinRight(bv);

		// Y = f(Hplus W)
		*Y = H->dot(*W).map(&transfer);

		
		// Now the backward propagation phase

		// Wdelta = (Y - T) * Y * (1 - Y)
		Matrix diff = Y->sub(*T);
		Matrix omY = Matrix(Y->numRows(), Y->numCols(), 1.0, "omY");
		omY = omY.sub(*Y);
		Matrix Wd = diff.mul(*Y).mul(omY);

		// Hdelta = Hplus * (1 - Hplus) * (Wdelta Wtrans)
		Matrix omH = Matrix(H->numRows(), H->numCols(), 1.0, "omH");
		omH = omH.sub(*H);
		Matrix WdW = Wd.dotT(*W);
		Matrix Hd = H->mul(omH).mul(WdW);

		// W -= eta(Hplus-trans Wdelta) 
		*W = W->sub(H->Tdot(Wd).scalarMul(eta));

		// Hdelta -> Hdelta-mnus
		Hd = Hd.extract(0, 0, Hd.numCols() - 1, 0);

		// V -= eta(Xplus-trans Hdelta-minus)
		Matrix d = X->Tdot(Hd).scalarMul(eta);
		*V = V->sub(d);
		
		
	}
}

void NeuralNet::guess() {

	// Add the bias to the inputs X by adding a -1 column. X -> Xplus.
	Matrix bv = Matrix(trI->numRows(), 1, bias, "bias");
	*X = trI->joinRight(bv);

	// Retrieve expected matrix
	*T = trO;

	// H = f(Xplus V)
	*H = X->dot(V).map(&transfer);

	// Add the bias to the hidden layer by adding a -1 column. H -> Hplus.
	bv = Matrix(H->numRows(), 1, bias, "bias");
	*H = H->joinRight(bv);

	// Y = f(Hplus W)
	*Y = H->dot(W).map(&transfer);

	//trO->print("Target:");
	Y->print("Predicted:");		

}

double transfer (double v) {
	return 1.0 / (1.0 + exp(-slope * v)); //sigmoid function
	//return *v < 0.5 ? 0.0 : 1.0; //square cutoff
}
