/*****************************
 * By:			Ethan Corgatelli
 * File:		pca.cpp
 * Project:		Assignment 3
 * Class:		cs475
 * Asn. Pg:		http://marvin.cs.uidaho.edu/Teaching/CS475/pas03.pdf
 *
 *****************************/
#include "pca.h"

int main (int argc, char *argv[]) {

	// variables
	int R, C, K;
	Matrix O, X, M, V, W;
	bool color, trans;

	// check that args are good
	if (argc-1 != 1) {
		printf("%d is not a valid number of inputs. Please provide a single integer for the number of eigenvectors to be kept.\n", argc-1);
		return -1;
	}

	// get K from command line args
	K = atoi(argv[1]);
	if (K == 0) {
		printf("cannot reduce to 0 dimensions\n");
		return -1;
	}

	// check if we should transpose before and after (denoted by negative arg)
	if (K < 0) {
		K = -K;
		trans = true;
	} else trans = false;

	// 0 Read in a picture
	X.readImagePixmap("", "Image", color);
	if (trans) X = X.transpose();
	O = Matrix(X); // save original
	R = X.numRows();
	C = X.numCols();
	X.printSize();

	// 1 Center the data
	Matrix meanByCol("Mean");
	meanByCol = X.meanRowVectors();
	meanByCol.printSize();

	X.subRowVector(meanByCol);

	// 2 Compute covariance matrix
	M = X.Tdot(X).scalarMul(1.0f / (float)R);

	// 3 Compute eigenvalues and eigenvectors of M
	V = Matrix(M);
	W = V.eigenSystem();
	V.setName("EigenVectors");
	W.setName("EigenValues");
	V.printSize();
	W.printSize();

	// 4 Normalize the eigenvectors
	// (should be done because of eigenSystem()

	// 5 Sort the eigenvectors by eigenvalue
	// should already be sorted. so we should just truncate
	V = V.extract(0, 0, K, 0);
	W = W.extract(0, 0, 0, K);
	
	// 6 Translate the normalized data
	X = X.dotT(V);
	X.setName("Encoded");
	X.printSize();

	// --- this is where machine learning stuff goes ---


	// 7 Recovering Data from Compressed Data
	X = X.dot(V);				// rotate data back
	X.addRowVector(meanByCol);	// move data back (uncenter)
	X.setName("Decoded");
	X.printSize();

	// 8 The Component Matrix
	// zzz
	
	// calculate RMSE
	float rmse = 0.0f;
	for (int r = 0; r < R; r++) {
		for (int c = 0; c < C; c++) {
			rmse += pow(X.get(r, c) - O.get(r, c), 2);
		}
	}
	rmse /= (float)(R * C);
	rmse = sqrt(rmse);


	// Output
	if (trans) {
		O = O.transpose();
		X = X.transpose();
	}
	if (color) X.writeImagePpm("z.ppm", "output");
	else X.writeImagePgm("z.pgm", "output");

	printf("Min/Max/Mean Value per Pixel of Original: %f %f %f\n", O.min(), O.max(), O.mean());
	printf("Min/Max/Mean Value per Pixel of Recovered: %f %f %f\n", X.min(), X.max(), X.mean());
	printf("Root Mean Squared Error: %f\n", rmse);
}
