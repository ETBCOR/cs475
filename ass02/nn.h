/*****************************
 * By:			Ethan Corgatelli
 * File:        nn.h
 * Project:     Assignment 1 (improved!)
 * Class:       cs475
 * Asn. Pg:     http://marvin.cs.uidaho.edu/Teaching/CS475/pas01.pdf
 *
 *****************************/
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "mat.h"
#include "rand.h"

class NeuralNet {
//friend class Layer;

public:
	static bool debug; // debugging flag

private:
		// constants
	static double eta;
	static double bias;
	static double spread;

	int N;	// size of 1 input vector
	int M;	// size of 1 output vector

	int nodes;

		// read-in matrices
	Matrix *trI; // training input
	Matrix *trO; // training output

		// for training
	Matrix *X;
	Matrix *V;
	Matrix *H;
	Matrix *W;
	Matrix *Y;
	Matrix *T;

		// matrices output
	Matrix *r; // recursivly built row by row
	Matrix *m; // final output matrix

public:
	// constructor
	NeuralNet(int _nodes);

	// read/prep input from stdin
	void read();

	// train for a specified amount of times
	void train(int passes = 1);

	// test and print results
	void guess();

};

double transfer (double v);
