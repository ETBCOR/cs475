/*****************************
 * By:			Ethan Corgatelli
 * File:		kmeans.cpp
 * Project:		Assignment 7
 * Class:		cs475
 * Asn. Pg:		http://marvin.cs.uidaho.edu/Teaching/CS475/pas07.pdf
 *
 *****************************/
#include "kmeans.h"

int main (int argc, char **argv) {
	// check number of args
	if (argc-1 != 2) {
		printf("Not enough arguments. Provide two integers: k (amount of clusters), and t (amount of trials).\n");
		return -1;
	}

	// get k and t
	const int k = atoi(argv[1]);
	const int t = atoi(argv[2]);

	// check k
	if (k < 1) {
		printf("k (the first argument) is too small. We must use at least one clusture.\n");
		return -1;
	}

	// check t
	if (t < 1) {
		printf("t (the second argument) is too small. We must do at least one trial.\n");
		return -1;
	}
	
	initRand();

	// read in the data
	Matrix data("kdtree");
	data.read();

	// add col to the data for which class it belongs to
	data.widen(data.numCols() + 1, -1.0f);

	// arrays to store each trial in
	Matrix trials[t];
	float scores[t];
	
	// do multiple trials
	for (int trial = 0; trial < t; trial++) {

		// make a version of the data that doesn't have the extra col
		Matrix dataClean("dataClean");
		dataClean = data.subMatrix(0, 0, 0, data.numCols() - 1);
		Matrix means(k, data.numCols()-1, "Pts");
		means.sampleWithoutRows(dataClean);

		//means.print(); data.print();

		bool moving = true; int n = 0;
		do { // main loop
			//placeholders for points
			Matrix point(1, data.numCols(), "point");
			Matrix mean(1, data.numCols(), "mean");
			// assign each point to a clustre
			for (int i = 0; i < data.numRows(); i++) {
				point = data.subMatrix(i, 0, 1, means.numCols());
				float best = FLT_MAX; int index = -1;
				for (int j = 0; j < means.numRows(); j++) {
					mean = means.subMatrix(j, 0, 1, 0);
					float d = point.dist(mean);
					if (d < best) {
						best = d;
						index = j;
					}
				}
				data.set(i, data.numCols() - 1, index);
			}

			// create a spot for the new means
			Matrix newMeans(k, means.numCols(), 0.0f, "newMeans");

			// move means
			Matrix c(1, data.numCols()-1, "c");
			for (int i = 0; i < k; i++) {
				c = data.subMatrixEq(data.numCols()-1, i);
				if (c.numRows() == 0) {
					c.sampleWithoutRows(data);
				}
				else c = c.meanRowVectors();
				newMeans.insert(c, i, 0);
			}

			if (means.dist(newMeans) == 0) {
				moving = false;
			} else n++;

			means = Matrix(newMeans);
		
		} while (moving);

		Matrix points(1, data.numCols(), "points");
		Matrix point(1, data.numCols(), "point");

		float avgDist = 0.0f;
		for (int i = 0; i < means.numRows(); i++) {
			points = data.subMatrixEq(data.numCols()-1, i);

			float d = 0.0f;
			for (int j = 0; j < points.numRows(); j++) {
				point = points.subMatrix(j, 0, 1, points.numCols()-1);
				d += point.dist2(means.extract(i, 0, 1, 0));
			}
			avgDist += d / points.numRows();
		}
		//grp.sub(mean).dist/magnitude.sum / numPoints

		means.sortRows();
		trials[trial] = means;
		scores[trial] = avgDist;
		
		printf("Num Tries: %d\n", n);
		means.print();
	
		printf("total average dist: %f\n", avgDist);
	}

	float bD = FLT_MAX;
	int bT = -1;
	for (int i = 0; i < t; i++) {
		if (scores[i] < bD) {
			bD = scores[i];
			bT = i;
		}
	}

	trials[bT].print();
	printf("best average dist: %f\n", bD);
}
