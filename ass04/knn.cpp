/*****************************
 * By:			Ethan Corgatelli
 * File:		knn.cpp
 * Project:		Assignment 4
 * Class:		cs475
 * Asn. Pg:		http://marvin.cs.uidaho.edu/Teaching/CS475/pas04.pdf
 *
 *****************************/
#include "knn.h"

int main (int argc, char *argv[]) {

	// check that args are good
	if (argc-1 != 1) {
		printf("%d is not a valid number of inputs. Please provide a single integer for the number of nearest neighbors to search for.\n", argc-1);
		return -1;
	}

	// get K from command line args
	int K = atoi(argv[1]);
	if (K < 1) {
		printf("%d is not a valid input value for K. We must return at least one nearest neighbor.\n", K);
		return -1;
	}

	Matrix list("list");
	Matrix test("test");
	SymbolNumMap lables;

	list.readLabeledRow(&lables);
	test.read();

	for (int i = 0; i < test.numRows(); i++) {
		printf("\n");
		Matrix dist(0, 1, 0.0f);
		Matrix t = test.extract(i, 0, 1, 0);
		t.print("SOLVE:", false);

		for (int j = 0; j < list.numRows(); j++) {
			Matrix l = list.extract(j, 1, 1, 0);
			double d = t.dist(l);
			dist.lengthen(dist.numRows() + 1);
			dist.set(dist.numRows()-1, 0, d);
		}

		Matrix c("dist");
		c = dist.joinRight(list);

		c.sortRowsByCol(0);
		c.shorten(K);
		c.printLabeledRow(&lables, 1, "Sorted", true);
	}
}
