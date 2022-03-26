/*****************************
 * By:			Ethan Corgatelli
 * File:		kdtree.cpp
 * Project:		Assignment 5
 * Class:		cs475
 * Asn. Pg:		http://marvin.cs.uidaho.edu/Teaching/CS475/pas05.pdf
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

// BUILD THE KD-TREE
// sort by column c
// median is middle of sorted list
// call build with left side and c+1
// call build with right side and c+1


//SEARCH THE KD-TREE
//kdtree(int bestrow, int bestdistance) returns bestrow and bestdistance
// if leaf node
// if better bestrow save that as new best node
// else if parent node case
// if item left of split point?
// do left then right
// call search(best) on left -> bestrow and best
// if dist(item, split point) in dimension c > best then return
// call search(best) on right
// else
// do right then left
// call search(best) on right -> bestrow and best
// if dist(item, split point) in dimension c > best then return
// call search(best) on left
// do split point (parent)
// if better bestrow save that as new bestnode

	}
}
