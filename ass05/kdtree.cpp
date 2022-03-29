/*****************************
 * By:			Ethan Corgatelli
 * File:		kdtree.cpp
 * Project:		Assignment 5
 * Class:		cs475
 * Asn. Pg:		http://marvin.cs.uidaho.edu/Teaching/CS475/pas05.pdf
 *
 *****************************/
#include "kdtree.h"

int main (int argc, char *argv[]) {

	Matrix tree("tree");
	Matrix test("test");
	SymbolNumMap lables;

	tree.readLabeledRow(&lables);
	test.read();

// BUILD THE KD-TREE
	build(&tree, 1, 0, tree.numRows() - 1, 0);
	//tree.printLabeledRow(&lables);

//SEARCH THE KD-TREE
	int *bestrow;
	float *bestdist;
//kdtree(int bestrow, int bestdistance) returns bestrow and bestdistance
	search(&tree, bestrow, bestdist, 0, tree.numRows() - 1);
}

// BUILD THE KD-TREE
void build (Matrix *t, int c, int lower, int upper, int i) {
	int size = upper-lower+1;
	//printf("#%d (%d-%d) size: %d\n", i, lower, upper, size);

	// base case
	if (size <= 2 || i > 50) return;

	// sort by column c
	t->sortRowsByCol(c, lower, upper);

	// median is middle of sorted list
	int m = size / 2.0f + lower;
	//printf("m: %d\n", m);

	c = (c == 3) ? 1 : c + 1; // increment c

	// call build with left side and c+1
	//printf("do left (%d-%d):\n", lower, m - 1);
	build(t, c, lower, m - 1, i + 1);
	// call build with right side and c+1
	//printf("do right (%d-%d):\n", m + 1, upper);
	build(t, c, m + 1, upper, i + 1);
}

//SEARCH THE KD-TREE
void search(Matrix *t, Matrix x, int *bestrow, float *bestdist, int lower, int upper) {
	int size = upper-lower+1;
	// if leaf node
	if (size <= 0) {
		return;
	} else if (size <= 1) {
		check();
	} else if (size <= 2)
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

int check (Matrix *t, Matrix x, int *bestrow, float *bestdist, int r) {
	Matrix ch(); // best

	b = 
	x = t.subMatrix(r, 1, 1, 0);
	if (x.dist(b) < *bestrow) {
		
	}
	return 0;
}
