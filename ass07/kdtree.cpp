/*****************************
 * By:			Ethan Corgatelli
 * File:		kdtree.cpp (improved output)
 * Project:		Assignment 6
 * Class:		cs475
 * Asn. Pg:		http://marvin.cs.uidaho.edu/Teaching/CS475/pas06.pdf
 *
 *****************************/
#include "kdtree.h"

/*
int main (int argc, char *argv[]) {

	Matrix tree("kdtree");
	Matrix test("test");
	SymbolNumMap lables;

	tree.readLabeledRow(&lables);
	test.read();

// BUILD THE KD-TREE
	build(&tree, 1, 0, tree.numRows() - 1, 0);
	tree.printLabeledRow(&lables);

//SEARCH THE KD-TREE
	for (int i = 0; i < test.numRows(); i++) {
		
		Matrix x = test.subMatrix(i, 0, 1, 0);
		x.printFmt("\nFIND:", "%.0f ", false);

		int bestrow = -1;
		float bestdist = FLT_MAX;
		
		search(&tree, x, &bestrow, &bestdist, 0, tree.numRows() - 1, n);
		printf("Num Compares: %d\n", n);

		Matrix r(1, 1, bestdist, "r");
		x = tree.extract(bestrow, 0, 1, 0);
		r = r.joinRight(x);
		r.printLabeledRow(&lables, 1, "", false);
	}
*/
}

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

	c = (c >= 3) ? 1 : c + 1; // increment c

	// call build with left side and c+1
	//printf("do left (%d-%d):\n", lower, m - 1);
	build(t, c, lower, m - 1, i + 1);
	// call build with right side and c+1
	//printf("do right (%d-%d):\n", m + 1, upper);
	build(t, c, m + 1, upper, i + 1);
}

void search(Matrix *t, Matrix x, int *bestrow, float *bestdist, int lower, int upper, int c) {
	int size = upper-lower+1;
	int m = size / 2.0f + lower;
	c = (c >= 3) ? 1 : c + 1; // increment c

	//printf("br: %d, bd: %lf, l: %d, u: %d, size: %d, m: %d, c: %d\n", *bestrow, *bestdist, lower, upper, size, m, c);
	
	// if leaf node
	if (size <= 0) {
		return;
	} else if (size <= 1) {
		check(t, x, bestrow, bestdist, upper);
	}
	else if (size <= 2) {
		check(t, x, bestrow, bestdist, upper);
		check(t, x, bestrow, bestdist, lower);
	}
	// else if parent node case
	else {
		// if item left of split point?
		if (x.get(0, c-1) < t->get(m, c)) {
			// do left then right
			// call search(best) on left -> bestrow and best
			search (t, x, bestrow, bestdist, lower, m - 1, c);
			// if dist(item, split point) in dimension c > best then return
			//if (check(t, x, bestrow, bestdist, m)) return;
			if (abs(t->get(m, c) - x.get(0, c-1)) > *bestdist) return;
			// call search(best) on right
			search (t, x, bestrow, bestdist, m + 1, upper, c);
		} else {
			// do right then left
			// call search(best) on right -> bestrow and best
			search (t, x, bestrow, bestdist, m + 1, upper, c);
			// if dist(item, split point) in dimension c > best then return
			if (abs(t->get(m, c) - x.get(0, c-1)) > *bestdist) return;
			// call search(best) on left
			search (t, x, bestrow, bestdist, lower, m - 1, c);
		}
		// do split point (parent)
		// if better bestrow save that as new bestnode
		check(t, x, bestrow, bestdist, m);
	}
	return;
}

int check (Matrix *t, Matrix x, int *bestrow, float *bestdist, int r) {
	Matrix ck = t->subMatrix(r, 1, 1, 0);

	// if better bestrow save that as new best node
	if (x.dist(ck) < *bestdist) {
		*bestrow = r;
		*bestdist = x.dist(ck);
		return 1;
	}
	return 0;
}
