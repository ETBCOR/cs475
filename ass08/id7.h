/*****************************
 * (Modified) By:       Ethan Corgatelli
 * File:                id7.h
 * Project:             Assignment 8
 * Class:               cs475
 * Asn. Pg:             http://marvin.cs.uidaho.edu/Teaching/CS475/pas08.pdf
 *****************************/

#include <iostream>
#include <stdlib.h>

#include "tree.h"
#include "mat.h"
#include "rand.h"

// // // // // // // // // // // // // // // // // // // // // // // // 
//
// ID3 Algorithm to build a decision tree (header)
//

// returns the ONE common label (a symbol number) if there is one else returns -1
// returns symbol number
int isOneLabel(int c, const Matrix &d);

// given data and features -> entropy in the answer column.
// The less diverse the answers in the set of data the lower the entropy.
// The maximum value of entropy is Log2(n) where n is the number of
// different answers in the Ans column of the given data.
double entropy(const Matrix &data, const Matrix &features);

// given data, features, and a feature col number return information gain
double gain(const Matrix &data, const Matrix &features, int fcol);

// given the data find the most popular answer value searching
// in the order of the feature values.  Return symbol number of featurevalue
int vote(const Matrix &data, const Matrix &features);

// returns a tree. Tree has nodes and edges each of which can be
// labeled with a string.
Tree *build(const Matrix &data,
            const Matrix &features,
            const SymbolNumMap *syms,   
            const Matrix &availCol);   // available features to ask about

// look up in decision tree
string find(Tree *tree,
            const SymbolNumMap *syms,
            const SymbolNumMap *fNames,
            const Matrix &query,
            int r);
