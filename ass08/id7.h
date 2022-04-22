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
//  DATA FOR CONSTRUCTING DECISION TREE
//
//  The answer is always the last row in features and last column in data.
//  The data is in two matrices: features and then data which looks like:
//  The feature matrix begins with the name of the feature in column 0.
//  The last feature (row) in the feature matrix is the answer.
//  If there are fewer featurevalues than the maximum number of columns
//      then the row is filled with "-" strings.
//  Terms:
//     feature is the name of the feature.
//           Appears in col 0 of feature matrix.
//     featurevalue is the values that the feature can take on
//           Appears in col 1, 2, ... of feature matrix.
//     ansvalue is the featurevalues for the last row of the feature matrix
//     featureColumn is the row in the feature matrix and column in the data
//           where a feature can be found.
//
//
//  Data in the COL k of the data matrix must be from the list of featurevalues
//     for ROW k in the feature matrix.  Therefore the lastcol is
//     is the answer column and is composed of featurevalues from
//     the last row of the feature matrix.
//
// numFeatures maxFeatureValues
// feature fval1 fval2 fval3...
//  ...
// feature fval1 fval2 fval3...
// Ans fval1 fval2 fval3...
// numrows numFeatures
// data data data data data data 
// data data data data data data 
// data data data data data data 
// data data data data data data 

// // // // // // // // // // // // // // // // // // // // // // // // 
//
// ID3 Algorithm to build a decision tree
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


// BUILD DECISION TREE
// 
// returns a tree. Tree has nodes and edges each of which can be
// labeled with a string.
// 
// If a node is a parent node then
//     label of node is the name of a
//     feature and edges are labeled with featurevalues
// else node is a leaf
//     label of node is an ansvalue
//
// The algorithm has 4 cases:
// A) all the answers in data are the same -> produce leaf with ans
// B) run out of features -> produce leaf with voting on ans
// C) if feature has only a single featurevalue exhibited in data -> remove that feature from availCol
// Da) for each of the featurevalues produce a subtree
// Db) if featurevalues not present then provide a default by voting

Tree *build(const Matrix &data,
            const Matrix &features,
            const SymbolNumMap *syms,   
            const Matrix &availCol);         // available features to ask about



// look up in decision tree
string find(Tree *tree,
            const SymbolNumMap *syms,
            const SymbolNumMap *fNames,
            const Matrix &query,
            int r);
