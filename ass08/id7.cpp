/*****************************
 * (Modified) By:       Ethan Corgatelli
 * File:                id7.cpp
 * Project:             Assignment 8
 * Class:               cs475
 * Asn. Pg:             http://marvin.cs.uidaho.edu/Teaching/CS475/pas08.pdf
 *****************************/

#include "id7.h"

// GLOBAL feature stuff
int ansCol;                      // column with answer in it (FOR READABILITY)
int missingValue;                // num in SymbolNumMap of missingValue string
int maxFeatureValues;            // the most number of different feature values+1


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
int isOneLabel(int c, const Matrix &d)
{
    double ans;

    ans = d.minCol(c);
    if (ans == d.maxCol(c)) return int(ans);
    else return -1;
}



// given data and features -> entropy in the answer column.
// The less diverse the answers in the set of data the lower the entropy.
// The maximum value of entropy is Log2(n) where n is the number of
// different answers in the Ans column of the given data.
double entropy(const Matrix &data, const Matrix &features)
{
    double sum, p;
    int v;  // symbol number of feature value

    sum = 0.0;
    for (int i=1; i<maxFeatureValues; i++) {
        v = features.get(ansCol, i);   // get feature label i
        if (v==missingValue) break;    // if "-" then quit

        p = data.countEqCol(ansCol, v)/double(data.numRows());
        if (p>0) sum += -p * log2(p);
    }

    return sum;
}


// given data, features, and a feature col number return information gain
double gain(const Matrix &data, const Matrix &features, int fcol)
{
    double sum, p;
    int v;  // symbol number of feature value

    // look through feature values for feature fcol starting with col=1
    // since col=0 is the name of the feature.  Data is subseted by feature f
    sum = 0.0;
    for (int i=1; i<maxFeatureValues; i++) {
        v = features.get(fcol, i);
        if (v==missingValue) break;

        p = data.countEqCol(fcol, v)/double(data.numRows());
        if (p>0) sum += p * entropy(data.subMatrixEq(fcol, v), features);
    }

    return entropy(data, features) - sum;  // compute the gain
}



// given the data find the most popular answer value searching
// in the order of the feature values.  Return symbol number of featurevalue
int vote(const Matrix &data, const Matrix &features) {
    int bestCnt, bestV, cnt;

    printf("VOTE\n");

    // look through all featuresvalue of Ans start with col=1
    bestCnt = -1;
    for (int i=1; i<maxFeatureValues; i++) {
        int v;  // symbol number of feature value

        v = features.get(ansCol, i);
        if (v==missingValue) break;

        cnt = data.countEqCol(ansCol, v);
        if (cnt>bestCnt) {
            bestCnt = cnt;
            bestV = v;
        }
    }

    return bestV;
}


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
            const Matrix &availCol)          // available features to ask about
{
    double ans;
    Tree *result;
    
    // IF NOTHING LEFT -> system error
    if (data.numRows()==0) {                        
        features.print("SYSTEM ERROR. FEATURES");
        result = new Tree("NONE");  // should NEVER get here
    }

    // ELSE IF ONLY ONE TYPE OF ANSWER FOR THIS BRANCH
    else if ((ans = isOneLabel(ansCol, data)) != -1) {      
        printf("TREE A\n");
        // BUILD LEAF
        // [ ADD CODE HERE ]
    }

    // ELSE IF RUN OUT OF FEATURES TO TEST -> VOTE
    else if (availCol.numRows()==0) {           
        printf("TREE B\n");
        // BUILD LEAF WITH VOTE OVER data
        // [ ADD CODE HERE ]
    }

    // ELSE PICK BEST GAIN AND PROCEED WITH SUBTREES
    else {
        double g, bestGain;
        int bestFeatureCol;
        int fcol;
        
        // LOOP THROUGH FEATURE COLUMNS find bestFeatureCol with biggest gain
        // [ ADD CODE HERE ]
        printf("[ ADD CODE HERE ]");

        // IGNORE FEATURE WITH ONLY ONE FEATURE VALUE PRESENTED
        if (isOneLabel(bestFeatureCol, data)>=0) {
            printf("TREE C\n");
            // remove the bestFeatureCol and RECURSE
    	    // [ ADD CODE HERE ]
	        printf("[ ADD CODE HERE ]");
        }
        else {
            printf("TREE D\n");
            // BUILD PARENT NODE FOR TWO OR MORE CHILDREN FEATURE VALUES
    	    // [ ADD CODE HERE ]
	        printf("[ ADD CODE HERE ]");

            // ADD TWO OR MORE CHILDREN
            bool allFound=true;
            for (int i=1; i<maxFeatureValues; i++) {
                int v; // feature value
    			    // [ ADD CODE HERE ]
			        printf("[ ADD CODE HERE ]");

                    // ADD FEATUREVALUE CHILD IF THERE ARE NONZERO AMOUNT OF DATA
    			    // [ ADD CODE HERE ]
			        printf("[ ADD CODE HERE ]");
            }

            // IF NOT ALL FEATUREVALUES REPRESENTED ADD IN DEFAULT CHILD AS VOTE FROM DATA
            if (! allFound) {
                printf("TREE E\n");
   			    // [ ADD CODE HERE ]
		        printf("[ ADD CODE HERE ]");
            }
        }
    }

    return result;
}



// look up in decision tree
string find(Tree *tree,
            const SymbolNumMap *syms,
            const SymbolNumMap *fNames,
            const Matrix &query,
            int r) {
    int c;
    Tree *next;

    // LEAF: THEN RETURN NODE NAME
    if (tree->isLeaf()) {
	    // [ ADD CODE HERE ]
        printf("[ ADD CODE HERE ]");
    }

    // NOT LEAF
    else {
        // TRY TO FOLLOW EDGE TO SUBTREE
        c = fNames->getNum(tree->getName());  // get column of the feature
        next = tree->getChild(syms->getStr(query.get(r, c)));
        if (next!=NULL) {
		    // [ ADD CODE HERE ]
	        printf("[ ADD CODE HERE ]");
        }
        // IF NO EDGE FOR FEATURE VALUE THEN GET DEFAULT VALUE
        else {
		    // [ ADD CODE HERE ]
	        printf("[ ADD CODE HERE ]");
        }
    }
}


int main()
{
    Matrix features("Features");
    Matrix data("Data");
    SymbolNumMap *syms = new SymbolNumMap("Symbol Map");
    Tree *tree;

    // READ FEATUES
    missingValue = syms->add("-");      // missing values are symbol missingValue.  Save globally
    syms = features.readStrings(syms); // features list in column 0

    // READ DATA
    syms = data.readStrings(syms);     // each feature in a different column
    data.assertOtherRhs(features, "main");     // check dimensions match
    ansCol = data.numCols()-1;         // NOTE: answer guaranteed in last column!  Save globally
    maxFeatureValues = features.numCols(); // maximum number of feature values anywhere.  Save globally

    // CREATE LIST OF UNUSED FEATURES
    Matrix availCol(ansCol, 1, "Available Columns");  // a column vector
    availCol.initLinear(1, 0, 0);             // list of available feature columns


	features.printLabeledRow(syms);
	data.print();


/*
    // BUILD DECISION TREE
    tree = build(data, features, syms, availCol);
    tree->printWithEdges();

    // SEARCH DECISION TREE
    Matrix query("Query");
    query.readStrings(syms);

    // a map starting with number 0 for first symbol
    SymbolNumMap *fNames = new SymbolNumMap("Feature Names", 0);
    for (int i=0; i<features.numRows(); i++) {
        // map  i <--> ith feature name
        fNames->add(syms->getStr(features.get(i, 0))); 
    }

    for (int r=0; r<query.numRows(); r++) {
        printf("%d ", r);
        query.writeLineStrings(syms, r);
        printf("%s\n", find(tree, syms, fNames, query, r).c_str());
    }
*/
    return 0;
}
