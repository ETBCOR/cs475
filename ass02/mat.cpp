// // // // // // // // // // // // // // // //
//
// Some simple matrix operators
//
// These routines are self contained except for random number
// generation.   Comment these out if you want completely stand-alone
// functionality.  Uses some code from Numerical Recipes (see comments).
// WARNING: if you use the included random number generator then be sure
// to call initRand() before using any random number function!
//

// // // // // // // // // // // // // // // // // // // // // // // // // 
// Author: Dr. Robert B. Heckendorn, University of Idaho
// Version: 7.2WX
// Date: Feb 1, 2022
//
// You are free to use, change, or redistribute the code in any way
// you wish for non-commercial purposes, but please maintain the name
// of the original author. This code is not guaranteed to function
// correctly and comes with no warranty of any kind for any purpose
// what-so-ever and is NOT SUPPORTED.   Good luck and have fun.

// Matrices are allocated as an array of rows where are arrays.
// The rows and the array of rows can be allocated separately
// Matrices can be:

//  1) UNALLOCATED in which case the dimensions maxr, maxc are both -1.

//  2) The array of rows alone or the array of rows and all the rows
//  can be ALLOCATED in which case the dimensions in maxr, maxc are the
//  size allocated.

//  3) The array can be allocated and DEFINED, that is have values for
//  all its elements.

#include "mat.h"
// can do this in class in C++11
char *Matrix::realFormat=(char *)"%8.4lf ";
char *Matrix::intFormat=(char *)"%8d ";
char *Matrix::shortIntFormat=(char *)"%5d ";
char *Matrix::binaryFormat=(char *)"%s ";

// the following are routines taken from the book Numerical Recipes in C
static void householder(double **a, int n, double d[], double e[]);
static void eigen(double *d, double *e, int n, double **z);
static bool gaussj(double **a, int n, double **b, int m);

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
//
//  class SymbolNumMap
//
//  a class for a 1-1 mapping std::string to int. Useful if you want
//  to read in data as strings and force them into a numeric matrix.

// constructor
SymbolNumMap::SymbolNumMap(string namex, int startnum, string defaultValue) :
    name(namex), startNum(startnum), defaultMissing(defaultValue), next(startnum){
};


// add a symbol to the map if it isn't there otherwise look it up
int SymbolNumMap::add(string s) {
    if (strToNum.find(s) == strToNum.end()) {
        strToNum[s] = next;
        numToStr[next] = s;
        return next++;
    }
    return strToNum.find(s)->second;
}

// add a symbol to the map if it isn't there otherwise look it up
// captures a char* for error proofing of NULL
int SymbolNumMap::add(const char *s) {
    if (s==NULL) {
        printf("ERROR(SymbolNumMap::add): Trying to create a std:string from a char * that is NULL\n");
        exit(1);
    }
    else {
        return add(std::string(s));
    }
}


// get the int map to the given string if it is there otherwise return -1
int SymbolNumMap::getNum(const string s) const {
    if (strToNum.find(s) == strToNum.end()) {
        return -1;
    }
    return strToNum.find(s)->second;
}


// get the int map to the given char* if it is there otherwise return -1
// error if the char* is NULL
int SymbolNumMap::getNum(const char *s) const {
    if (s==NULL) {
        printf("ERROR(SymbolNumMap::getNum): Trying to create a std:string from a char * that is NULL\n");
        exit(1);
    }
    else {
        return getNum(std::string(s));
    }
}


// get the string associate with the int or return default string
string SymbolNumMap::getStrDefault(int n) const {
    if (n<startNum || n>=next) {
        printf("ERROR(SymbolNumMap::getStrDefault): index into labels must be an integer in the range %d to %d but is instead %d.\n", startNum, next-1, n);
        exit(1);
    }
    if (numToStr.find(n) == numToStr.end()) {
        return defaultMissing;
    }
    return numToStr.find(n)->second;
}


// get the string associate with the double or error for one of several reasons
string SymbolNumMap::getStr(double x) const {
    int n = x;
    if (n != x) {
        printf("ERROR(SymbolNumMap::getStr): index into labels must be an integer but is %lg instead.\n", x);
    exit(1);
    }
    else if (n<startNum || n>=next) {
        printf("ERROR(SymbolNumMap::getStr): index into labels must be an integer in the range %d to %d but is instead %d.\n", startNum, next-1, n);
        exit(1);
    }
    else if (numToStr.find(n) == numToStr.end()) {  
        return defaultMissing;   // should never get here
    }
    return numToStr.find(n)->second;
}

void SymbolNumMap::clear() {
    strToNum.clear();
    numToStr.clear();
    next = startNum;
}

void SymbolNumMap::print(string msg){
    if (msg.length()) {
        printf("%s ", msg.c_str());
    }

    if (name.length()) {
        printf("(size of SymbolNumMap %s: %d)\n", name.c_str(), next-startNum);
    }
    else {
        printf("(size of SymbolNumMap: %d\n", next-startNum);
    }

    for (int i=startNum; i<next; i++) {
        printf("%3d %s\n", i, numToStr[i].c_str());
    }
    fflush(stdout);
}

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
//
//  class MatrixRowIter
//
//  a friend class to Matrix that iterates through a Matrix: a row iterator
//
//  Note: although it uses std::string it uses printf because I hate cout's
//  awkwardness. :-)
//  Note: errors are reported to stdout
//

// this allocates the space for the row
MatrixRowIter::MatrixRowIter(Matrix *newmat)
{
    mat = newmat;
    r = 0;
    arow = new Matrix(1, mat->maxc, "row of " + newmat->name);  // allocate the space for the row matrix
    more = true;
}


// this deallocates the row space
MatrixRowIter::~MatrixRowIter()
{
    mat = NULL;
    r = 0;
    delete arow;                  // deallocate the row pointed to by arow
    more = false;
}



// by row iterator
Matrix *MatrixRowIter::rowBegin()
{
    mat->assertDefined("MatrixRowIter");

    r = 0;
    for (int i=0; i<mat->maxc; i++) arow->m[0][i] = mat->m[r][i];
    more = true;
    arow->defined = true;
    arow->submatrix = false;

    return arow;
}


Matrix *MatrixRowIter::rowNext()
{
    if (r < mat->maxr-1) {
        r++;
        for (int i=0; i<mat->maxc; i++) arow->m[0][i] = mat->m[r][i];
        arow->defined = true;
        arow->submatrix = false;
    }
    else {
        more = false;
    }
    return arow;
}


bool MatrixRowIter::rowNotEnd()
{
    return more;
}


int MatrixRowIter::row()
{
    return r;
}



// // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
//
//  class Matrix
//
// A simple class for matrix operations.   It has some nice debugging features like
// trying hard to check that the proper dimensions are used.   It is very draconian
// about this so there is an important difference between row vectors and column vectors.
// I find this helps students get the math to work out correctly if you pay attention to
// this difference.   The routines allow you to name a matrix.  The name is then used
// in debug output.   Other things checked include referencing out of bounds.
// If a row is a nonnegative number then r row pointers will be allocated.
// If a row and col are nonnegative numbers then c columns will be allocated for each row.
// If a row is nonnegative but col is negative then only row pointers allocated which is used
// for submatrices.
// Note: if row<0 no space at all will be allocated.
//
// Default: isSubMatrix=false
//

// turn on debugging here.   It shows when a matrix is allocated and then deallocated.
bool Matrix::debug = false;

void Matrix::allocate(int r, int c, std::string namex, bool isSubMatrix)
{
    if (isSubMatrix) c=-1;    // how you signal you are allocating a submatrix
    maxr = r;
    maxc = c;
    name = namex;
    m = NULL;

    if (maxr < 0 || maxc < 0) {
        if (maxr != -1 && maxc !=-1) {
            if (name.length()==0)
                printf("ERROR(allocate): Trying to create a matrix of size %d X %d\n", r, c);
            else
                printf("ERROR(allocate): Trying to create matrix \"%s\" of size %d X %d\n",
                       name.c_str(), r, c);
            exit(1);
        }
    }

    if (maxr>=0) {
        m = new double * [maxr];
        if (maxc>=0) for (int i=0; i<maxr; i++) m[i] = new double [maxc];  // allocate if not a submatrix
    }

    defined = false;
    submatrix = isSubMatrix;
    if (debug) printf("DEBUG(  allocate): name \"%s\", size %d X %d\n", name.c_str(), r, c);
}


bool Matrix::deallocate()
{
    bool allocated;

    allocated = (m!=NULL);
    if (allocated) {
        if (!submatrix) for (int i=0; i<maxr; i++) delete [] m[i];
        delete [] m;
        m = NULL;   // to be sure
    }

    if (debug) printf("DEBUG(deallocate): name \"%s\", size %d X %d\n", name.c_str(), maxr, maxc);

    defined = false;
    submatrix = false;
    maxr = -1;
    maxc = -1;

    return allocated;
}


// zzz can be made more efficient!!!
void Matrix::reallocate(int otherMaxr, int otherMaxc, std::string namex)
{
    if (maxr!=otherMaxr || maxc!=otherMaxc) {
        if (debug) printf("DEBUG(reallocate): name \"%s\", size %d X %d\n", name.c_str(), otherMaxr, otherMaxc);
        deallocate();
        allocate(otherMaxr, otherMaxc, namex);
    }
    if (namex!="") name = namex;
}


Matrix::Matrix(std::string namex)
{
    allocate(-1, -1, namex);   // allocate no size
    m = NULL;
}



// This creates a subMatrix so only allocated the rows!
Matrix::Matrix(int r, std::string namex)
{
    allocate(r, -1, namex, true);       // WARNING: allocate as subMatrix!!!
}


// creates a matrix with a given name
Matrix::Matrix(int r, int c, std::string namex)
{
    allocate(r, c, namex);
}


// creates a matrix with a given name and initialized to a constant value
Matrix::Matrix(int r, int c, double value, std::string namex)
{
    allocate(r, c, namex);
    constant(value);
}



// creates a matrix with a given name and initialized to values given from
// array of data used by first filling each row before moving to the next row
Matrix::Matrix(int r, int c, const double *data, std::string namex)
{
    allocate(r, c, namex);

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = *data++;
        }
    }

    defined = true;
}


// creates a matrix with a given name and initialized to values given from
// array of data used by first filling each row before moving to the next row
Matrix::Matrix(int r, int c, int *data, std::string namex)
{
    allocate(r, c, namex);

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = *data++;
        }
    }

    defined = true;
}


// real copy constructor
Matrix::Matrix(const Matrix &other, std::string namex)
{
//    printf("Matrix Copy Constructor\n");
    other.assertDefined("Matrix Copy Constructor");
    allocate(other.maxr, other.maxc, (namex=="" ? other.name : namex));

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = other.m[r][c];
        }
    }

    defined = true;
}


// pointer version
Matrix::Matrix(Matrix *other)
{
    printf("Warning(Matrix(Matrix *)): Matrix::Matrix(Matrix *other) is depricated.  Consider using a Matrix declaration\n");
    printf("rather than using new to create a matrix pointer.\n");
    other->assertDefined("Matrix Copy Constructor (pointer version)");
    allocate(other->maxr, other->maxc, "");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = other->m[r][c];
        }
    }

    defined = true;
}



Matrix::~Matrix()
{
    deallocate();
}



// WARNING: as written allows direct access to matrix name (not copy)
const std::string &Matrix::getName(const std::string &defaultName) const
{
    if (name.length()==0) return defaultName;
    else return name;
}


// get the value of an element of the matrix
double Matrix::get(int r, int c) const
{
    assertIndexOK(r, c, "get");

    return m[r][c];
}


// increment a given element returning the new value
double Matrix::inc(int r, int c)
{
    assertIndexOK(r, c, "inc");

    return ++m[r][c];
}


// decrement a given element returning the new value
double Matrix::dec(int r, int c)
{
    assertIndexOK(r, c, "dec");
    return --m[r][c];
}


// set the value of an element of the matrix returning the new value
double Matrix::set(int r, int c, double v)
{
    assertIndexOK(r, c, "set");
    m[r][c] = v;

    return v;
}


// set the value of an element of the matrix  (use this carefully!)
void Matrix::setDefined()
{
    defined = true;
}


// set the name of a matrix
void Matrix::setName(std::string newName)
{
    name = newName;
}


// remove trailing columns WITHOUT actually giving up the space.
void Matrix::narrow(int newc)
{
    assertDefined("narrow");
    assertColIndexOK(newc-1, "narrow");  // allow newc to equal maxc

    maxc = newc;
}


void Matrix::widen(int newc, double fill)
{
    assertDefined("widen");
    if (newc<=maxc) {
        if (name.length()==0)
            printf("ERROR(widen): new width %d for matrix is less than or equal to old withd %d\n", newc, maxc);
        else
            printf("ERROR(widen): new width %d for matrix \"%s\" is less than or equal to old withd %d\n", newc, name.c_str(), maxc);
        exit(1);
    }

    for (int r=0; r<maxr; r++) {
        double *newrow;

        newrow = new double [newc];
        for (int c=0; c<maxc; c++) {
            newrow[c] = m[r][c];
        }
        for (int c=maxc; c<newc; c++) {
            newrow[c] = fill;
        }
        delete m[r];
        m[r] = newrow;
    }
    maxc = newc;
}


void Matrix::lengthen(int newr, double fill)
{
    double **newrows;

    assertDefined("lengthen");
    if (newr<=maxr) {
        if (name.length()==0)
            printf("ERROR(lengthen): new length %d for matrix is less than or equal to old withd %d\n", newr, maxr);
        else
            printf("ERROR(lengthen): new length %d for matrix \"%s\" is less than or equal to old withd %d\n", newr, name.c_str(), maxr);
        exit(1);
    }

    newrows = new double * [newr];

    for (int r=0; r<maxr; r++) {
        newrows[r] = m[r];
    }
    for (int r=maxr; r<newr; r++) {
        newrows[r] = new double [maxc];
        for (int c=0; c<maxc; c++) {
            newrows[r][c] = fill;
        }
    }
    delete m;
    m = newrows;
    maxr = newr;
}



double *Matrix::getRowPtr(int r)
{
    assertDefined("getRowPtr");
    assertRowIndexOK(r, "getRowPtr");

    return m[r];
}


double **Matrix::getMatPtr()
{
    assertDefined("getMatPtr");

    return m;
}




// remove trailing rows without actually giving up the space.
// The argument give is the new number of rows in the matrix.
// DANGER: this trims rows.   Not a memory leak but old row length
// not "visibibly" saved.   Freeing up the space should work correctly.
void Matrix::shorten(int newr)
{
    assertIndexOK(newr-1, 0, "shorten");   // allow newr to equal maxr

    maxr = newr;
}


// // // // // // // // // // // // // // // // // // // // //
//
// Assertions
//

void Matrix::assertAllocated(std::string msg) const
{
    if (maxr == -1) {
        if (name.length()==0)
            printf("ERROR(%s): matrix is expected to be allocated (given a size) but is not.\n", msg.c_str());
        else
            printf("ERROR(%s): matrix \"%s\" is expected to be allocated (given a size) but is not.\n", msg.c_str(), name.c_str());
        exit(1);
    }
}


void Matrix::assertDefined(std::string msg) const
{
    if (!defined) {
        if (name.length()==0)
            printf("ERROR(%s): matrix is undefined\n", msg.c_str());
        else
            printf("ERROR(%s): matrix \"%s\" is undefined\n", msg.c_str(), name.c_str());
        exit(1);
    }
}



// assert matrix doesn't have negative dimensions
void Matrix::assertUsableSize(std::string msg) const
{
    if (maxr < 0 || maxc < 0) {
        if (name.length()==0)
            printf("ERROR(%s): Matrix is of unusable size %d X %d\n", msg.c_str(), maxr, maxc);
        else
            printf("ERROR(%s): Matrix \"%s\" is of unusable size %d X %d\n",
                   msg.c_str(), name.c_str(), maxr, maxc);
        exit(1);
    }
}


void Matrix::assertSquare(std::string msg) const
{
    assertDefined(msg);

    if (maxr != maxc) {
        if (name.length()==0)
            printf("ERROR(%s): the matrix is %dX%d and not square as expected!\n",
               msg.c_str(), maxr, maxc);
        else
            printf("ERROR(%s): the matrix \"%s\" is %dX%d and not square as expected!\n",
               msg.c_str(), name.c_str(), maxr, maxc);
        exit(1);
    }
}


// assert size is rxc
void Matrix::assertSize(int r, int c, std::string msg) const
{
    assertDefined(msg);
    if (maxr != r || maxc != c) {
        if (name.length()==0)
            printf("ERROR(%s): the matrix is %dX%d and not %dX%d as expected!\n",
               msg.c_str(), maxr, maxc, r, c);
        else
            printf("ERROR(%s): matrix \"%s\" is %dX%d and not %dX%d as expected!\n",
               msg.c_str(), name.c_str(), maxr, maxc, r, c);
        exit(1);
    }
}


// assert the row is in bounds for matrix
void Matrix::assertRowIndexOK(int r, std::string msg) const
{
    if (r<0 || r>=maxr) {
        if (name.length()==0) {
            printf("ERROR(%s): row index %d out of bounds.  Matrix size is %d X %d\n",
                   msg.c_str(), r, maxr, maxc);
        }
        else {
            printf("ERROR(%s): row index %d out of bounds.  Matrix \"%s\" is %d X %d\n",
                   msg.c_str(), r, name.c_str(), maxr, maxc);
        }
        exit(1);
    }
}

// assert the column is in bounds for matrix
void Matrix::assertColIndexOK(int c, std::string msg) const
{
    if (c<0 || c>=maxc) {
        if (name.length()==0) {
            printf("ERROR(%s): column index %d out of bounds.  Matrix size is %d X %d\n",
                   msg.c_str(), c, maxr, maxc);
        }
        else {
            printf("ERROR(%s): column index %d out of bounds.  Matrix \"%s\" is %d X %d\n",
                   msg.c_str(), c, name.c_str(), maxr, maxc);
        }
        exit(1);
    }
}


// assert r,c is in matrix
void Matrix::assertIndexOK(int r, int c, std::string msg) const
{
    if (r<0 || r>=maxr || c<0 || c>=maxc) {
        if (name.length()==0) {
            printf("ERROR(%s): index out of bounds: asking for (%d, %d) but size is %d X %d\n",
                   msg.c_str(), r, c, maxr, maxc);
        }
        else {
            printf("ERROR(%s): index out of bounds: asking for (%d, %d) but size of matrix \"%s\" is %d X %d\n",
                   msg.c_str(), r, c, name.c_str(), maxr, maxc);
        }
        exit(1);
    }
}



// assert other is the same size
// assert other can be rhs of matrix product op
void Matrix::assertOtherRhs(const Matrix &other, std::string msg) const
{
    if (maxc!=other.maxr) {
        if (name.length()==0 && other.name.length()==0) {
            printf("ERROR(%s): Dimensions do not match: self: %d X %d other: %d X %d\n",  msg.c_str(), maxr, maxc, other.maxr, other.maxc);
        }
        else {
            printf("ERROR(%s): Dimensions do not match: self \"%s\": %d X %d other \"%s\": %d X %d\n", msg.c_str(), name.c_str(), maxr, maxc, other.name.c_str(), other.maxr, other.maxc);
        }
        exit(1);
    }
}


// assert other can be lhs of mult op
void Matrix::assertRowsEqual(const Matrix &other, std::string msg) const
{
    if (maxr!=other.maxr) {
        if (name.length()==0 && other.name.length()==0) {
            printf("ERROR(%s): Row dimensions do not match: self: %d X %d other: %d X %d\n", msg.c_str(), maxr, maxc, other.maxr, other.maxc);
        }
        else {
            printf("ERROR(%s): Row dimensions do not match: self \"%s\": %d X %d other \"%s\": %d X %d\n", msg.c_str(), name.c_str(), maxr, maxc, other.name.c_str(), other.maxr, other.maxc);
        }
        exit(1);
    }
}


// assert other can be lhs of mult op
void Matrix::assertColsEqual(const Matrix &other, std::string msg) const
{
    if (maxc!=other.maxc) {
        if (name.length()==0 && other.name.length()==0) {
            printf("ERROR(%s): Column dimensions do not match: self: %d X %d other: %d X %d\n", msg.c_str(), maxr, maxc, other.maxr, other.maxc);
        }
        else {
            printf("ERROR(%s): Column dimensions do not match: self \"%s\": %d X %d other \"%s\": %d X %d\n", msg.c_str(), name.c_str(), maxr, maxc, other.name.c_str(), other.maxr, other.maxc);
        }
        exit(1);
    }
}


// assert two matrices have the same size
void Matrix::assertOtherSizeMatch(const Matrix &other, std::string msg) const
{
    assertRowsEqual(other, msg);
    assertColsEqual(other, msg);
}


// assert is a row vector
void Matrix::assertRowVector(std::string msg) const
{
    if (maxr!=1) {
        if (name.length()==0) {
            printf("ERROR(%s): expecting matrix is row vector but size is %d X %d\n",
                   msg.c_str(), maxr, maxc);
        }
        else {
            printf("ERROR(%s): expecting matrix %s is row vector but size is %d X %d\n",
                   msg.c_str(), name.c_str(), maxr, maxc);
        }
        exit(1);
    }
}


// assert is a column vector
void Matrix::assertColVector(std::string msg) const
{
    if (maxc!=1) {
        if (name.length()==0) {
            printf("ERROR(%s): expecting matrix is column vector but size is %d X %d\n",
                   msg.c_str(), maxr, maxc);
        }
        else {
            printf("ERROR(%s): expecting matrix %s is column vector but size is %d X %d\n",
                   msg.c_str(), name.c_str(), maxr, maxc);
        }
        exit(1);
    }
}


// assert row size is a power of 2
void Matrix::assertRowPower2(std::string msg) const
{
    assertDefined(msg);

    if ((maxr & (maxr-1)) != 0) {
        if (name.length()==0)
            printf("ERROR(%s): the number of rows in the matrix is %d and not a power of 2 as expected!\n",
               msg.c_str(), maxr);
        else
            printf("ERROR(%s): the number of the rows in the matrix \"%s\" is %d and not a power of 2 as expected!\n",
               msg.c_str(), name.c_str(), maxr);
        exit(1);
    }
}


// assert col size is a power of 2
void Matrix::assertColPower2(std::string msg) const
{
    assertDefined(msg);

    if ((maxc & (maxc-1)) != 0) {
        if (name.length()==0)
            printf("ERROR(%s): the number of columns in the matrix is %d and not a power of 2 as expected!\n",
               msg.c_str(), maxc);
        else
            printf("ERROR(%s): the number of the columns in the matrix \"%s\" is %d and not a power of 2 as expected!\n",
               msg.c_str(), name.c_str(), maxc);
        exit(1);
    }
}


void Matrix::assertRandInitialized(std::string msg) const
{
    if (!isInitRand()) {
        if (name.length()==0)
            printf("ERROR(%s): requires initialization first with a call to initRand()\n", msg.c_str());
            
        else
            printf("ERROR(%s): on matrix \"%s\" requires initialization first with a call to initRand()\n",
                   msg.c_str(), name.c_str());
        exit(1);
    }
}


// helper function that swaps two rows and does not check matrix for validity
void Matrix::swapRows(int i, int j)
{
    double *tmp;

    tmp = m[i]; m[i] = m[j]; m[j] = tmp;
}


// helper function that compares two rows and does no assertions
bool Matrix::isLessRows(int i, int j) const
{
//    assertDefined("isLessRows");

    for (int c=0; c<maxc; c++) {
        if (m[i][c] > m[j][c]) return false;
        if (m[i][c] < m[j][c]) return true;
    }

    return false;
}


// assign a matrix.   Size does NOT have to match
// NOTE: this will copy the matrix values over
// NOTE: this will make a FULL copy of a submatrix (fix?)
Matrix &Matrix::operator=(const Matrix &other)
{
    other.assertDefined("rhs of operator=");

    if (this==&other) return *this;       // avoid self copy

    // allocate if a new size
    reallocate(other.maxr, other.maxc, name);

    // copy
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = other.m[r][c];
        }
    }
    defined = true;

    return *this;
}


#ifdef WALSH
#include "matwalsh.cpp"
#endif


// extracts a matrix from another starting at (minr, minc) and of
// the size given.
// NOTE: zero size means "to the end of row or column"!
// WARNING: allocates new matrix for answer
Matrix Matrix::extract(int minr, int minc, int sizer, int sizec)
{
    if (sizer==0) sizer = maxr - minr;
    if (sizec==0) sizec = maxc - minc;

    assertIndexOK(minr, minc, "lower bounds extract");
    assertIndexOK(minr+sizer-1, minc+sizec-1, "upper bounds extract");

    Matrix out(sizer, sizec);

    for (int r=minr; r<minr+sizer; r++) {
        for (int c=minc; c<minc+sizec; c++) {
            out.m[r-minr][c-minc] = m[r][c];
        }
    }
    out.defined = true;

    return out;
}


// extracts a matrix from another starting at (minr, minc) and of
// the stride length given for rows and cols.
// WARNING: allocates new matrix for answer
Matrix Matrix::extractStride(int minr, int minc, int stepr, int stepc)
{
    int newr, newc;
    assertIndexOK(minr, minc, "lower bounds extract");

    newr = (maxr-minr + (stepr - 1))/stepr + 1;
    newc = (maxc-minc + (stepc - 1))/stepc + 1;
    Matrix out(newr, newc);

    for (int r=minr; r<maxr; r+=stepr) {
        for (int c=minc; c<maxc; c+=stepc) {
            out.m[(r-minr)/stepr][(c-minc)/stepc] = m[r][c];
        }
    }
    out.defined = true;

    return out;
}


// does the same extraction as above but requires that the out Matrix
// be correctly allocated beforehand!  <--- WARNING!
Matrix &Matrix::extract(int minr, int minc, int sizer, int sizec, Matrix &out)
{
    if (sizer==0) sizer = maxr - minr;
    if (sizec==0) sizec = maxc - minc;

    for (int r=minr; r<minr+sizer; r++) {
        for (int c=minc; c<minc+sizec; c++) {
            out.m[r-minr][c-minc] = m[r][c];
        }
    }
    out.defined = true;

    return out;
}


// insert Matrix other at location (minr, minc) in self
// NOTE: anything outside of allocated space will issue a warning but
// not be copied!
Matrix &Matrix::insert(const Matrix &other, int minr, int minc)
{
    assertIndexOK(minr, minc, "insert");

    for (int r=0; r<other.maxr; r++) {
        if (r>=maxr) break;
        for (int c=0; c<other.maxc; c++) {
            if (c>=maxc) break;
            m[r+minr][c+minc] = other.m[r][c];
        }
    }

    return *this;
}


// insert at the specified row the other matrix which is a row vector matrix
Matrix &Matrix::insertRowVector(int loc, const Matrix &other)
{
    other.assertRowVector("insertRowVector");
    assertColsEqual(other, "insertRowVector");
    assertIndexOK(loc, 0, "insertRowVector");

    for (int c=0; c<other.maxc; c++) {
        m[loc][c] = other.m[0][c];
    }

    return *this;
}


// indices is a list of indices to select from data
// A new Matrix is created
Matrix Matrix::indexCols(int *indices, int sizeList)
{
    assertDefined("indexRowVector");

    Matrix out(maxr, sizeList, "");

    for (int i=0; i<sizeList; i++) {
        if (0 <= indices[i] && indices[i] < maxc) {
            for (int r=0; r<maxr; r++) {
                out.m[r][i] = m[r][indices[i]];
            }
        }
        else {
            printf("ERROR(indexCols): index %d is out of bounds for input row vector of size %d\n", indices[i], maxc);

            exit(1);
        }
    }

    out.defined = true;

    return out;
}


// rowOfIndices is a Matrix of indices to select from data stored in a Matrix
// A new Matrix is created
Matrix Matrix::indexCols(Matrix &rowOfIndices)
{
    int sizeList;
    double *indices;
    
    assertDefined("indexRowVector first arg");
    rowOfIndices.assertDefined("indexRowVector second arg");
    rowOfIndices.assertRowVector("indexRowVector second arg");
    // more assertions

    sizeList = rowOfIndices.maxc;
    indices = rowOfIndices.m[0];
    
    Matrix out(maxr, sizeList, "");

    for (int i=0; i<sizeList; i++) {
        if (0 <= int(indices[i]) && int(indices[i]) < maxc) {
            for (int r=0; r<maxr; r++) {
                out.m[r][i] = m[r][int(indices[i])];
            }
        }
        else {
            printf("ERROR(indexCols): index %d is out of bounds for input row vector of size %d\n", int(indices[i]), maxc);

            exit(1);
        }
    }

    out.defined = true;

    return out;
}



// returns answer in arguments
void Matrix::argMax(int &rr, int &cc) const
{
    double max;

    assertDefined("argMax");

    max = m[0][0];
    rr = cc = 0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (m[r][c] > max) {
                max = m[r][c];
                rr = r;
                cc = c;
            }
        }
    }
}


// returns answer in arguments
void Matrix::argMin(int &rr, int &cc) const
{
    double min;

    assertDefined("argMin");

    min = m[0][0];
    rr = cc = 0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (m[r][c] < min) {
                min = m[r][c];
                rr = r;
                cc = c;
            }
        }
    }
}


// WARNING: allocates new matrix for answer
Matrix Matrix::argMinRow() const
{
    int cc;
    double min;
    assertDefined("argMinRow");

    Matrix out(maxr, 1);

    for (int r=0; r<maxr; r++) {
        min = m[r][0];
        cc = 0;
        for (int c=0; c<maxc; c++) {
            if (m[r][c] < min) {
                min = m[r][c];
                cc = c;
            }
        }
        out.m[r][0] = cc;
    }

    out.defined = true;

    return out;
}


// WARNING: allocates new matrix for answer
Matrix Matrix::argMaxRow() const
{
    int cc;
    double max;
    assertDefined("argMaxRow");

    Matrix out(maxr, 1);

    for (int r=0; r<maxr; r++) {
        max = m[r][0];
        cc = 0;
        for (int c=0; c<maxc; c++) {
            if (m[r][c] > max) {
                max = m[r][c];
                cc = c;
            }
        }
        out.m[r][0] = cc;
    }

    out.defined = true;

    return out;
}


// WARNING: allocates new matrix for answer
Matrix Matrix::minRow() const
{
    double min;
    assertDefined("minRow");

    Matrix out(maxr, 1);

    min = m[0][0];
    for (int r=0; r<maxr; r++) {
        min = m[r][0];
        for (int c=0; c<maxc; c++) {
            if (m[r][c] < min) {
                min = m[r][c];
            }
        }
        out.m[r][0] = min;
    }

    out.defined = true;

    return out;
}


// max of whole array
double Matrix::max() const
{
    double max;

    assertDefined("max");

    max = m[0][0];
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (m[r][c] > max) max = m[r][c];
        }
    }

    return max;
}


// min of whole array
double Matrix::min() const
{
    double min;

    assertDefined("min");

    min = m[0][0];
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (m[r][c] < min) min = m[r][c];
        }
    }

    return min;
}



// mean of whole array
double Matrix::mean() const
{
    double sum;

    assertDefined("mean");

    sum = 0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            sum += m[r][c];
        }
    }

    return sum/maxr/maxc;
}



// variance of whole array
// this does it brute force to avoid some problems with numerical stability
double Matrix::var() const
{
    double sum, sumd, mean;
    int size;

    assertDefined("var");

    size = maxr * maxc;

    sum = 0.0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            sum += m[r][c];
        }
    }
    mean = sum/size;

    sumd = 0.0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            sumd += (m[r][c] - mean) * (m[r][c] - mean);
        }
    }

    return sumd/size;
}


// standard deviation of whole array
double Matrix::stddev() const
{
    return sqrt(var());
}


// covariance of whole array
// this does it brute force to avoid some problems with numerical stability
double Matrix::cov(Matrix &other) const
{
    double sum, sumo, mean, meano;
    int size;

    assertDefined("lhs of cov");
    other.assertDefined("rhs of cov");
    assertOtherSizeMatch(other, "cov");

    size = maxr * maxc;

    sum = sumo = 0.0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            sum += m[r][c];
            sumo += other.m[r][c];
        }
    }
    mean = sum/size;
    meano = sumo/size;

    sum = 0.0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            sum += (m[r][c] - mean) * (other.m[r][c] - meano);
        }
    }

    return sum/size;
}



double Matrix::minCol(int c) const
{
    double min;

    assertDefined("min");
    assertColIndexOK(c, "minCol");

    min = m[0][c];
    for (int r=1; r<maxr; r++) {
        if (m[r][c] < min) min = m[r][c];
    }

    return min;
}


double Matrix::maxCol(int c) const
{
    double max;

    assertDefined("max");
    assertColIndexOK(c, "maxCol");

    max = m[0][c];
    for (int r=1; r<maxr; r++) {
        if (m[r][c] > max) max = m[r][c];
    }

    return max;
}


double Matrix::meanCol(int c) const
{
    double sum;

    assertDefined("sum");
    assertColIndexOK(c, "meanCol");

    sum = 0.0;
    for (int r=0; r<maxr; r++) sum += m[r][c];

    return sum/maxr;
}


// should do in a more numerically stable way
double Matrix::stddevCol(int c) const
{
    double sum, sumd, mean;
    int size;

    assertDefined("sum");
    assertColIndexOK(c, "stddevCol");

    size = maxr * maxc;

    sum = 0.0;
    for (int r=0; r<maxr; r++) {
        sum += m[r][c];
    }
    mean = sum/size;

    sumd = 0.0;
    for (int r=0; r<maxr; r++) {
        sumd += (m[r][c] - mean) * (m[r][c] - mean);
    }

    return sqrt(sumd/size);
}


// count number of items in column c equal to value
int Matrix::countEqCol(int c, double value) const
{
    int count;

    count = 0;
    for (int r=0; r<maxr; r++) {
        if (m[r][c]==value) count++;
    }

    return count;
}


// count number of items in column c not equal to value
int Matrix::countNeqCol(int c, double value) const
{
    int count;

    count = 0;
    for (int r=0; r<maxr; r++) {
        if (m[r][c]!=value) count++;
    }

    return count;
}



// normalize each row vector to have a length of 1
Matrix &Matrix::normalizeRowVectors()
{
    double sum;

    print("ME");

    // for each row
    for (int r=0; r<maxr; r++) {
        sum = 0;
        for (int c=0; c<maxc; c++) {
            sum += m[r][c]*m[r][c];
        }
        
        if (sum==0.0) {
            if (name.length()==0)
                printf("ERROR(normalizeRowVectors): row %d has zero length and so can't be normalized!\n",
                   r);
        else
            printf("ERROR(normalizeRowVectors): row %d in the matraix \"%s\" has zero length and so can't be normalized!\n",
                   r,
                   name.c_str());
            exit(1);
        }

        printf("S: %lf\n", sum);
        sum = 1.0/sqrt(sum);
        printf("S: %lf\n", sum);
        for (int c=0; c<maxc; c++) {
            m[r][c] *= sum;
        }
    }

    return *this;
}


// Normalizes the whole array in place so the values fall between 0 and 1.
void Matrix::normalize(double newmax)
{
    double min, max, span;

    assertDefined("normalize");

    min = max = m[0][0];
    for (int c=0; c<maxc; c++) {
        for (int r=0; r<maxr; r++) {
            if (m[r][c] < min) min = m[r][c];
            if (m[r][c] > max) max = m[r][c];
        }
    }
    span = max - min;

    // rescale
    if (span==0) {
        if (name.length()==0)
            printf("ERROR(normalize): all elements in the matrix are %lg and so the matrix can't be normalized!\n",
                   min);
        else
            printf("ERROR(normalize): all elements in the matrix \"%s\" are %lg and so the matrix can't be normalized!\n",
                   name.c_str(), min);
        exit(1);
    }

    span = newmax/span;
    for (int c=0; c<maxc; c++) {
        for (int r=0; r<maxr; r++) {
            m[r][c] = (m[r][c] - min)*span;
        }
    }
}


// returns new matrix for a matrix of min and max of each column
// normalizes within each column according to the min and max in that column
// so the range is now between 0 and 1 in each column
// NOTE: it will not rescale a column that is a constant!!
// WARNING: allocates a new matrix for answer AND alters matrix self
Matrix Matrix::normalizeCols()
{
    double min, max;

    assertDefined("normalizeCols");

    Matrix minMax(2, maxc, "minMax for " + name);

    for (int c=0; c<maxc; c++) {

        // find min and max
        min = max = m[0][c];
        for (int r=0; r<maxr; r++) {
            if (m[r][c] < min) min = m[r][c];
            if (m[r][c] > max) max = m[r][c];
        }

        // remember it
        minMax.m[0][c] = min;
        minMax.m[1][c] = max;

        // rescale column
        if (max!=min) {
            for (int r=0; r<maxr; r++) {
                m[r][c] = (m[r][c] - min)/(max - min);
            }
        }
    }
    minMax.defined = true;

    return minMax;
}


// normalizes within each column according to the min and max in that
// column supplied in the minMax matrix.  NOTE: This is used to scale
// two matrices the same way in the same columns.   This is useful
// for scaling training data and testing data.
Matrix &Matrix::normalizeCols(Matrix &minMax)
{
    double min, max;

    assertDefined("normalizeCols");

    for (int c=0; c<maxc; c++) {

        // recover min and max
        min = minMax.m[0][c];
        max = minMax.m[1][c];

        // rescale column
        if (min != max) {
            for (int r=0; r<maxr; r++) {
                m[r][c] = (m[r][c] - min)/(max - min);
            }
        }
    }

    return *this;
}

// undo the normalization done by normalizeCols using the minMax matrix
void Matrix::unnormalizeCols(Matrix &minMax)
{
    double min, max;

    assertDefined("normalizeCols");

    for (int c=0; c<maxc; c++) {
        // recover min and max
        min = minMax.m[0][c];
        max = minMax.m[1][c];

        // rescale column
        if (min != max) {
            for (int r=0; r<maxr; r++) {
                m[r][c] = m[r][c] * (max - min) + min;
            }
        }
    }
}


// this tests for true equality of all elements
bool Matrix::equal(const Matrix &other) const
{
    assertDefined("lhs of equal");
    other.assertDefined("rhs of equal");
    assertOtherSizeMatch(other, "equal");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (m[r][c] != other.m[r][c]) return false;
        }
    }

    return true;
}




// Return true if the differnce function is less than epsilon for
// all elements of the array.
bool Matrix::nearEqual(double epsilon, const Matrix &other) const
{
    assertDefined("lhs of nearEqual");
    other.assertDefined("rhs of nearEqual");
    assertOtherSizeMatch(other, "nearEqual");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (fabs(m[r][c]-other.m[r][c])>epsilon) return false;
        }
    }

    return true;
}


// the number of elements in self greater than the elements in other
int Matrix::countGreater(const Matrix &other) const
{
    int count;

    assertDefined("lhs of countGreater");
    other.assertDefined("rhs of countGreater");
    assertOtherSizeMatch(other, "countGreater");

    count = 0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (m[r][c] > other.m[r][c]) count++;
        }
    }

    return count;
}


// the number of elements in self greater than the elements in other
int Matrix::countGreater(const double value) const
{
    int count;

    assertDefined("countGreater");

    count = 0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (m[r][c] > value) count++;
        }
    }

    return count;
}


// sums up all the elements in the matrix
double Matrix::sum() const
{
    double sum;

    assertDefined("dist2");

    sum = 0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            sum += m[r][c];
        }
    }

    return sum;
}



// sum of the squares of the whole matrix
double Matrix::dist2() const
{
    double sum;

    assertDefined("dist2");

    sum = 0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            double tmp;

            tmp = m[r][c];
            sum += tmp * tmp;
        }
    }

    return sum;
}


// SQUARE of distance between two matrices
// this is an element by element operation and not like matrix multiply
double Matrix::dist2(const Matrix &other) const
{
    double sum;

    assertDefined("lhs of dist2");
    other.assertDefined("rhs of dist2");
    assertOtherSizeMatch(other, "dist2");

    sum = 0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            double tmp;

            tmp = m[r][c] - other.m[r][c];
            sum += tmp * tmp;
        }
    }

    return sum;
}



// SQUARE of distance between two matrices
// this is an element by element operation and not like matrix multiply
double Matrix::dist(const Matrix &other) const
{
    double sum;

    assertDefined("lhs of dist");
    other.assertDefined("rhs of dist");
    assertOtherSizeMatch(other, "dist");

    sum = 0;
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            double tmp;

            tmp = m[r][c] - other.m[r][c];
            sum += tmp * tmp;
        }
    }

    return sqrt(sum);
}


// magnitudes of the row vectors -> col vector
Matrix Matrix::distRow() const
{
    assertDefined("distRow");

    Matrix dist(maxr, 1);

    for (int r=0; r<maxr; r++) {
        double sum;

        sum = 0;
        for (int c=0; c<maxc; c++) {
            sum += m[r][c] * m[r][c];
        }
        dist.m[r][0] = sqrt(sum);
    }

    dist.defined = true;

    return dist;
}



// SQUARE of length of each row as a vector.  Same as x.Tdot(x).
Matrix Matrix::dist2Row() const
{
    assertDefined("dist2Row");

    Matrix dist(maxr, 1);

    for (int r=0; r<maxr; r++) {
        double sum;

        sum = 0;
        for (int c=0; c<maxc; c++) {
            sum += m[r][c] * m[r][c];
        }
        dist.m[r][0] = sum;
    }

    dist.defined = true;

    return dist;
}


// SQUARE of distance of row of this with col of other -> double
double Matrix::dist2(int r, int c, const Matrix &other) const
{
    double sum;

    assertDefined("lhs of dist2");
    other.assertDefined("rhs of dist2");
    assertOtherRhs(other, "dist2");

    sum = 0;
    for (int i=0; i<maxc; i++) {
        double tmp;

        tmp = m[r][i] - other.m[i][c];
        sum += tmp * tmp;
    }

    return sum;
}



// matrix multiply
// dot of row of this with col of other -> double
double Matrix::dot(int r, int c, const Matrix &other) const
{
    double sum;

    assertDefined("lhs of rowdot");
    other.assertDefined("rhs of rowdot");
    assertOtherRhs(other, "dot of row by col");

    sum = 0;
    for (int i=0; i<maxc; i++) {
        sum += m[r][i] * other.m[i][c];
    }

    return sum;
}



// +=
Matrix &Matrix::add(const Matrix &other)
{
    assertDefined("lhs of add");
    other.assertDefined("rhs of add");
    assertOtherSizeMatch(other, "add");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] += other.m[r][c];
        }
    }

    return *this;
}


// -=
Matrix &Matrix::sub(const Matrix &other)
{
    assertDefined("lhs of sub");
    other.assertDefined("rhs of sub");
    assertOtherSizeMatch(other, "sub");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] -= other.m[r][c];
        }
    }

    return *this;
}



// IMPORTANT: this is x - self   not   self - x
// can be used to negate
Matrix &Matrix::scalarPreSub(double x)
{
    assertDefined("scalarPreSub");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = x - m[r][c];
        }
    }

    return *this;
}


// IMPORTANT: this  self - x
Matrix &Matrix::scalarPostSub(double x)
{
    assertDefined("scalarPostSub");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = m[r][c] - x;
        }
    }

    return *this;
}



// add each column by a column vector
// the given default value.
Matrix &Matrix::addColVector(const Matrix &other)
{
    assertDefined("lhs of addColVector");
    other.assertDefined("rhs of addColVector");
    other.assertColVector("addColVector");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] += other.m[r][0];
        }
    }

    return *this;
}



// subtract each column by a column vector
// the given default value.
Matrix &Matrix::subColVector(const Matrix &other)
{
    assertDefined("lhs of subColVector");
    other.assertDefined("rhs of subColVector");
    other.assertColVector("subColVector");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] -= other.m[r][0];
        }
    }

    return *this;
}



// multiply each column by a column vector
// the given default value.
Matrix &Matrix::mulColVector(const Matrix &other)
{
    assertDefined("lhs of mulColVector");
    other.assertDefined("rhs of mulColVector");
    other.assertColVector("mulColVector");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] *= other.m[r][0];
        }
    }

    return *this;
}



// divide each column by a column vector
Matrix &Matrix::divColVector(const Matrix &other)
{
    assertDefined("lhs of divColVector");
    other.assertDefined("rhs of divColVector");
    other.assertColVector("divColVector");

    for (int r=0; r<maxr; r++) {
        if (other.m[r][0]==0.0) {
            if (name.length()==0) {
                printf("ERROR(divColVector): Trying to divide by element [%d, %d] of a matrix which is zero\n", r, 0);
            }
            else {
                printf("ERROR(divColVector): Trying to divide by element [%d, %d] of matrix named \"%s\" which is zero\n", r, 0, name.c_str());
            }
            exit(1);
        }

        for (int c=0; c<maxc; c++) {
            m[r][c] /= other.m[r][0];
        }
    }

    return *this;
}



// divide one matrix by another.  If denominator = 0 for an element use
// the given default value.
Matrix &Matrix::divRowVector(const Matrix &other)
{
    assertDefined("lhs of divRowVector");
    other.assertDefined("rhs of divRowVector");
    other.assertRowVector("divRowVector");

    for (int c=0; c<maxc; c++) {
        if (other.m[0][c]==0.0) {
            if (other.name.length()==0) {
                printf("ERROR(divRowVector): Trying to divide by element [%d, %d] of a matrix which is zero\n", 0, c);
            }
            else {
                printf("ERROR(divRowVector): Trying to divide by element [%d, %d] of matrix named \"%s\" which is zero\n", 0, c, other.name.c_str());
            }
            exit(1);
        }

        for (int r=0; r<maxr; r++) {
            m[r][c] /= other.m[0][c];
        }
    }

    return *this;
}



// multiply each row in self by row vector matrix in other
Matrix &Matrix::mulRowVector(const Matrix &other)
{
    assertDefined("mulRowVector");
    assertColsEqual(other, "mulRowVector");
    other.assertRowVector("mulRowVector");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] *= other.m[0][c];
        }
    }

    return *this;
}


// add a row vector matrix in other to each row of self
Matrix &Matrix::addRowVector(const Matrix &other)
{
    assertDefined("addRowVector");
    assertColsEqual(other, "addRowVector");
    other.assertRowVector("addRowVector");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] += other.m[0][c];
        }
    }

    return *this;
}

// add a row vector matrix in other to the given row of self
Matrix &Matrix::addRowVector(int r, const Matrix &other)
{
    assertDefined("addRowVector");
    assertColsEqual(other, "addRowVector");
    other.assertRowVector("addRowVector");

    for (int c=0; c<maxc; c++) {
        m[r][c] += other.m[0][c];
    }

    return *this;
}



// subtract row matrix to each row of self
Matrix &Matrix::subRowVector(const Matrix &other)
{
    assertDefined("subRowVector");
    assertColsEqual(other, "subRowVector");
    other.assertRowVector("subRowVector");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] -= other.m[0][c];
        }
    }

    return *this;
}



// element by element absolute value in place
Matrix &Matrix::Matrix::abs()
{
    assertDefined("abs");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = fabs(m[r][c]);
        }
    }

    return *this;
}




// element by element multiply
Matrix &Matrix::mul(const Matrix &other)
{
    assertDefined("lhs of mul");
    other.assertDefined("rhs of mul");
    assertOtherSizeMatch(other, "mul");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] *= other.m[r][c];
        }
    }

    return *this;
}


// element by element divide
Matrix &Matrix::div(const Matrix &other)
{
    assertDefined("lhs of div");
    other.assertDefined("rhs of div");
    assertOtherSizeMatch(other, "div");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (other.m[r][c]==0.0) {
                if (name.length()==0) {
                    printf("ERROR(div): Trying to divide by element [%d, %d] of a matrix which is zero\n", r, c);
                }
                else {
                    printf("ERROR(div): Trying to divide by element [%d, %d] of matrix named \"%s\" which is zero\n", r, c, name.c_str());
                }
                exit(1);
            }
            m[r][c] /= other.m[r][c];
        }
    }

    return *this;
}


// increment the values in a given row by 1
Matrix &Matrix::rowInc(int r)
{
    for (int c=0; c<maxc; c++) {
        m[r][c]++;
    }

    return *this;
}



// swap two matrices
Matrix &Matrix::swap(Matrix &other)
{
    assertDefined("lhs of swap");
    other.assertDefined("rhs of swap");
    assertOtherSizeMatch(other, "swap");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            double tmp;

            tmp = m[r][c];
            m[r][c] = other.m[r][c];
            other.m[r][c] = tmp;
        }
    }

    return *this;
}



// pick from SELF those rows where LIST has an element equal to
// MATCH. can be used to do a complex selection from an array based on
// classes for each row stored in another array. Essentially find all
// elements in the list argument that are equal to match and select
// the corresponding row in self. The number of matches is equal to
// the number of rows in the returned Matrix.
// WARNING: allocates new matrix for answer.
// WARNING: number of rows could be zero.  User should check.
// Also see submatrix version.
Matrix Matrix::pickRows(int match, const Matrix &list, int matchCol)
{
    int num;

    assertDefined("lhs of pickRows");
    list.assertDefined("rhs of pickRows");
    list.assertColVector("rhs of pickRows");
    assertRowsEqual(list, "pickRows");

    num = 0;
    for (int r=0; r<maxr; r++) if (list.m[r][0]==match) num++;

    // WARNING: number of rows could be zero!!

    Matrix out(num, maxc);

    { int rr=0;
        for (int r=0; r<maxr; r++) {
            if (list.m[r][matchCol]==match) {
                for (int c=0; c<maxc; c++) {
                    out.m[rr][c] = m[r][c];
                }
                rr++;
                if (rr>=num) break;
            }
        }
    }

    out.defined = true;

    return out;
}



// join (append) the other matrix on the right side of self returning a new matrix
// NOTE: if self is undefined then it simply copies other.
// WARNING: allocates new matrix for answer
Matrix Matrix::joinRight(Matrix &other)
{
    other.assertDefined("joinRight");

    if (!isDefined()) {
        Matrix out(other);
        return out;
    }
    else {
        Matrix out(maxr, maxc + other.maxc);
    
        assertRowsEqual(other, "joinRight");

        // copy left matrix
        for (int r=0; r<maxr; r++) {
            for (int c=0; c<maxc; c++) {
                out.m[r][c] = m[r][c];
            }
        }

        // copy right matrix
        for (int r=0; r<other.maxr; r++) {
            for (int c=0; c<other.maxc; c++) {
                out.m[r][c+maxc] = other.m[r][c];
            }
        }

        out.defined = true;
        return out;
    }
}


// join (append) the other matrix on the right side of self returning a new matrix
// NOTE: if self is undefined then it simply copies other.
// WARNING: allocates new matrix for answer
Matrix Matrix::joinBottom(Matrix &other)
{
    other.assertDefined("joinBottom");

    if (!isDefined()) {
        Matrix out(other);
        return out;
    }
    else {
        Matrix out(maxr + other.maxr, maxc);
    
        assertColsEqual(other, "joinBottom");

        // copy top matrix
        for (int r=0; r<maxr; r++) {
            for (int c=0; c<maxc; c++) {
                out.m[r][c] = m[r][c];
            }
        }

        // copy bottom matrix
        for (int r=0; r<other.maxr; r++) {
            for (int c=0; c<other.maxc; c++) {
                out.m[r+maxr][c] = other.m[r][c];
            }
        }

        out.defined = true;

        return out;
    }
}


// dot or inner product or classic matrix multiply
// WARNING: allocates new matrix for answer
Matrix Matrix::dot(const Matrix &other)
{
    assertDefined("lhs of dot");
    other.assertDefined("rhs of dot");
    assertOtherRhs(other, "dot");

    Matrix out(maxr, other.maxc);
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<other.maxc; c++) {
            double sum;

            sum = 0;
            for (int i=0; i<maxc; i++) {
                sum += m[r][i] * other.m[i][c];
            }
            out.m[r][c] = sum;
        }
    }

    out.defined = true;

    return out;
}



// dot or inner product or classic matrix multiply BUT
// the SECOND argument is transposed!
// WARNING: allocates new matrix for answer!
Matrix Matrix::dotT(const Matrix &other)
{
    assertDefined("lhs of dotT");
    other.assertDefined("rhs of dotT");
    assertColsEqual(other, "dotT");

    Matrix out(maxr, other.maxr);

    for (int r=0; r<maxr; r++) {               // use columns from first
        for (int c=0; c<other.maxr; c++) {
            double sum;

            sum = 0;
            for (int i=0; i<maxc; i++) {       // sum over columns
                sum += m[r][i] * other.m[c][i];  // go down the transpose
            }
            out.m[r][c] = sum;
        }
    }

    out.defined = true;

    return out;
}



// dot or inner product or classic matrix multiply BUT
// the FIRST argument is transposed!
// WARNING: allocates new matrix for answer
Matrix Matrix::Tdot(const Matrix &other)
{
    assertDefined("lhs of Tdot");
    other.assertDefined("rhs of Tdot");
    assertRowsEqual(other, "Tdot");

    Matrix out(maxc, other.maxc);        // use columns from first
    for (int r=0; r<maxc; r++) {               // use columns from first
        for (int c=0; c<other.maxc; c++) {
            double sum;

            sum = 0;
            for (int i=0; i<maxr; i++) {       // sum over rows
                sum += m[i][r] * other.m[i][c];  // go down the transpose
            }
            out.m[r][c] = sum;
        }
    }

    out.defined = true;

    return out;
}



// this is the dot product of two row vectors which is the same as
// a.dotT(b).get(0, 0)
double Matrix::dotRowVector(const Matrix &other)
{
    assertDefined("lhs of dotRowVector");
    assertRowVector("dotRowVector");
    other.assertDefined("rhs of dotRowVector");
    other.assertRowVector("dotRowVector");
    assertColsEqual(other, "dotRowVector");
    
    double sum;

    sum = 0;
    for (int c=0; c<maxc; c++) {         // sum over columns
        sum += m[0][c] * other.m[0][c];
    }

    return sum;
}


// this is the dot product of two col vectors which is the same as
// a.Tdot(b).get(0, 0)
double Matrix::dotColVector(const Matrix &other)
{
    assertDefined("lhs of dotColVector");
    assertColVector("dotColVector");
    other.assertDefined("rhs of dotColVector");
    other.assertColVector("dotColVector");
    assertColsEqual(other, "dotColVector");
    
    double sum;

    sum = 0;
    for (int r=0; r<maxr; r++) {         // sum over rows
        sum += m[r][0] * other.m[r][0];
    }

    return sum;
}



// This computes the mean of every column and puts it into a row vector
// which is the same as computing the mean of the row vectors in a matrix
// WARNING: allocates new matrix for answer
Matrix Matrix::meanRowVectors()
{
    Matrix mean(1, maxc);

    for (int c=0; c<maxc; c++) {
        double sum;

        sum = 0;
        for (int r=0; r<maxr; r++) {
            sum += m[r][c];
        }
        mean.m[0][c] = sum/maxr;
    }

    mean.defined = true;

    return mean;
}


// This computes the mean of every column and puts it into a row vector
// WARNING: allocates new matrix for answer
Matrix Matrix::stddevRowVectors() {
    assertDefined("stddevRowVectors");

    Matrix stddev(1, maxc);
    for (int c=0; c<maxc; c++) {
        double sum, sumd, mean;
        int size;

        size = maxr * maxc;

        sum = 0.0;
        for (int r=0; r<maxr; r++) {
            sum += m[r][c];
        }
        mean = sum/size;

        sumd = 0.0;
        for (int r=0; r<maxr; r++) {
            sumd += (m[r][c] - mean) * (m[r][c] - mean);
        }

        stddev.m[0][c] = sqrt(sumd/size);
    }

    stddev.defined = true;

    return stddev;
}



// covariance matrix with each column being a random variable and the matrix
// being the covariance between columns.
// WARNING: This is NOT the unbiased covariance in which
// you divide by (n - 1)!  In this routine we divide by n.
//
// WARNING: allocates new matrix for answer
Matrix Matrix::covMatrix()
{
    assertDefined("covMatrix");

    double inv;
    double *mean;

    // get the mean of each column
    mean = new double [maxc];
    for (int c=0; c<maxc; c++) {
        double sum;

        sum = 0;
        for (int r=0; r<maxr; r++) {
            sum += m[r][c];
        }
        mean[c] = sum/maxr;
    }

    Matrix out(maxc, maxc);
    inv = 1.0/maxr;
    for (int r=0; r<maxc; r++) {
        for (int c=r; c<maxc; c++) {  // note: starts at r
            double sum;

            sum = 0;
            for (int i=0; i<maxr; i++) {       // sum over rows
                sum += (m[i][r]-mean[r]) * (m[i][c]-mean[c]);  // go down the transpose
            }
            out.m[r][c] = out.m[c][r] = sum * inv;
        }
    }

    out.defined = true;
    delete [] mean;

    return out;
}


// returns the covariance between two matrices.  Each column
// is treated as a random variable so does a column by column
// comparision.
// WARNING: This is NOT the unbiased covariance in which
// you divide by (n - 1)!  In this routine we divide by n.
// This matrix is not necessarily symmetric!
//
// WARNING: allocates new matrix for answer
Matrix Matrix::covMatrix(Matrix &other)
{
    assertDefined("covMatrix");
    assertRowsEqual(other, "covMatrix");

    double inv;
    double *mean, *meano;

    // get mean in each column
    mean = new double [maxc];
    for (int c=0; c<maxc; c++) {
        double sum;

        sum = 0;
        for (int r=0; r<maxr; r++) {
            sum += m[r][c];
        }
        mean[c] = sum/maxr;
    }

    // get mean in each column of other
    meano = new double [other.maxc];
    for (int c=0; c<other.maxc; c++) {
        double sum;

        sum = 0;
        for (int r=0; r<other.maxr; r++) {
            sum += other.m[r][c];
        }
        meano[c] = sum/maxr;
    }

    Matrix out(maxc, other.maxc);

    inv = 1.0/maxr;
    for (int r=0; r<maxc; r++) {
        for (int c=0; c<other.maxc; c++) {
            double sum;

            sum = 0;
            for (int i=0; i<maxr; i++) {       // sum over rows
                sum += (m[i][r]-mean[r]) * (other.m[i][c]-meano[c]);  // go down the transpose
            }
//            out.m[r][c] = out.m[c][r] = sum * inv;
            out.m[r][c] = sum * inv;
        }
    }

    out.defined = true;
    delete [] mean;
    delete [] meano;

    return out;
}



Matrix &Matrix::zero()
{
    constant(0.0);

    return *this;
}


Matrix &Matrix::identity()
{
    assertSquare("identity");

    constant(0.0);
    constantDiagonal(1.0);

    return *this;
}


// scalar multiply
Matrix &Matrix::scalarMul(double x)
{
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] *= x;
        }
    }

    return *this;
}


// // scalar divide
Matrix &Matrix::scalarDiv(double x)
{
    if (x==0) {
        if (name.length()==0) {
            printf("ERROR(scalarDiv): Trying to divide matrix by zero\n");
        }
        else {
            printf("ERROR(scalarDiv): Trying to divide matrix \"%s\" by zero\n", name.c_str());
        }
        exit(1);
    }

    return scalarMul(1.0/x);   // mul is faster than div
}



// scalar add
Matrix &Matrix::scalarAdd(double x)
{
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] += x;
        }
    }

    return *this;
}



// performs the function on each row of the matrix producing a column vector of answers
// WARNING: allocates new matrix for answer
Matrix Matrix::mapRow(double (*f)(int size, double *x))
{
    assertDefined("mapRow");

    Matrix out(maxr, 1);

    for (int r=0; r<out.maxr; r++) {
        out.m[r][0] = f(maxc, m[r]);
    }
    out.defined = true;

    return out;
}


// performs the function on each col of the matrix producing a row vector of answers
// WARNING: allocates new matrix for answer
Matrix Matrix::mapCol(double (*f)(int size, double *x))
{
    assertDefined("mapCol");

    Matrix out(1, maxc);
    double *tmp = new double [maxr];

    for (int c=0; c<out.maxc; c++) {
        for (int r=0; r<out.maxr; r++) tmp[r] = m[r][c];
        out.m[0][c] = f(maxr, tmp);
    }
    out.defined = true;

    delete [] tmp;

    return out;
}



// performs the function over the cartesian product over the rows of
// the two matrices.  WARNING: allocates new matrix for answer
Matrix Matrix::cartesianRow(double (*f)(int size, double *x, double *y), const Matrix &other)
{
    assertDefined("cartesianRow");
    other.assertDefined("cartesianRow");
    assertColsEqual(other, "cartesianRow");

    Matrix out(maxr, other.maxr);

    for (int r=0; r<out.maxr; r++) {
        for (int c=0; c<out.maxc; c++) {
            out.m[r][c] = f(maxc, m[r], other.m[c]);
        }
    }
    out.defined = true;

    return out;
}


// Sample a column as if it is a time series. sampling is done
// in column specified with given stride and numsteps numsteps.
// The NUMBER OF COLUMNS of the output matrix is numsteps+1.
// NOTE: this is NOT a random sample
// A new matrix is allocated by this routine.
// Example:   x...x...x...y   were there are numstep x's and y is the target.
Matrix Matrix::seriesSampleCol(int col, int numsteps, int stride)
{
    assertColIndexOK(col, "seriesSampleCol");

    if (maxr - numsteps*stride<1) {
        if (name.length()==0)
            printf("ERROR(seriesSample): Trying to sample matrix with numsteps: %d and stride: %d but matrix is only %d X %d\n", numsteps, stride, maxr, maxc);
        else
            printf("ERROR(seriesSample): Trying to sample matrix \"%s\" with numsteps: %d and stride: %d but matrix is only %d X %d\n", name.c_str(), numsteps, stride, maxr, maxc);
        exit(1);
    }

    Matrix out(maxr - numsteps*stride, numsteps+1, "out");

    for (int r=0; r+numsteps*stride<maxr; r++) {
        for (int c=0; c<=numsteps; c++) {
            out.m[r][c] = m[r+c*stride][col];
        }
    }
    out.defined = true;

    return out;
}



// apply a function to every element
// WARNING: overwrites self
Matrix &Matrix::map(double (*f)(double x))
{
    assertDefined("map");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = f(m[r][c]);
        }
    }

    return *this;
}


// apply a function to every element in a column
// WARNING: overwrites self
Matrix &Matrix::mapCol(int c, double (*f)(double x))
{
    assertDefined("mapCol");
    assertColIndexOK(c, "mapCol");

    for (int r=0; r<maxr; r++) {
        m[r][c] = f(m[r][c]);
    }

    return *this;
}


// apply a function to each row feeding the whole row to a function to let
// the function do what it wants to the row
// WARNING: may overwrite self
Matrix &Matrix::mapEachRowSelf(void (*f)(int size, double *x))
{
    assertDefined("mapEachRowSelf");

    for (int r=0; r<maxr; r++) {
        f(maxc, m[r]);
    }

    return *this;
}


// apply a function to each col feeding the whole col to a function to let
// the function do what it wants to the col
Matrix &Matrix::mapEachCol(void (*f)(int size, double *x))
{
    assertDefined("mapEachCol");

    double *arg = new double [maxr];

    for (int c=0; c<maxc; c++) {
        for (int r=0; r<maxr; r++) {
            arg[r] = m[r][c];
        }
        f(maxr, arg);
        // zzz could copy back if mapeachcolself version
    }

    return *this;
}


// apply a function to each row feeding the whole row to a function to let
// the function do what it wants to the row
// WARNING: may overwrite self
Matrix &Matrix::mapEachRowIndexSelf(void (*f)(int size, int r, double *x))
{
    assertDefined("mapEachRowIndexSelf");

    for (int r=0; r<maxr; r++) {
        f(maxc, r, m[r]);
    }

    return *this;
}


// apply a function to every element and its index
// NOTE: it does not check if the array is undefined or not
// so the function is free to use only the index pair.
// WARNING: overwrites self
Matrix &Matrix::mapIndex(double (*f)(int r, int c, double x))
{
    assertDefined("mapIndex");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = f(r, c, m[r][c]);
        }
    }

    defined = true;  // this may not be true if function uses undefined values

    return *this;
}



// initializes the matrix to a constant
Matrix &Matrix::constant(double x)
{
    assertUsableSize("constant");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = x;
        }
    }

    defined = true;

    return *this;
}



// initializes a column in a matrix to a constant
Matrix &Matrix::constantCol(int c, double x)
{
    assertColIndexOK(c, "constantCol");

    for (int r=0; r<maxr; r++) {
        m[r][c] = x;
    }

    defined = true;

    return *this;
}


// INITIALIZES A COLUMN in a matrix to a constant
// WARNING: even though this only sets one column it marks the matrix as defined!
//  This is because this is often used to init a matrix.
Matrix &Matrix::constantColRange(int c, double start, double step)
{
    assertColIndexOK(c, "constantColRange");

    for (int r=0; r<maxr; r++) {
        m[r][c] = start;
        start += step;
    }

    defined = true;

    return *this;
}



// INITIALIZES A ROW in a matrix to a constant
// WARNING: even though this only sets one row it marks the matrix as defined!
//  This is because this is often used to init a matrix.
Matrix &Matrix::constantRowRange(int r, double start, double step)
{
    assertRowIndexOK(r, "constantRowRange");

    for (int c=0; c<maxc; c++) {
        m[r][c] = start;
        start += step;
    }

    defined = true;

    return *this;
}





// initializes the diagonal of an existing matrix to a constant
Matrix &Matrix::constantDiagonal(double x)
{
    int len;

    len = maxr;
    if (maxc<maxr) len = maxc;

    for (int r=0; r<len; r++) {
            m[r][r] = x;
    }

//TOOBOLD    defined = true;   // makes bold assumption that the rest of matrix set

    return *this;
}


// init m[r][c] = A*r + B*c + C;
Matrix &Matrix::initLinear(double A, double B, double C) 
{
    assertAllocated("initLinear");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = A*r + B*c + C;
        }
    }

    defined = true;

    return *this;
}





 // fill with random doubles in the given range: [min, max)
Matrix &Matrix::rand(double min, double max)
{
    assertRandInitialized("rand");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = randUnit()*(max-min) + min;
        }
    }

    defined = true;

    return *this;
}


// fill the given column with random doubles in the given range: [min, max)
// does not set the state of undefined
Matrix &Matrix::randCol(int c, double min, double max)
{
    assertRandInitialized("randCol");

    for (int r=0; r<maxr; r++) {
        m[r][c] = randUnit()*(max-min) + min;
    }

    return *this;
}


// random reals in normal distribution
Matrix &Matrix::randNorm(double mean, double stddev)
{
    assertRandInitialized("randNorm");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = ::randNorm(stddev) + mean;
        }
    }

    defined = true;

    return *this;
}


// fill with random integers in the given range: [min, max)
Matrix &Matrix::rand(int min, int max)
{
    assertRandInitialized("rand");

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = randMod(max-min) + min;
        }
    }

    defined = true;

    return *this;
}


// extracts a random sample WITH REPLACEMENT of rows FROM DATA and
// PUTS IT INTO SELF. The size of the sample is the number of rows
// already defined to be in self. The values are copied from data to
// self. This is designed for repeated sampling. The routine returns
// the sample matrix. NOTE: the number of rows extracted are the number of
// rows defined in out.
Matrix &Matrix::sample(Matrix &data)
{
    data.assertDefined("sample");
    assertColsEqual(data, "sample");
    assertRandInitialized("sample");

    for (int r=0; r<maxr; r++) {
        int otherRow;

        otherRow = randMod(data.maxr);
        for (int c=0; c<maxc; c++) {
            m[r][c] = data.m[otherRow][c];
        }
    }

    defined = true;

    return *this;
}



// extracts a random sample WITHOUT REPLACEMENT OF ROWS from Matrix data and
// PUTS IT INTO SELF. The size of the sample is the number of rows in
// self! The values are copied from data Matrix to self Matrix. This is designed for
// repeated sampling. The routine returns the sample matrix.

// WARNING: the sample by construction is semi-sorted. Don't assume
// randomly ordered. This is because it uses Floyd's algorithm for
// sampling without replacement.
Matrix &Matrix::sampleWithoutRows(Matrix &data)
{
    data.assertDefined("sampleWithoutRows");
    assertColsEqual(data, "sampleWithoutRows");
    assertRandInitialized("sampleWithoutRows");

    int s, n, next;

    s = maxr;  // m
    n = data.maxr;

    // floyd's algorithm
    int list[s];

    next = 0;
    for (int j=n-s; j<n; j++) {
        int t;
        bool found;
        
        t = randMod(j+1);

        found = false;
        for (int i=0; i<next; i++) {
            if (t == list[i]) {
                found = true;
                break;
            }
        }
        if (found) list[next++] = j;
        else list[next++] = t;
    }

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = data.m[list[r]][c];
        }
    }
    defined = true;

    return *this;
}


// extracts a random sample WITHOUT REPLACEMENT of cols from RowVector data and
// PUTS IT INTO SELF. The size of the sample is the number of cols in
// self! The values are copied from data Matrix to self Matrix. This is designed for
// repeated sampling. The routine returns the sample matrix. WARNING:
// the sample by construction is semi-sorted. Don't assume randomly
// ordered. This is becaue it uses Floyd's algorithm for sampling
// without replacement.
Matrix &Matrix::sampleWithoutCols(Matrix &data)
{
    data.assertDefined("sampleWithoutCols");
    assertRowsEqual(data, "sampleWithoutCols");
    assertRandInitialized("sampleWithoutCols");

    int s, n, next;

    s = maxc;
    n = data.maxc;

    // floyd's algorithm
    int list[s];    // array is size of self (the requested number of samples)

    next = 0;
    for (int j=n-s; j<n; j++) {
        int t;
        bool found;
        
        t = randMod(j+1);

        found = false;
        for (int i=0; i<next; i++) {
            if (t == list[i]) {
                found = true;
                break;
            }
        }
        if (found) list[next++] = j;
        else list[next++] = t;
    }

//    for (int r=0; r<maxr; r++) {
//        for (int c=0; c<maxc; c++) {
//            m[r][c] = data.m[list[r]][c];
//        }
//    }

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            m[r][c] = data.m[r][list[c]];
        }
    }
    defined = true;

    return *this;
}


// randomly shuffles the rows of the matrix in place
Matrix &Matrix::shuffle()
{
    int rr;

    assertDefined("shuffle");
    assertRandInitialized("shuffle");

    if (!isInitRand()) {
        if (name.length()==0)
            printf("ERROR(shuffle): requires initialization first with a call to initRand()\n");
            
        else
            printf("ERROR(shuffle): on matrix \"%s\" requires initialization first with a call to initRand()\n", name.c_str());
        exit(1);
    }

    // loop through selecting the element to place at location r
    for (int r=0; r<maxr-1; r++) {
	double *tmp;

	rr=r+randMod(maxr-r);
	tmp = m[r];
	m[r] = m[rr];
	m[rr] = tmp;
    }

    return *this;
}


// WARNING: allocates new matrix for answer
Matrix Matrix::transpose() const
{
    assertDefined("transpose");

    Matrix out(maxc, maxr);

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            out.m[c][r] = m[r][c];
        }
    }
    out.defined = true;

    return out;
}




// transposes in place, but will reallocate and copy if a nonsquare matrix!
// WARNING: overwrites self
Matrix &Matrix::transposeSelf()
{
    assertDefined("transposeSelf");

    // handle square matrix in place
    if (maxr == maxc) {
        for (int r=0; r<maxr; r++) {
            for (int c=r+1; c<maxc; c++) {
                double tmp;
                tmp = m[r][c]; m[r][c] = m[c][r]; m[c][r] = tmp;
            }
        }
    }
    // handle non-square matrix with reallocation
    else {
        double **newm;

        newm = new double * [maxc];
        for (int i=0; i<maxc; i++) newm[i] = new double [maxr];

        for (int r=0; r<maxr; r++) {
            for (int c=0; c<maxc; c++) {
                newm[c][r] = m[r][c];
            }
        }

        deallocate();  // deallocate AFTER copying

        { int tmp; tmp = maxr; maxr = maxc; maxc = tmp; }
        m = newm;
        defined = true;
    }

    return *this;
}



// LU decomposition IN PLACE
// Uses simple Dolittle Algorithm
// Returns the permuation of the rows
int *Matrix::LU()
{
    int *perm;

    assertDefined("LU decomposition");

    perm = new int(maxr);
    for (int r=0; r<maxr; r++) perm[r] = r;

    for (int r=0; r<maxr; r++) {
        for (int rr=r+1; rr<maxr; rr++) {
            double l;

            if (m[r][r]==0) {
                for (int j=r+1; j<maxr; j++) {
                    if (m[j][r]!=0) {
                        int tmp;
                        swapRows(r, j);
                        tmp = perm[r]; perm[r] = perm[j]; perm[j] = tmp;
                        break;
                    }
                }
            }

            l = m[rr][r]/m[r][r];
            for (int c=r; c<maxc; c++) {
                m[rr][c] -= l * m[r][c];
            }
            m[rr][r] = l;
        }
    }

    for (int r=0; r<maxr; r++) printf("%d\n", perm[r]);

    return perm;
}


// solve Ax = B where A is this matrix object and B is a matrix in
// which each COLUMN is a vector to solve for.
// output: this matrix is replaced by its matrix inverse, and argument
// matrix B is replaced by the corresponding set of solution
// vectors.
Matrix &Matrix::solve(Matrix &B)
{
    assertSquare("solve");

    if (!gaussj(m, maxc, B.m, B.maxc)) {
        if (name.length()==0)
            printf("ERROR(solve): matrix is singular\n");
        else
            printf("ERROR(solve): matrix \"%s\" is singular\n", name.c_str());
        exit(1);
    }

    return B;
}


// replaces a matrix by it's inverse and also returns a reference to itself
Matrix &Matrix::inverse()
{
    assertSquare("inverse");

    if (!gaussj(m, maxc, NULL, 0)) {
        if (name.length()==0)
            printf("ERROR(solve): matrix is singular\n");
        else
            printf("ERROR(solve): matrix \"%s\" is singular\n", name.c_str());
        exit(1);
    }

    return *this;
}


// just print the size and name of the matrix
void Matrix::printSize(std::string msg, bool size) const
{
    if (msg.length()) {
        printf("%s ", msg.c_str());
    }

    if (size) {
        if (name.length()) {
            printf("(size of %s: %d X %d)\n", name.c_str(), maxr, maxc);
        }
        else {
            printf("(size: %d X %d)\n", maxr, maxc);
        }
    }

    fflush(stdout);
}


// does a formated printing of the matrix
//  msg is any leading string you want printed first
//  fmt is the format for the numbers in the matrix (fmt="" sets to default)
//  if size is true it will print the name and size of the matrix.  The default
//  is to print the matrix WITHOUT the name and size!
const Matrix &Matrix::printFmt(std::string msg, std::string fmt, bool size)
{
    assertDefined("printFmt");

    printSize(msg, size);

    if (fmt=="") fmt=Matrix::realFormat;

    const char *fmtStr = fmt.c_str();  // convert once
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            printf(fmtStr, m[r][c]);
        }
        printf("\n");
    }

    fflush(stdout);

    return *this;
}


const Matrix &Matrix::print(std::string msg, bool size, bool newline) const
{
    assertDefined("print");

    printSize(msg, size);

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            printf(realFormat, m[r][c]);
        }
        if (newline) printf("\n");
    }

    fflush(stdout);

    return *this;
}


// print the whole matrix including it's name and size as integers
void Matrix::printInt(std::string msg, bool size) const
{
    assertDefined("printInt");

    printSize(msg, size);

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (m[r][c] != int(m[r][c])) {
                printf("ERROR(printInt): Trying to print an integer matrix but element at position %d, %d is %10.5lg which is not an integer\n", r, c, m[r][c]);
                exit(1);
            }
            printf(shortIntFormat, int(m[r][c]));
        }
        printf("\n");
    }

    fflush(stdout);
}


const Matrix &Matrix::printMathematica(std::string msg) const
{
    assertDefined("printMathematica");

    if (msg.length()) {
        printf("\"%s\"\n", msg.c_str());
    }
    
    printf("{\n");
    for (int r=0; r<maxr; r++) {
        printf("{");

        printf("%g", m[r][0]);
        for (int c=1; c<maxc; c++) {
            printf(",%g", m[r][c]);
        }

        if (r<maxr-1) printf("},\n");
        else printf("}\n");
        fflush(stdout);
    }

    printf("}\n");
    fflush(stdout);

    return *this;
}


// print the whole matrix including it's name and size as one
// character per number
static int strlen(char *s)
{
    int i;

    i = 0;
    while (*s++) i++;

    return i;
}


void Matrix::printChar(char *code, std::string msg, bool size) const
{
    int len;

    assertDefined("printChar");

    printSize(msg, size);

    len = strlen(code);
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            int n;
            char ch;

            n = int(m[r][c]);
            if (n<0) ch = code[0];
            else if (n>=len) ch = code[len-1];
            else ch = code[n];

            printf("%c", ch);
        }
        printf("\n");
    }

    fflush(stdout);
}



void Matrix::printNZ(double epsilon, std::string msg, bool size) const
{
    assertDefined("printNZ");

    printSize(msg, size);

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (fabs(m[r][c]) < epsilon) {
                printf(intFormat, 0);
            }
            else {
                printf(realFormat, m[r][c]);
            }
        }
        printf("\n");
    }

    fflush(stdout);
}




// print the whole matrix including it's name and size
// included in this call is a list of row labels which are indexed by column 0.
// zzz WARNING: no attempt to make sure column 0 indices are in range.
void Matrix::printLabeledRow(const SymbolNumMap *labels, std::string msg, bool size) const
{
    assertDefined("printLabeledRow");

    printSize(msg, size);

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            if (c==0) {
                    printf("%s ", labels->getStr(m[r][0]).c_str());
            }
            else {
                printf("%lg ", m[r][c]);
            }
        }
        printf("\n");
    }

    fflush(stdout);
}


// print matrix containing only strings
void Matrix::printStrings(const SymbolNumMap *symbols, std::string msg, bool size) const
{
    assertDefined("printStrings");

    printSize(msg, size);

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            printf("%s ", symbols->getStr(m[r][c]).c_str());
        }
        printf("\n");
    }

    fflush(stdout);
}



// write out just the matrix data in a form that can be read back in
void Matrix::write() const
{
    assertDefined("write");

    printf("%d %d\n", maxr, maxc);
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {
            printf("%.15lg ", m[r][c]);
        }
        printf("\n");
    }
}




// this will print a row with a terminal blank but no terminal newline
void Matrix::writeLine(int r) const
{
    assertDefined("writeLine");
    assertIndexOK(r, 0, "writeLine");

    for (int c=0; c<maxc; c++) {
        printf(realFormat, m[r][c]);
    }
}



// this will print a row with a terminal blank but no terminal newline
void Matrix::writeLineStrings(const SymbolNumMap *symbols, int r) const
{
    assertDefined("writeLineStrings");
    assertIndexOK(r, 0, "writeLineStrings");

    for (int c=0; c<maxc; c++) {
        printf("%s ", symbols->getStr(m[r][c]).c_str());
    }
}



// read a matrix in.
// first two numbers are the number of rows and columns.
// then the matrix values by row.
// will deallocate old array if different size.
void Matrix::read()
{
    readAux(NUM, false, NULL);
}


// just like read() but matrix will be transposed
// WARNING: does not work with labeled matrices at the moment
void Matrix::readT()
{
    readAux(NUM, true, NULL);
}


// read a possibly labeled matrix that is in the format
//  row col
//  label 1 2 3 .. col
//  label 1 2 3 .. col
//  label 1 2 3 .. col
//  .
//  .
//  .
//
// where label is a continuous collection of non-whitespace
// characters. (Uses %s in scanf to read the label).
// WARNING: returns an array of char * strings.  The caller
// is responsible for deallocating the array of strings.
SymbolNumMap *Matrix::readLabeledRow(SymbolNumMap *syms)
{
    return readAux(LABELEDROW, false, syms);
}


// read in a matrix which contains all strings.
// default for syms is NULL
SymbolNumMap *Matrix::readStrings(SymbolNumMap *syms)
{
    return readAux(STRINGS, false, syms);
}


// if char s is in C style string t then true else false.
// NOTE: this is NOT strchr.
static inline bool strisin(const char s, const char *t)
{
    while (*t) {
        if (*t++==s) return true;
    }
    return false;
}



// read in the matrix assuming the size of the matrix determines
// how many elements to read
SymbolNumMap *Matrix::readRaw(ElementType labeled, SymbolNumMap *syms)
{
    int numread;

    // allocate a SymbolNumMap if needed but not supplied
    if (syms==NULL && (labeled==LABELEDROW || labeled==STRINGS)) {
        syms = new SymbolNumMap("", 1000);
    }

    // get matrix
    const int bufferSize=4096;   // buffer
    char buffer[bufferSize];

    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc; c++) {

            // read in a string?
            if ((labeled==LABELEDROW && c==0) || labeled==STRINGS) {
                numread = scanf("%s", buffer);
                if (numread==EOF) {
                    if (name.length()==0) {
                        printf("ERROR(read): Trying to read element [%d, %d] of a matrix but end of file was found\n", r, c);
                    }
                    else {
                        printf("ERROR(read): Trying to read element [%d, %d] of matrix named \"%s\" but end of file was found\n", r, c, name.c_str());
                    }
                    exit(1);
                }
                m[r][c] = syms->add(buffer);   // effectively numeric pointer to the syms
            }

            // read in a number?
            else {
                numread = scanf("%lf", &(m[r][c]));
                if (numread==EOF) {
                    if (name.length()==0) {
                        printf("ERROR(read): Trying to read element [%d, %d] of a matrix but end of file was found\n", r, c);
                    }
                    else {
                        printf("ERROR(read): Trying to read element [%d, %d] of matrix named \"%s\" but end of file was found\n", r, c, name.c_str());
                    }
                    exit(1);
                }
                if (numread!=1) {
                    printf("ERROR(read): invalid number when trying to read row: %d and col: %d.  First character is '%c'\n", r, c, getchar());
                    exit(1);
                }
            }
        }
    }

    defined = true;

    // return the SymbolNumMap
    return syms;
}


// support function for read functions
SymbolNumMap *Matrix::readAux(ElementType labeled, bool transpose, SymbolNumMap *syms)
{
    int r, c;
    int numread;

    // try to read in the number of rows
    numread = scanf("%d", &r);
    if (numread==EOF) {
        if (name.length()==0) {
            printf("ERROR(read): Trying to read a matrix from stdin but end of file was found\n");
        }
        else {
            printf("ERROR(read): Trying to read matrix named \"%s\" from stdin but end of file was found\n", name.c_str());
        }
        exit(1);
    }
    if (numread!=1) {
        printf("ERROR(read): number of rows was not a valid integer.  First character is '%c'\n", getchar());
        exit(1);
    }

    // try to read in the number of columns
    numread = scanf("%d", &c);
    if (numread==EOF) {
        if (name.length()==0) {
            printf("ERROR(read): Trying to read a matrix from stdin but end of file was found\n");
        }
        else {
            printf("ERROR(read): Trying to read matrix named \"%s\" from stdin but end of file was found\n", name.c_str());
        }
        exit(1);
    }
    if (numread!=1) {
        printf("ERROR(read): number of columns was not a valid integer.  First character is '%c'\n", getchar());
        exit(1);
    }


    // transpose?
    if (transpose) {
        int tmp;
        tmp = r; r = c; c = tmp;
    }

    // space allocation
    if (maxr!=r || maxc!=c) {
        reallocate(r, c, name);
    }

    return readRaw(labeled, syms);
}    


void Matrix::readRaw()
{
    readRaw(NUM, NULL);
}


// tri-diagonalize a symmetric matrix.  The matrix will be destroyed and
// the diagonal will be returned in d and off diagonal in e.   It uses
// the Householder transformation
// WARNING: ALLOCATES SPACE for d and e!!!
void Matrix::tridiagonalize(double *&d, double *&e)
{
    d = new double [maxc];  // the diagonal elements
    e = new double [maxc];  // the off-diagonal elements

    householder(m, maxc, d, e);
}


// compute the eigenvalues and eigenvectors of a SYMMETRIC matrix
//
// input:  A symmetric matrix (WARNING: this fact is NOT verified!)
//         Input matrix is destroyed.
// output: Returns a row matrix of eigen values and
//         transforms the matrix into a matrix of eigenvectors (one vector in each row)
//         in the same order as the eigenvalues.
// NOTE: The matrix is NOT the set of eigen vectors in columns!
// NOTE: The vectors are SORTED by decreasing magnitude of eigenvalue
// NOTE: The vectors are normalized (norm of each vector is 1)
// NOTE: The input matrix is destroyed

// WARNING: allocates a row matrix of eigenvalues

// quick insertion sort for the eigen values AND corresponding vectors
// in DECREASING order of MAGNITUDE
void isort(double a[], double *b[], int len)
{
    for (int i=1; i<len; i++) {
        double aa, *bb;
        int j;

        aa = a[i];
        bb = b[i];
        for (j = i-1; (j>=0) && (fabs(a[j])<fabs(aa)); j--) {
            a[j+1] = a[j];
            b[j+1] = b[j];
        }
        a[j+1] = aa;
        b[j+1] = bb;
    }
}

// Destroys self by replacing self with EIGENVECTORS in rows.
// Returns a new matrix (a row vector) with the EIGENVALUES in it.
// Eigenvalues and vectors returned sorted from largest magnitude to smallest
// WARNING: allocates new matrix for answer
Matrix Matrix::eigenSystem()
{
    assertDefined("eigenSystem");
    assertSquare("eigenSystem");

    Matrix values(1, maxc);   // allocates space for eigen values

    {
        double *d, *e;

        tridiagonalize(d, e);         // allocates space for 2 double arrays
        eigen(d, e, maxc, m);         // returns eigen values in d

        delete values.m[0];           // save the eigen values from above routines
        values.m[0] = d;

        delete [] e;
    }

    transposeSelf();              // result was returned in columns so put it in rows
    isort(values.m[0], m, maxc);  // sort both eigenvalues and vectors in rows so vector with max eigenvalue is first
    values.defined = true;        // mark eigen value matrix as defined

    return values;
}


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
//
// The following routines are from "Numerical Recipes in C"
//
// REF: Eigenvalue solvers, tred2 (householder) and tqli (eigen), from
// "Numerical Recipes in C" (Cambridge Univ. Press) by W.H. Press,
// S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery. Translated from
// 1 based array (originally from FORTRAN) by Robert Heckendorn,
// University of Idaho
//
// Householder reduction of a real, symmetric matrix a[0..n-1][0..n-1] to a
// symmetric tridiagonal matrix.
//
// On output, a is replaced by the orthogonal matrix Q effecting the
// transformation. d[0..n-1] returns the diagonal elements of the tridiagonal matrix,
// and e[0..n-1] the off-diagonal elements, with e[0]=0. Several statements, as noted
// in comments, can be omitted if only eigenvalues are to be found, in which case a
// contains no useful information on output. Otherwise they are to be included.
//
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static void householder(double **a, int n, double d[], double e[])
{
    int l, k, j, i;
    double scale, hh, h, g, f;

    for (i=n-1; i>=1; i--) {
        l=i-1;
        h=scale=0.0;
        if (l > 0) {
            for (k=0; k<=l; k++) {
                scale += fabs(a[i][k]);
            }

            if (scale == 0.0) {             // Skip transformation.
                e[i]=a[i][l];
            }
            else {
                for (k=0; k<=l; k++) {
                    a[i][k] /= scale;       // Use scaled a's for transformation.
                    h += a[i][k]*a[i][k];   // Form sigma in h.
                }
                f=a[i][l];
                g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i]=scale*g;
                h -= f*g;                   // Now h is equation (11.2.4).
                a[i][l]=f-g;                // Store u in the ith row of a.
                f=0.0;
                for (j=0; j<=l; j++) {
                    // Next statement can be omitted if eigenvectors not wanted
                    a[j][i]=a[i][j]/h;      // Store u/H in ith column of a.

                    g=0.0;                  // Form an element of A.u in g.
                    for (k=0; k<=j; k++) {
                        g += a[j][k]*a[i][k];
                    }
                    for (k=j+1; k<=l; k++) {
                        g += a[k][j]*a[i][k];
                    }

                    e[j]=g/h;               // Form element of p in temporarily unused element of e.
                    f += e[j]*a[i][j];
                }
                hh=f/(h+h);                 // Form K, equation (11.2.11).

                // Form q and store in e overwriting p.
                for (j=0; j<=l; j++) {
                    f=a[i][j];
                    e[j]=g=e[j]-hh*f;

                    // Reduce a, equation (11.2.13).
                    for (k=0; k<=j; k++) {
                        a[j][k] -= (f*e[k] + g*a[i][k]);
                    }
                }
            }
        }
        else {
            e[i]=a[i][l];
        }
        d[i]=h;
    }

    // Next statement can be omitted if eigenvectors not wanted
    d[0]=0.0;
    e[0]=0.0;

    // Contents of this loop can be omitted if eigenvectors not
    //   wanted except for statement d[i]=a[i][i];

    // Begin accumulation of transformation matrices.
    for (i=0; i<n; i++) {
        l=i-1;
        if (d[i]) {                        // This block skipped when i=0.
            for (j=0; j<=l; j++) {
                // Use u and u/H stored in a to form P.Q.
                g=0.0;
                for (k=0; k<=l; k++) {
                    g += a[i][k]*a[k][j];
                }
                for (k=0; k<=l; k++) {
                    a[k][j] -= g*a[k][i];
                }
            }
        }
        d[i]=a[i][i];                       // This statement remains.
        a[i][i]=1.0;                        // Reset row and column of a to identity matrix for next iteration.
        for (j=0; j<=l; j++) {
            a[j][i]=a[i][j]=0.0;
        }
    }
}


// Compute the eigen values and vectors of a symmetric tridiagonal matrix
//
// QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors
// of a real, symmetric, tridiagonal matrix, or of a real, symmetric matrix
// previously reduced by householder sec. 11.2. On input, d[0..n-1] contains the diagonal
// elements of the tridiagonal matrix. On output, it returns the eigenvalues. The
// vector e[0..n-1] inputs the subdiagonal elements of the tridiagonal matrix, with
// e[0] arbitrary. On output e is destroyed. When finding only the eigenvalues,
// several lines may be omitted, as noted in the comments. If the eigenvectors of
// a tridiagonal matrix are desired, the matrix z[0..n-1][0..n-1] is input as the
// identity matrix. If the eigenvectors of a matrix that has been reduced by householder
// are required, then z is input as the matrix output by householder. In either case,
// the kth column of z returns the normalized eigenvector corresponding to d[k].
//
// input: d - diagonal of symmetric tridiagonal matrix
//        e - offdiagonal of symmetric tridiagonal matrix
//        z - identity if you want eigensystem of symmetric tridiagonal matrix
//          - OR the householder reduction of a symmetric matrix
// output: d - eigenvalues
//         z - the corresponding eigen vectors in the COLUMNS!!!
static void eigen(double *d, double *e, int n, double **z)
{
    double pythag(double a, double b);
    int m, l, iter, i, k;
    double s, r, p, g, f, dd, c, b;

      // Convenient to renumber the elements of e.
    for (i=1; i<n; i++) e[i-1]=e[i];
    e[n-1]=0.0;

    for (l=0; l<n; l++) {
        iter=0;
        do {
            // Look for a single small subdiagonal element to split the matrix.
            for (m=l; m<n-1; m++) {
                dd=fabs(d[m])+fabs(d[m+1]);
                if ((double)(fabs(e[m])+dd) == dd) break;
            }

            if (m != l) {
                if (iter++ == 30) printf("Too many iterations in tqli");
                g=(d[l+1]-d[l])/(2.0*e[l]);       // Form shift.
                r=pythag(g, 1.0);
                g=d[m]-d[l]+e[l]/(g+SIGN(r, g));       // This is dm - ks.
                s=c=1.0;
                p=0.0;
                for (i=m-1; i>=l; i--) {      // A plane rotation as in the original QL, followed by Givens
                    f=s*e[i];                // rotations to restore tridiagonal form.
                    b=c*e[i];
                    e[i+1]=(r=pythag(f, g));
                    if (r == 0.0) {      // Recover from underflow.
                        d[i+1] -= p;
                        e[m]=0.0;
                        break;
                    }
                    s=f/r;
                    c=g/r;
                    g=d[i+1]-p;
                    r=(d[i]-g)*s+2.0*c*b;
                    d[i+1]=g+(p=s*r);
                    g=c*r-b;
                    // Next loop can be omitted if eigenvectors not wanted
                    // Form eigenvectors.
                    for (k=0; k<n; k++) {
                        f=z[k][i+1];
                        z[k][i+1]=s*z[k][i]+c*f;
                        z[k][i]=c*z[k][i]-s*f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l]=g;
                e[m]=0.0;
            }
        } while (m != l);
    }
}


//******************************************************************************
// Computes (a2 + b2)1/2 without destructive underflow or overflow.
//
double pythag(double a, double b)
{
    double absa, absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
//
// translated from gaussj in Numerical Recipes in C, pg 36
//
// Linear equation solution by Gauss-Jordan elemination. a[0..n-1][0..n-1]
// is an input matrix of n by n elements. b[0..n-1][0..m-1] is an input
// matrix of n by m containing the m right-hand side vectors. On output, a is
// replaced by its matrix inverse, and b is replaced by the corresponding set
// of solution vectors.
//
// returns true if successful and returns false if matrix is singular
#define SWAP(a,b) {double temp=(a); (a)=(b); (b)=temp; }
static bool gaussj(double **a, int n, double **b, int m)
{
    int *ipiv;
    int *indxr, *indxc;
    int i, j, k, l, ll;
    int icol=-666, irow=-666;  // init these so uninitialized warning does not occur
    double big, pivinv;

    // allocate some temp space for pivoting (in C++11 we could allocate in declaration)
    // Jack A. suggests using a vector anyway.
    ipiv = new int[n];
    indxr = new int[n];
    indxc = new int[n];

    for (j=0; j<n; j++) ipiv[j] = 0;

    // This is the main loop over the columns to be reduced
    for (i=0; i<n; i++) {
        big = 0.0;

        // select the largest a[irow][icol]
        for (j=0; j<n; j++) {
            if (ipiv[j] != 1)  {
                for (k=0; k<n; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else {
                        if (ipiv[k] > 1) {
                            printf("SYSTEM ERROR(gaussj): Singular Matrix type 1\n");
                            fflush(stdout);
                            return false;
                        }
                    }
                }
            }
        }
        ipiv[icol]++;


        // We now have the pivot element, so we intechange rows, if needed
        // to put the pivot element on the diagonal. The columns are not
        // physically interchanged, only relabeled: indxc[i], the column
        // of the ith pivot element, is the ith column that is reduced, while
        // indxr[i] is the row in which that pivot element was originally located.
        // If indxr[i] != indxc[i] there is an implied column interchange.
        // With this form of bookkeeping, the solution b's will end up in the
        // correct order, and the inverse matrix will be scrambled by columns
        if (irow != icol) {
            for (l=0; l<n; l++) SWAP(a[irow][l],a[icol][l]);
            for (l=0; l<m; l++) SWAP(b[irow][l],b[icol][l]);
        }

        // We are now ready to divide the pivot row by the pivot element
        // located at irow and icol
        indxr[i] = irow;
        indxc[i] = icol;
// zzz    catastrophic subtraction can happen with
//        double data[] = {0, 1, 2, 3,
//                 4, 5, 6, 7,
//                 8, 9, 10, 11};
//        printf("A: %lg\n", a[icol][icol]);
//            fflush(stdout);
        if (a[icol][icol] == 0.0) {
            return false;   // matrix is singular!!
        }

        pivinv = 1.0/a[icol][icol];
        a[icol][icol] = 1.0;
        for (l=0; l<n; l++) a[icol][l] *= pivinv;
        for (l=0; l<m; l++) b[icol][l] *= pivinv;

        // next, we reduce the rows ...
        // .. except for the pivot one, of course.
        // DANGER: beware of catastrophic subtraction!!! -RH
        for (ll=0; ll<n; ll++) {
            if (ll != icol) {
                double dum;

                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l=0; l<n; l++) {
                    a[ll][l] -= a[icol][l]*dum;
                }
                for (l=0; l<m; l++) {
                    b[ll][l] -= b[icol][l]*dum;
                }
            }

            //RH zzz put *real* test for catastrophic subtraction here.  This is bogus.
            if (fabs(a[ll][ll]) < 1E-14) a[ll][ll]=0;
            //RH printf("LLA: %lg\n", a[ll][ll]);
        }
    }

    // This is the end of the main loop over columns of the reduction.
    // it only remains to unscramble the solution in view of the column
    // interchanges. We do this by interchanging pairs of columns in the reverse
    // order that the permutation was built.
    for (l=n-1; l>=0; l--) {
        if (indxr[l] != indxc[l]) {
            for (k=0; k<n; k++) SWAP(a[k][indxr[l]], a[k][indxc[l]]);
        }
    }

    return true;
}




// use this n^2 sort for small numbers of elements
// to sort the rows of a matrix from rows numbered:
// lower to upper inclusive
void Matrix::selectSort(int lower, int upper)
{
    int bestLoc;

    for (int l=lower; l<upper; l++) {
        bestLoc = l;
        for (int u=l+1; u<=upper; u++) {
            if (isLessRows(u, bestLoc)) bestLoc=u;
        }
        if (bestLoc!=l) swapRows(l, bestLoc);
    }
}


// use this n^2 sort for small numbers of elements
// to sort the rows of a matrix from rows numbered:
// lower to upper inclusive
void Matrix::selectSortCol(int c, int lower, int upper)
{
    int bestLoc;

    for (int l=lower; l<upper; l++) {
        bestLoc = l;
        for (int u=l+1; u<=upper; u++) {
            if (m[u][c] < m[bestLoc][c]) bestLoc=u;
        }
        if (bestLoc!=l) swapRows(l, bestLoc);
    }
}


// do a quick sort of elements a[lower]...a[upper]
void Matrix::qs(int lower, int upper)
{
    int save;

    if (upper-lower<32) {
        selectSort(lower, upper);
        return;
    }

    save = lower;      // start scanning up with ptr from lower to upper
    for (int ptr=lower; ptr<upper; ptr++) {
        // keep shuffling the values less that split to be in the
        // array indexed below save.
        if (isLessRows(ptr, upper)) {
            swapRows(save, ptr);
            save++;
        }
    }

    swapRows(upper, save);

    if (save-lower>1) qs(lower, save-1);
    if (upper-save>1) qs(save+1, upper);
}


// do a quick sort of elements a[lower]...a[upper]
void Matrix::qsCol(int c, int lower, int upper)
{
    int save;

    if (upper-lower<32) {
        selectSortCol(c, lower, upper);
        return;
    }

    save = lower;      // start scanning up with ptr from lower to upper
    for (int ptr=lower; ptr<upper; ptr++) {
        // keep shuffling the values less that split to be in the
        // array indexed below save.
        if (m[ptr][c] < m[upper][c]) {
            swapRows(save, ptr);
            save++;
        }
    }

    swapRows(upper, save);

    if (save-lower>1) qsCol(c, lower, save-1);
    if (upper-save>1) qsCol(c, save+1, upper);
}


// sort the rows of a matrix.   WARNING: SORTS IN PLACE
Matrix &Matrix::sortRows() {
    assertDefined("sortRows");
    if (maxr>1) qs(0, maxr-1);

    return *this;
}


// do a quick sort of elements a[lower]...a[upper]
// k is the smallest index of the row that is to be sorted into.
// That is, the maximum numRows-k values will be put in the upper
// elements of the matrix
void Matrix::qsMaxK(int k, int lower, int upper)
{
    int save;

    if (upper<k) return;   // don't sort what you don't need

    if (upper-lower<4) {
        selectSort(lower, upper);
        return;
    }

    save = lower;      // start scanning up with ptr from lower to upper
    for (int ptr=lower; ptr<upper; ptr++) {
        // keep shuffling the values less that split to be in the
        // array indexed below save.
        if (isLessRows(ptr, upper)) {
            swapRows(save, ptr);
            save++;
        }
    }

    swapRows(upper, save);

    if (save-lower>1) qsMaxK(k, lower, save-1);
    if (upper-save>1) qsMaxK(k, save+1, upper);
}


void Matrix::qsMaxKCol(int k, int c, int lower, int upper)
{
    int save;

    if (upper<k) return;   // don't sort what you don't need

    if (upper-lower<4) {
        selectSortCol(c, lower, upper);
        return;
    }

    save = lower;      // start scanning up with ptr from lower to upper
    for (int ptr=lower; ptr<upper; ptr++) {
        // keep shuffling the values less that split to be in the
        // array indexed below save.
        if (m[ptr][c] < m[upper][c]) {
            swapRows(save, ptr);
            save++;
        }
    }

    swapRows(upper, save);

    if (save-lower>1) qsMaxKCol(k, c, lower, save-1);
    if (upper-save>1) qsMaxKCol(k, c,  save+1, upper);
}


// sort the rows of a matrix so the max k values are in the k
// highest indexes.   That is, it does a partial sort to find
// the largest k values.   WARNING: SORTS IN PLACE
Matrix &Matrix::maxKRows(int k) {
    assertDefined("maxKRows");
// zzz test that k is reasonable
    if (maxr>1) qsMaxK(maxr-k, 0, maxr-1);

    return *this;
}


Matrix &Matrix::maxKRowsByCol(int k, int c) {
    assertDefined("maxKRowsByCol");
// zzz test that k is reasonable
    if (maxr>1) qsMaxKCol(maxr-k, c, 0, maxr-1);

    return *this;
}



Matrix &Matrix::sortRows(int startRow, int endRow) {
    assertDefined("sortRows");
    assertRowIndexOK(startRow, "sortRows");
    assertRowIndexOK(endRow, "sortRows");
    if (maxr>1) qs(startRow, endRow);

    return *this;
}


// sort the rows of a matrix using column c as the key.
// Column numbering starts at 0.
// WARNING: sorts in place
Matrix &Matrix::sortRowsByCol(int c) {
    assertDefined("sortRowsCol");
    assertColIndexOK(c, "sortRowsByCol");
    if (maxr>1) qsCol(c, 0, maxr-1);

    return *this;
}


// sort rows in place in a range of rows
// Column numbering starts at 0.
// WARNING: sorts in place
Matrix &Matrix::sortRowsByCol(int c, int startRow, int endRow)
{
    assertDefined("sortRowsByCol");
    assertColIndexOK(c, "sortRowsByCol");
    assertRowIndexOK(startRow, "sortRowsByCol");
    assertRowIndexOK(endRow, "sortRowsByCol");
    if (maxr>1) qsCol(c, startRow, endRow);

    return *this;
}



// Create a subMatrix.   DANGER: This bit of evil is a matrix that POINTS
// INTO ANOTHER MATRIX!   DANGER: Do not use the subMatrix after you
// deallocate the other matrix!!   In a sense this is not a real matrix.
// If you want this matrix to persist then you have to make a copy of it.
// It is great for efficiency.  BE SURE to deallocate this matrix before you
// deallocate the "mother matrix"
Matrix Matrix::subMatrix(int minr, int minc, int sizer, int sizec) const
{
    if (sizer==0) sizer = maxr - minr;
    if (sizec==0) sizec = maxc - minc;

    assertIndexOK(minr, minc, "lower bounds submatrix extract");
    assertIndexOK(minr+sizer-1, minc+sizec-1, "upper bounds submatrix extract");

    Matrix out(sizer);                         // allocate a subMatrix!
    out.maxc = sizec;                          // fix internal column width

    for (int r=0; r<sizer; r++) {
        out.m[r] = &(m[minr][minc]);               // DANGER: we are copying pointers into other Matrix!!!
        minr++;
    }

    out.defined = true;

    return out;
}



// Create a subMatrix.   DANGER: This bit of evil is a matrix that POINTS
// INTO ANOTHER MATRIX!   DANGER: Do not use the subMatrix after you
// deallocate the other matrix!!   In a sense this is not a real matrix.
// If you want this matrix to persist then you have to make a full copy of it.
Matrix Matrix::subMatrixEq(int c, double value) const
{
    assertColIndexOK(c, "subMatrixEq");

    std::vector<double *> rowList;        // this is a retrofit of using an array originally when vector better

    for (int r=0; r<maxr; r++) {
        if (m[r][c]==value) rowList.push_back(m[r]);
    }

    Matrix out(rowList.size(), "sub" + name);                         // allocate a subMatrix!
    out.maxc = maxc;

    for (unsigned int r=0; r<rowList.size(); r++) {
        out.m[r] = rowList[r];                          // DANGER: we are copying pointers into other Matrix!!!
    }

    out.defined = true;

    return out;
}


// Create a subMatrix.   DANGER: This bit of evil is a matrix that POINTS
// INTO ANOTHER MATRIX!   DANGER: Do not use the subMatrix after you
// deallocate the other matrix!!   In a sense this is not a real matrix.
// If you want this matrix to persist then you have to make a full copy of it.
Matrix Matrix::subMatrixNeq(int c, double value) const
{
    assertColIndexOK(c, "subMatrixNeq");

    std::vector<double *> rowList;        // this is a retrofit of using an array originally when vector better

    for (int r=0; r<maxr; r++) {
        if (m[r][c]!=value) rowList.push_back(m[r]);
    }

    Matrix out(rowList.size(), "sub" + name);                         // allocate a subMatrix!
    out.maxc = maxc;

    for (unsigned int r=0; r<rowList.size(); r++) {
        out.m[r] = rowList[r];                          // DANGER: we are copying pointers into other Matrix!!!
    }

    out.defined = true;

    return out;
}


// Create a subMatrix.   DANGER: This bit of evil is a matrix that POINTS
// INTO ANOTHER MATRIX!   DANGER: Do not use the subMatrix after you
// deallocate the other matrix!!   In a sense this is not a real matrix.
// If you want this matrix to persist then you have to make a full copy of it.
Matrix Matrix::subMatrixPickRows(int match, const Matrix &list)
{
    assertDefined("lhs of subMatrixPickRows");
    list.assertDefined("rhs of subMatrixPickRows");
    list.assertColVector("rhs of subMatrixPickRows");
    assertRowsEqual(list, "subMatrixPickRows");

    std::vector<double *> rowList;        // this is a retrofit of using an array originally when vector better
    for (int r=0; r<maxr; r++) {
        if (list.m[r][0]==match) rowList.push_back(m[r]);
    }

    if (rowList.size()==0) {
        Matrix out("subMatrixPickRows Found No Rows In "+name);

        return out;
    }

    Matrix out(rowList.size());                         // allocate a subMatrix!
    out.maxc = maxc;

    for (unsigned int r=0; r<rowList.size(); r++) {
        out.m[r] = rowList[r];                          // DANGER: we are copying pointers into other Matrix!!!
    }

    out.defined = true;

    return out;
}




// // // // // // // // // // // // // // // // // // // //
//
// image (picture) support (currently just pgm files)
//
// image (picture) support (currently only supports 8 bit pgm and ppm formats)
// output is in ascii formats (zzz: fix someday to use more compressed output)
// 8 bit gray is one integer in the range 0-255 for each pixel
// 8 bit color is three integers in a row in the range 0-255 for RGB in each pixel.
// That is an 8 bit color square 100x100 pixels gens a 100x300 dimensional array

// helper routine for writing images
int Matrix::byteValue(double x)
{
    int z;

    z = int(x);
    if (z<0) z = 0;
    if (z>255) z = 255;

    return z;
}

// helper routine for reading images
Matrix &Matrix::readImage(std::string expectedType,
                         std::string caller,
                         std::string filename,
                         std::string namex,
                         bool &isColor)
{
    char magic[3];               // magic number
    const int bufferSize=4096;   // buffer
    char buffer[bufferSize];
    FILE *IN;                    // input file
    int newr, newc, max;               // picture parms

    if (filename.length()>0) {
        IN = fopen(filename.c_str(), "r");
        if (IN==NULL) {
            printf("ERROR(%s): Trying to open file \"%s\" but failed.\n", caller.c_str(), filename.c_str());
            exit(1);
        }
    }
    else {
        IN = stdin;
    }

    // get magic number of file
    if (fscanf(IN, "%2s", magic)!=1) {
        if (filename.length()>0) {
            printf("ERROR(%s): Unable to read file magic number for file named \"%s\".\n", caller.c_str(), filename.c_str());
        }
        else {
            printf("ERROR(%s): Unable to read file magic number for file from stdin.\n", caller.c_str());
        }

        exit(1);
    }

    if (! (magic[0]=='P' && strisin(magic[1], expectedType.c_str()))) {
        if (filename.length()>0) {
            printf("ERROR(%s): Wrong magic number for file named \"%s\".  Expecting P%c or P%c format but had magic number: \"%c%c\".\n",
                   caller.c_str(),
                   filename.c_str(),
                   expectedType[0],
                   expectedType[1],
                   magic[0],
                   magic[1]);
        }
        else {
            printf("ERROR(%s): Wrong magic number for file on stdin.  Expecting P%c or P%c format but had magic number: \"%c%c\".\n",
                   caller.c_str(),
                   expectedType[0],
                   expectedType[1],
                   magic[0],
                   magic[1]);
        }

        exit(1);
    }

    // read comment lines (Warning: assumes comments come right after magic number)
    fscanf(IN, "%s", buffer);
    while (*buffer=='#') {
            fgets(buffer, bufferSize, IN);
//            printf("# %s", buffer);
            fscanf(IN, "%s", buffer);
    }

    // read picture parameters
    newc = atoi(buffer);      // number of cols
    fscanf(IN, "%d", &newr);  // number of rows
    fscanf(IN, "%d", &max);   // maximum value for pixel

    // is color?
    isColor = (magic[1]=='3') || (magic[1]=='6');
    if (isColor) newc *= 3;

    // reallocate myself
    reallocate(newr, newc, namex);

    // is ascii numbers?
    if (magic[1]=='2' || magic[1]=='3') {
        for (int r=0; r<maxr; r++) {
            for (int c=0; c<maxc; c++) {
                int tmp;
                if (fscanf(IN, "%d", &tmp)!=1) {
                    printf("ERROR(%s): Trying to read ascii pixel value at position (%d, %d) from file \"%s\" but failed.\n",
                           caller.c_str(),
                           r, c,
                           filename.c_str());
                    exit(1);
                }
                m[r][c] = tmp;
            }
        }
    }

    // is binary numbers?
    if (magic[1]=='5' || magic[1]=='6') {
        getc(IN);
        for (int r=0; r<maxr; r++) {
            for (int c=0; c<maxc; c++) {
                int byte;
                byte = getc(IN);
                if (byte==EOF) {
                    printf("ERROR(%s): Trying to read a byte of pixel value at position (%d, %d) from file \"%s\" but got EOF.\n",
                           caller.c_str(),
                           r, c,
                           filename.c_str());
                    exit(1);
                }

//                    if (feof(IN)) printf("ERROR() EOF\n");
//                    if (ferror(IN)) {
//                        printf("ERROR() ferror\n");
//                        perror("What's up:");
//                    }
                m[r][c] = byte;
            }
        }
    }

    defined = true;

    return *this;
}



// read in a netpbm format file (Only P2, P3, P5, or P6 magic numbers
// will be accepted) from file named filename. If filename is the
// empty string, "" then the data will be read from stdin. namex is
// the name of the matrix created. isColor is a bool variable passed
// in by reference and returns the whether the file was in ppm (color
// which is P3 or P6) or pgm (gray-scale which is P2 or P5). If a
// grayscale file is read, a matrix will be created with the numeric
// value of each pixel and isColor will be set to false. If color, RGB
// will be read in as 3 consecutive numbers in the matrix and isColor
// will be set to true.
Matrix &Matrix::readImagePixmap(std::string filename, std::string namex, bool& isColor)
{
    return readImage("2356", "readImagePixmap", filename, namex, isColor);  // accept types P2, P3, P5, or P6
}


// Read a pgm file  (8 bit gray scale) in P2 or P5 format.
// WARNING: crudely assumes comments are less than 4K bytes
Matrix &Matrix::readImagePgm(std::string filename, std::string namex)
{
    bool isColor;  // space holder
    return readImage("25", "readImagePgm", filename, namex, isColor);  // accept types P2 or P5
}


// Read specifically a ppm (8 bit color) format file
// WARNING: crudely assumes comments are less than 4K bytes
Matrix &Matrix::readImagePpm(std::string filename, std::string namex)
{
    bool isColor;  // space holder
    return readImage("36", "readImagePpm", filename, namex, isColor);  // accept types P3 or P6
}



// Write a pgm file  (8 bit gray scale)
// It uses the readable character P2 representation of a picture, rather than
// the binary P5 representation.  Line length is unrestricted.
// WARNING: the user is entrusted with the task of using the pgm file extension
// in the filename
void Matrix::writeImagePgm(std::string filename, std::string comment)
{
    FILE *OUT;

    assertDefined("writeImagePgm");
    if (filename.length()>0) {
        OUT = fopen(filename.c_str(), "w");
        if (OUT==NULL) {
            printf("ERROR(writeImagePgm): Trying to open file \"%s\" but failed.\n", filename.c_str());
            exit(1);
        }
    }
    else {
        OUT = stdout;
    }

    fprintf(OUT, "P2\n");
    if (name.length()>0) fprintf(OUT, "# Name: %s\n", name.c_str());
    if (comment.length()>0) fprintf(OUT, "# %s\n", comment.c_str());
    fprintf(OUT, "# 8 bit gray scale\n");
    fprintf(OUT, "%d %d\n", maxc, maxr);        // NOTE: columns then rows!
    fprintf(OUT, "255\n");                      // maximum level of gray
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc-1; c++) {
            fprintf(OUT, "%d ", byteValue(m[r][c]));
        }
        fprintf(OUT, "%d\n", byteValue(m[r][maxc-1]));
    }

    fclose(OUT);
}



// Write a ppm file  (8 bit gray scale)
// It uses the readable character P3 representation of a picture, rather than
// the binary P6 representation.  Line length is unrestricted.
// WARNING: the user is entrusted with the task of using the ppm file extension
// in the filename
void Matrix::writeImagePpm(std::string filename, std::string comment)
{
    FILE *OUT;

    assertDefined("writeImagePpm");
    if (maxc%3 != 0) {
        if (name.length()==0) {
            printf("ERROR(writeImagePpm): Number of columns %d not divisible by three but supposed to be a matrix of RGB values.\n", maxc);
        }
        else {
            printf("ERROR(writeImagePpm): Number of columns %d in matrix named \"%s\" is not divisible by three but supposed to be a matrix of RGB values .\n", maxc, name.c_str());
        }
        exit(1);
    }

    if (filename.length()>0) {
        OUT = fopen(filename.c_str(), "w");
        if (OUT==NULL) {
            printf("ERROR(writeImagePpm): Trying to open file \"%s\" but failed.\n", filename.c_str());
            exit(1);
        }
    }
    else {
        OUT = stdout;
    }

    fprintf(OUT, "P3\n");
    if (name.length()>0) fprintf(OUT, "# Name: %s\n", name.c_str());
    if (comment.length()>0) fprintf(OUT, "# %s\n", comment.c_str());
    fprintf(OUT, "# 8 bit color\n");
    fprintf(OUT, "%d %d\n", maxc/3, maxr);      // NOTE: columns then rows!
    fprintf(OUT, "255\n");                      // maximum level of color channels
    for (int r=0; r<maxr; r++) {
        for (int c=0; c<maxc-1; c++) {
            fprintf(OUT, "%d ", byteValue(m[r][c]));
        }
        fprintf(OUT, "%d\n", byteValue(m[r][maxc-1]));
    }

    fclose(OUT);
}


// // // // // // // // // // // // // // // // // // // // // // // //
//
// Some random tests for the matrix code
//


/* UNCOMMENT TO HAVE A MAIN FOR TESTING */
/*

int main()
{
    Matrix x(10, 20);
    Matrix z("dogs");
    Matrix y(2, 3, yvalues, "cats");
    Matrix a, b;

    Matrix xx = new Matrix(3, 4);

    initRand();

    y.print("matrix y");
//    z.print(); // undefined
//    x.print("matrix x");  // undefined
//    y.dot(y);   // wrong sizes

    printf("Supply a 3x3 or larger matrix to read:\n");
    x.read();

    printf("Read Matrix\n");
    x.write();

    a = x.transpose();
    printf("Transpose\n");
    a.write();

    a.mul(a);
    printf("Squared\n");
    a.write();

    b = a.extract(1, 1, 2, 1);
    printf("Extracted\n");
    b.write();

    a.insert(b, 0, 0);
    printf("Inserted\n");
    a.write();

    printf("\n");
    x.print();
    printf("\n");
    a.print();
    printf("\n");
    b.print();
    printf("\n");

    (x.extract(0, 2, 0, 0)).print();
    printf("\n");

    printf("\n");
    (x.extract(0, 0, 0, 2)).print();
    printf("\n");

    printf("Map\n");
    x.map(f);
    x.write();
    printf("\n");

    printf("Random -1.0 to 1.0\n");
    x.rand(-1.0, 1.0);
    x.write();
    printf("\n");

    printf("Normalize\n");
    x.write();
    printf("\n");
    a = x.normalizeCols();
    x.write();
    printf("\n");
    a.write();
    printf("\n");
    b.normalizeCols(a);
    b.write();
    printf("\n");

    printf("Random -5.0 to 10.0\n");
    x.rand(-5.0, 10.0);
    x.write();
    printf("\n");

    printf("Random -5 to 10\n");
    x.rand(-5, 10);
    x.write();
    printf("\n");

    printf("Random 0 to 6\n");
    x.rand(0, 6);
    printf("X\n");
    x.write();
    printf("\n");

    // print b
    b = x.transpose();
    printf("B\n");
    b.write();
    printf("\n");

    // print x . b
    printf("X.B\n");
    a = x.dot(b);
    a.write();
    printf("\n");

    // print b . x
    printf("B.X\n");
    a = b.dot(x);
    a.write();

    // constant
    a.constant(3.14159265);
    a.write();

    {
        Matrix m(5, 3);

        initRand();

        m.rand(0, 10);
        m.print("");

        MatrixRowIter a(&m);
        for (Matrix *i = a.rowBegin(); a.rowNotEnd(); a.rowNext()) {
            i->print("");
        }

        return 0;
    }
}
*/
/*
double yvalues[] = {2, 3, 5, 7, 11, 13};

int main()
{
    Matrix y(2, 3, yvalues, "cats");
    Matrix x(2, 3, yvalues, "newCats");
    Matrix z(10, 3, "sample");

    initRand();

    y.print();
    x.print();

    (x.Tdot(y)).print();
    x.print();
    x.write();

    x.sample(z);
    z.print();

    return 0;
}

*/

/*
double yvalues[] =  {5, 0, 3, 7, 1, -5, 7, 3, 4, 9, 8, 10};
double avalues[] = {2, 0, -9, 3, 4, 1};
double bvalues[] = {5, 2, 6, -4, 4, 9};
double cvalues[] = {2, 0, -9, 3, 4, 1};
double dvalues[] = {5, 2, 6, 8, -4, 4, 9, 7};
double evalues[] = {1, 1, 2, 2, 3, 1, 4, 2, 5, 1};

int main()
{
//    Matrix::debug = true;
    Matrix x(3, 4,  yvalues, "x");
    Matrix a(2, 3, avalues, "a");
    Matrix b(2, 3, bvalues, "b");
    Matrix c(2, 3, cvalues, "c");
    Matrix d(2, 4, dvalues, "d");
    Matrix e(5, 2, evalues, "e");
    Matrix y("y");

    x.print();
    x.subMatrix(1, 1, 2, 2).transposeSelf();
    x.print();

    e.print();
    e.subMatrixEq(1, 1).print();
    e.subMatrixEq(1, 3).print();
    e.subMatrixNeq(0, 3).print();
    printf("%d\n", e.countNeqCol(0, 1));
    printf("%d\n", e.countNeqCol(0, 2));
    printf("%d\n", e.countNeqCol(0, 3));
    printf("%d\n", e.countNeqCol(0, 8));

    x.covMatrix().print("covMatrix(x)");
    {
        y = e.covMatrix();
    }
    y.print();

    a.print();
    b.print();
    a.covMatrix(b).print("covMatrix(a, b)");
    b.covMatrix(a).print("covMatrix(b, a)");

    c.print();
    d.print();
    c.covMatrix(d).print("covMatrix(c, d)");
    d.covMatrix(c).print("covMatrix(d, c)");

    return 0;
}
*/
/*
int main()
{
    Matrix pic(300, 300, "picture");

    for (int c=0; c<pic.numCols(); c++) {
        pic.constantColRange(c, 0, 1);
    }
    pic.writeImagePgm("zfade.pgm", "fade.pgm");

    pic.readImagePgm("znano.pgm", "nano");
    pic.print();

    pic.readImagePgm("mondrianRedBlueAndYellow.pgm", "mondrian");

    pic.readImagePpm("mondrianRedBlueAndYellow.ppm", "mondrian");
    pic.writeImagePgm("zm.pgm", "mondrian.pgm");

    pic.readImagePpm("girlWithPearlEarringSm.ppm", "mondrian");
    pic.writeImagePgm("zg.pgm", "girl with pearl earring color");

//    pic.readImagePgm("z.pgm", "mondrian").printInt();
}
*/
