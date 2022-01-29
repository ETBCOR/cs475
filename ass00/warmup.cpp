/*****************************
 * File:	warmup.cpp
 * Project:	assignment 0
 * Class:	cs475
 * Asn. Pg:	http://marvin.cs.uidaho.edu/Teaching/CS475/pas00.pdf
 *
*****************************/
#include "mat.h"
#include "rand.h"

const bool debug = false;

int main (int argc, char *argv[]) {

	if(debug) printf("1 --->\n"); // 1
	Matrix a = Matrix("A");
	Matrix b = Matrix("B");

	if(debug) printf("2 --->\n"); // 2
	int mArr[12];// = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
	for (int i = 0; i < 12; i++) mArr[i] = i;

	Matrix m = Matrix(3, 4, mArr, "M");

	if(debug) printf("3 --->\n"); // 3
	int sArr[9] = {1, 2, 4, 2, 2, 3, 4, 3, 3};
	Matrix s = Matrix(3, 3, sArr, "S");

	if(debug) printf("4 --->\n"); // 4
	s.print();

	if(debug) printf("5 --->\n"); // 5
	a = s;

	if(debug) printf("6 --->\n"); // 6
	s.inverse().print("inverse:");
	
	if(debug) printf("7 --->\n"); // 7
	s.print(); // It has changed. inverse() replaces self
	
	if(debug) printf("8 --->\n"); // 8
	a.dot(s).printFmt("", "%12.4lg", true);

	if(debug) printf("9 --->\n"); // 9
	a.print();
	
	if(debug) printf("10 --->\n"); // 10
	printf("---\n");
	
	if(debug) printf("11 --->\n"); // 11
	a.read();

	if(debug) printf("12 --->\n"); // 12
	a.print();

	if(debug) printf("13 --->\n"); // 13
	m.print();

	if(debug) printf("14 --->\n"); // 14
	printf("---\n");
	
	if(debug) printf("15 --->\n"); // 15
	b = m.dot(a).print("result:");	

	if(debug) printf("16 --->\n"); // 16
	b.print();
	
	if(debug) printf("17 --->\n"); // 17
	printf("---\n");

	if(debug) printf("18 --->\n"); // 18
	b = a.transpose();

	if(debug) printf("19 --->\n"); // 19
	b.print();

	if(debug) printf("20 --->\n"); // 20
	b.dot(b.transpose()).print();

}
