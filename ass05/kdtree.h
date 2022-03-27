/*****************************
 * By:			Ethan Corgatelli
 * File:        kdtree.h
 * Project:     Assignment 5
 * Class:       cs475
 * Asn. Pg:     http://marvin.cs.uidaho.edu/Teaching/CS475/pas05.pdf
 *
 *****************************/
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "mat.h"
#include "rand.h"

void build (Matrix *t, int c, int lower, int upper, int i);
