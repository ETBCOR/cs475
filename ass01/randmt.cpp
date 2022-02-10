#include "rand.h"

//   A version of the 64 bit Mersenne Twister random number generator
//   Robert Heckendorn      2015/7/10 version
//
//   This is a version modified to work with the random number package by
//   Dr. Robert Heckendorn of the University of Idaho.  Below is the
//   original comment block.   This file is mostly an intact version of mt19937-64.c
// 
//   --------------------------------------------------------------------------------
// 
//   A C-program for MT19937-64 (2004/9/29 version).
//   Coded by Takuji Nishimura and Makoto Matsumoto.
//
//   This is a 64-bit version of Mersenne Twister pseudorandom number
//   generator.
//
//   Before using, initialize the state by using init_genrand64(seed)  
//   or init_by_array64(init_key, key_length).
//
//   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
//   All rights reserved.                          
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//
//     3. The names of its contributors may not be used to endorse or promote 
//        products derived from this software without specific prior written 
//        permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//   References:
//   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
//     ACM Transactions on Modeling and 
//     Computer Simulation 10. (2000) 348--357.
//   M. Matsumoto and T. Nishimura,
//     ``Mersenne Twister: a 623-dimensionally equidistributed
//       uniform pseudorandom number generator''
//     ACM Transactions on Modeling and 
//     Computer Simulation 8. (Jan. 1998) 3--30.
//
//   Any feedback is very welcome.
//   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
//   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)



///////////////////////////////////////////////////
//
// constants for the heckendorn rand library
//

// these constants are supplied to make the RNG return numbers
// on the open interval [0, 1)
#define PIDIV2	    1.570796326794896619231321691639751442099L
#define PI  	    3.141592653589793238462643383279502884197L

static double unitizer64 = 1.0/18446744073709551616.0;
static double unitizer64_2 = 2.0/18446744073709551616.0;
static double unitizer64_pi = PI/18446744073709551616.0;


///////////////////////////////////////////////////
//
// constants for the Mersenne Twister RNG
//

#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */


// The RNG state:
static unsigned long long mt[NN];   // The array for the state vector
static int mti=NN+1;   // mti is a state variable saying where the next random number is
                       // mti==NN+1 means mt[NN] is not initialized */



// initializes mt[NN] with a seed and sets mti using a single key
void init_genrand64(unsigned long long seed)
{
    mt[0] = seed;
    for (mti=1; mti<NN; mti++) 
        mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

// initializes mt[NN] with a seed and sets mti using an array with key_length as its length
void init_by_array64(unsigned long long init_key[], unsigned long long key_length)
{
    unsigned long long i, j, k;

    init_genrand64(19650218ULL);
    i=1; j=0;
    k = (NN>key_length ? NN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 3935559000370003845ULL))
          + init_key[j] + j; /* non linear */
        i++; j++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 2862933555777941757ULL))
          - i; /* non linear */
        i++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
    }

    mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */ 
}


void initRand()
{
    unsigned long long key[2];

#ifdef WINDOWS
    key[0] = 1405321245300013ULL;
#else
    key[0] = (unsigned long long int)getpid()*1405321245300013ULL;
#endif
    key[1] = (unsigned long long int)(time(NULL)<<1) | 0x1ULL;
    init_by_array64(key, 2ULL);
}


void initRand(unsigned long long int a, unsigned long long int b)
{
    unsigned long long key[2];

    key[0] = a;
    key[1] = b;
    init_by_array64(key, 2ULL);
}


// Generate a random number [0, 2^64-1]
unsigned long long randULL()
{
    int i;
    unsigned long long x;
    static unsigned long long mag01[2]={0ULL, MATRIX_A};

    if (mti >= NN) { /* generate NN words at one time */

        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == NN+1) 
            init_genrand64(5489ULL); 

        for (i=0;i<NN-MM;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        for (;i<NN-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        x = (mt[NN-1]&UM)|(mt[0]&LM);
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

        mti = 0;
    }
  
    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return x;
}

// return a uniformly distributed random number between 0 and 1
// NOTE: this can be relatively slow since the multiply is a 64 bit
// multiply.
double randUnit()
{
    return randULL()*unitizer64;
}



// return a uniformly distributed random number between -1 and 1
// NOTE: this can be relatively slow since the multiply is a 64 bit
// multiply.
double randPMUnit()
{
    return randULL()*unitizer64_2 - 1.0;
}


// return a uniformly distributed random number between 0 and m-1
int randMod(int m) {
    return randULL()%m;
}


// return two uniformly distributed random numbers in the rand 0 to m-1
// where they are not equal.
void randMod2(int m, int &a, int &b) 
{
    a = randMod(m);
    b = a + randMod(m-1) + 1;
    if (b>=m) b-=m;
}


// return a uniformly distributed random number whose bits lie in the 
//   masked bits.  If the mask is a right justified set of ones this
//   is far faster to generate mod of 2^k than using mod function.
int randMask(unsigned long long int mask) {
    return randULL()&mask;
}


// return true with a probability of prob
bool choose(double prob)
{
    return randUnit()<prob;
}


// returns zero or nonzero 50% of the time
unsigned long long randCoinToss()
{
    return (randULL()&0x40ULL);
}


// return true with a probability of eigth/8
// can be faster than using choose with a double probability
bool choose8(int eigth)
{
    return (randULL()&0x7ULL)<(unsigned long long int)eigth;
}

// for use in returning true with a probability of prob/(2^k)
// Mask must select the lower k bits
bool chooseMask(unsigned long long int mask, int prob)
{
    return (randULL()&mask)<(unsigned long long int)prob;
}


// Random number generator with normal (Gaussian) distribution
// from p 117 of Knuth vol 2 2nd ed.
// Note: this is slow.   Should use the Ziggurat method some day.
static bool gotSpare64=false;
static double spare64;
double randNorm(double stddev)
{
    double u, v, s;

    if (gotSpare64) {
        gotSpare64=false;
        return spare64;
    }
    else {
        do {
            u = 2*randUnit() - 1;
            v = 2*randUnit() - 1;
            s = u*u + v*v;
        } while (s>=1.0 || s==0.0);
    }

    s = sqrt(-2*log(s)/s)*stddev;
    spare64 = v*s;
    gotSpare64 = true;

    return u*s;
}



// Random number generators with a Cauchy distribution
// based on the inversion method using CDF F(x) = .5 + atan(x)/pi
// which yeilds  tan(pi F(x) - .5) as a the tranformation.
// To incorporate a mean and scale: scale*randCauchy()+mean
// 


double randCauchy()
{
    unsigned long long int r;

    do r=randULL(); while (r==0);

    return tan(unitizer64_pi*r - PIDIV2);
}


double randCauchy(double mean, double scale)
{
    unsigned long long int r;

    do r=randULL(); while (r==0);

    return scale*tan(unitizer64_pi*r - PIDIV2) + mean;
}





/*

#include <stdio.h>
#include <stdlib.h>

// Return a uniformly distributed random number between 0 and m-1.
// Suitable for 0 < m < 2^16.
// Uses rejection sampling to avoid the computation of a mod
// Requires a mask be set that of the form (1<<k)-1 that most
// closely contains the value m.   e.g.  randMod(63ULL, m);
// for 32<m<=64.
//
// IMPORTANT: if you are taking something mod m then the mask must be
// m-1 or contain m-1!  For example: if m=4 then mask=3 which is the
// smallest right justified all ones mask containing m-1.  If m=6 then
// mask=7.
int randMod(unsigned long long int mask, int m) {
    int ans;
    do {
        V = W;
        W = X;
        X = ((Y<<41) + (Y>>23)) + Z;
        Y = Z ^ W;
        Z = V + X;

        ans = Y&mask;
        if (ans<m) return ans;

        ans = (Y>>16)&mask;
        if (ans<m) return ans;

        ans = (Y>>32)&mask;
        if (ans<m) return ans;

        ans = (Y>>48)&mask;
    } while (ans>=m);

    return ans;
}


static unsigned long long int smallestMask(unsigned long long int x) {
    x |= x>>1;
    x |= x>>2;
    x |= x>>4;
    x |= x>>8;
    x |= x>>16;
    x |= x>>32;

    return x;
}


int main(int argc, char *argv[]) {
    initRand();
    int m;
    int z;
    unsigned long long int mask;
    z = 0;
    m = atoi(argv[1]);
    mask = smallestMask(m-1);
    for (int i=0; i<100000000; i++) {
//        printf("%d\n", randMod(63ULL, 17));
        z+=randMod(mask, m);
//       z+=randMod(m);
    }
    return z;
}
*/
