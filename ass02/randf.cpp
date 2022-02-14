#include "rand.h"
#include <stdio.h>
#include <stdlib.h>


// // // // // // // // // // // // // // // // // // // // // // // // //
//                                                                      //
//                   Portable Random Number Kernel
//

//  This random number kernel is designed to provide a very fast
//  random number generator fastRand and a slightly slower but
//  better generator called slowRand.  This is not to say that
//  the fast generator is not a good one (see Meysenburg and
//  Foster, "The Quality of Pseudo-Random Number Generators and
//  SImple Genetic Algorithm Performance", in The Seventh
//  International Conference on Genetic Algorithms, p276-281) We
//  want to be sure that there is no correlation between elements
//  in a population.  Both random number generators I have
//  modified come from references in the above article and were
//  chosen both for their speed and randomness.
//
//  Anecdote: I was running my larger general purpose genetic
//  algorithm program on day long tasks.  I decided to profile
//  it.  Interestingly I discovered that 40% of my time, for my
//  cached simple test problems, was spent in random number
//  generation.  I was using the compound random number generator
//  from Numerical Recipes.  So speed of a random number
//  generator can be important.
//                   Robert Heckendorn  Oct 20, 1998
// 
//  Revised Dec 31, 1999:  take two seeds instead of one
//                         add a second value for init vars a

//                                                                      //
// // // // // // // // // // // // // // // // // // // // // // // // //

//
//  Next are constants that are specific to *this* random number generator
// and so do not appear in the general purpose rand.h header file which will
// work for all my RNG.  A bit non-standard.
//
// // // // // // // // // // // // // // // // // // // // // // // // //

#define RANDSIZL   (8)  // I recommend 8 for crypto, (4) for less stringent 
#define RANDSIZ    (1<<RANDSIZL)

// size of random number table
#define RNDTABSIZ  997

// warning changed this from 147 to 149 (prime)
#define RNDLAG     211

#define RWORDSIZE   64
#define RHIGHBIT    0x8000000000000000ULL
#define RALLONES    0xffffffffffffffffULL
#define STEP       7

unsigned long long int nextBlockRandom();

// context for the ISAAC random number generator
class Randcontext
{
public:
  unsigned long long int randcnt;
  unsigned long long int randrsl[RANDSIZ];
  unsigned long long int randmem[RANDSIZ];
  unsigned long long int randa;
  unsigned long long int randb;
  unsigned long long int randc;
};



// // // // // // // // // // // // // // // // // // // // // // // // //
//                                                                      //
//                    The Slow Kernel                                   //

//  THIS KERNEL IS BASED ON THE Issac random number generator by
//  Bob Jenkins.  The state space is two arrays of 256 unsigned
//  ints.  It was just of best quality of those tested in the above
//  GA article.

//                                                                      //
// // // // // // // // // // // // // // // // // // // // // // // // //


#define ind(mm, x)  (*(unsigned long long int *)((unsigned char *)(mm) + ((x) & ((RANDSIZ-1)<<2))))
#define rngstep(mix, a, b, mm, m, m2, r, x) \
{ \
      x = *m;  \
      a = (a^(mix)) + *(m2++); \
      *(m++) = y = ind(mm,x) + a + b; \
      *(r++) = b = ind(mm,y>>RANDSIZL) + x; \
      }


static void isaac(Randcontext &context)
{
    unsigned long long int a,b,x,y,*m,*mm,*m2,*r,*mend;

    mm=context.randmem; r=context.randrsl;
    a = context.randa; b = context.randb + (++context.randc);
    for (m = mm, mend = m2 = m+(RANDSIZ/2); m<mend; /*NONE*/) {
	rngstep( a<<13, a, b, mm, m, m2, r, x);
	rngstep( a>>6 , a, b, mm, m, m2, r, x);
	rngstep( a<<2 , a, b, mm, m, m2, r, x);
	rngstep( a>>16, a, b, mm, m, m2, r, x);
    }

    for (m2 = mm; m2<mend; /*NONE*/) {
	rngstep( a<<13, a, b, mm, m, m2, r, x);
	rngstep( a>>6 , a, b, mm, m, m2, r, x);
	rngstep( a<<2 , a, b, mm, m, m2, r, x);
	rngstep( a>>16, a, b, mm, m, m2, r, x);
    }
    context.randb = b; context.randa = a;
}



unsigned long long int slowRand(Randcontext &r)
{
    if (!(r.randcnt)--) {
	isaac(r); 
	r.randcnt = RANDSIZ-1; 
    }

    return r.randrsl[r.randcnt];
}



#define mix(a,b,c,d,e,f,g,h) \
{ \
      a^=b<<23; d+=a; b+=c; \
      b^=c>>4;  e+=b; c+=d; \
      c^=d<<16;  f+=c; d+=e; \
      d^=e>>32; g+=d; e+=f; \
      e^=f<<20; h+=e; f+=g; \
      f^=g>>8;  a+=f; g+=h; \
      g^=h<<16;  b+=g; h+=a; \
      h^=a>>19;  c+=h; a+=b; \
      }



// if seedA is 0 then use the current time and process id
void setSlowRand(Randcontext &context, unsigned long long int &seedA, unsigned long long int &seedB)
{
    unsigned long long int i;
    unsigned long long int a,b,c,d,e,f,g,h;
    unsigned long long int *m,*r;

    context.randa = context.randb = context.randc = 0;
    m=context.randmem;
    r=context.randrsl;

//RH all constants set to different values rather than to same value
//RH Dec 31, 1999
    a=0xe76639d8;    // a random prime number - 1
    b=0xbb40e64d;    // 3141592653 pi
    c=0x90f41b15;    // a random prime number 
    d=0x676593ce;    // a random prime number + 1
    e=0x8c24fb61;    // a random prime number 
    f=0x9e3779b9;    // the golden ratio (float to unsigned long long int?)
    if (seedA==0) {
#ifdef WINDOWS
	seedA = 4294963987ULL;
#else
	seedA = (unsigned long long int)getpid()*4294963987ULL;
#endif
	seedB = ((unsigned long long int)time(NULL)<<1) | 0x1;
    }
    g = seedA;
    h = seedB;

    for (i=0; i<5; ++i) {          // scramble it
	mix(a,b,c,d,e,f,g,h);
    }

//RH the following part of the code was changed to fix that fact that
//RH much of m[] was not being set and i=5.  The fix was to add
//RH the loop.  The fix was approved by Bob Jenkins in email correspondence
//RH Sep 5, 1998
    // fill in m[] with messy stuff 
    for (i=0; i<RANDSIZ; i+=8) {
	mix(a,b,c,d,e,f,g,h);
	m[i  ]=a; m[i+1]=b; m[i+2]=c; m[i+3]=d;
	m[i+4]=e; m[i+5]=f; m[i+6]=g; m[i+7]=h;
    }

    isaac(context);            // fill in the first set of results 
    context.randcnt=RANDSIZ;  // prepare to use the first set of results 

}


// // // // // // // // // // // // // // // // // // // // // // // // //
//                                                                      //
//                    The Fast Kernel                                   //

//  This kernel is based on the r250 Random number generator.  This 
//  kernel must be bootstrapped from another (probably slower) random
//  number generator, in our case, we use the slow kernel.  This means
//  that kernel reinitializing is somewhat of a time sink.
//
//  The objective of this random number generator is to be fast and
//  pretty good for GAs.  see:
//
//  1. Kirkpatrick, S., and E. Stoll, 1981; "A Very Fast
//  Shift-Register Sequence Random Number Generator",
//  Journal of Computational Physics, V.40
//
//
//  2. W.L. Maier, Dr Dobbs Journal May 1991
//
//  The state space is an array of 997 unsigned ints


//  On my machine with compiler optimization on it takes 9.1s to
//  make 10^8 fastRand calls versus 13.9s for the same number of
//  slowRand calls.  This means slowRand is 52% slower.  I have
//  taken the liberty of making the buffer a little larger than
//  r250.  I have run some casual tests for different size
//  buffers and it seems to only make a difference in speed.  Too
//  large and the memory caching suffers.  997 offers about a 2%
//  improvement in speed over 250 and has a nice comfortablely
//  large internal state space.  997 is also prime.

//                                                                      //
// // // // // // // // // // // // // // // // // // // // // // // // //


// these constants are supplied to make the RNG return numbers
// on the open interval [0, 1)
#define PIDIV2	    1.570796326794896619231321691639751442099L
#define PI  	    3.141592653589793238462643383279502884197L

static double unitizer64 = 1.0/18446744073709551616.0;
static double unitizer64_2 = 2.0/18446744073709551616.0;
static double unitizer64_pi = PI/18446744073709551616.0;


// these are the internal state of the RNG
static unsigned long long int r250_buffer[RNDTABSIZ];
unsigned long long int *nextRnd;
unsigned long long int *end;
static unsigned long long int *rnda, *rndb, *start, *startPlusLag, *endMinusLag;


// setRandom: reinitializes the random number generators
// based on a seed value.  If seedA=0 then the time is used and
// the seed values will be returned.
//
static Randcontext slowRandContext;  // allocate this globally to save time

void initRand(unsigned long long int seedA, unsigned long long int seedB)
{
    unsigned long long int *ip, j, k;
    unsigned long long int mask, msb;

    setSlowRand(slowRandContext, seedA, seedB); // seed booting rnd num gen

    start = r250_buffer;
    startPlusLag = r250_buffer+ RNDLAG;
    endMinusLag = r250_buffer + RNDTABSIZ - RNDLAG;
    end = r250_buffer + RNDTABSIZ;

    // fill r250 buffer with RWORDSIZE-1 bit values
    for (ip = start; ip < end; ip++) {
	*ip = slowRand(slowRandContext);        // boot array
    }

    // set some RHIGHBITs to 1
    for (ip = start; ip < end; ip++) {
	if ( slowRand(slowRandContext) & 0x8000) *ip |= RHIGHBIT;
    }

    msb = RHIGHBIT;	        // turn on diagonal bit
    mask = RALLONES;	        // turn off the leftmost bits

    for (j=0; j < RWORDSIZE; j++) {
	k = (STEP * j + 3) % RNDTABSIZ;	// select a word to operate on
	r250_buffer[k] &= mask; // turn off bits left of the diagonal
	r250_buffer[k] |= msb;	// turn on the diagonal bit
	mask >>= 1;
	msb  >>= 1;
    }

    nextRnd = start;
    nextBlockRandom();          // mix up the numbers by refreshing the list
}



// set up fast rand using random system values as seeds.
// this will produced a different sequence every run.
void initRand()                        // WARNING: init required before use!!!
{
    initRand(0ULL, 0ULL);
}




// Generate the next block of random numbers
//
unsigned long long int nextBlockRandom() {
    for (rnda=start, rndb=startPlusLag; rnda<endMinusLag; rnda++, rndb++) {
	*rnda ^= *rndb;
    }
    for (rnda=endMinusLag, rndb=start; rnda<end; rnda++, rndb++) {
	*rnda ^= *rndb;
    }
    nextRnd = start;

    return *nextRnd++;
}


inline unsigned long long int randULL()
{
    return (nextRnd<end) ? *nextRnd++ : nextBlockRandom();
}


// return a uniformly distributed random number between 0 and 1
// NOTE: this can be relatively slow since the multiply is a 64 bit
// multiply.
inline double randUnit()
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


double randGamma(int ia, float scale)
{
    int j;
    double am, e, s, v1, v2, x, y;

    if (ia < 1) {
	printf("ERROR(randGamma): shape parm: %d is less than 1\n", ia);
	exit(1);
    }
    if (ia < 6) {
	x=1.0;
	for (j=1; j<=ia; j++) x *= randUnit();
	x = -log(x);
    } else {
	do {
	    do {
		do {
		    v1=2.0*randUnit()-1.0;
		    v2=2.0*randUnit()-1.0;
		} while (v1*v1+v2*v2 > 1.0);
		y=v2/v1;
		am=ia-1;
		s=sqrt(2.0*am+1.0);
		x=s*y+am;
	    } while (x <= 0.0);
	    e=(1.0+y*y)*exp(am*log(x/am)-s*y);
	} while (randUnit() > e);
    }

    return x/scale;
}
