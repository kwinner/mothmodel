/*
 * AUTHOR
 *   Dan Sheldon, sheldon@cs.umass.edu
 *   University of Massachusetts, Amherst
 *   December 10, 2012.
 *
 * Standalone version of dbinom from R math library, with lgamma
 * function from lightspeed toolbox. Original dbinom version by:
 *
 *   Catherine Loader, catherine@research.bell-labs.com.
 *   October 23, 2000
 *   ... and the R development team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *
 * DESCRIPTION
 *
 *   To compute the binomial probability, call dbinom(x,n,p).
 *   This checks for argument validity, and calls dbinom_raw().
 *
 *   dbinom_raw() does the actual computation; note this is called by
 *   other functions in addition to dbinom().
 *     (1) dbinom_raw() has both p and q arguments, when one may be represented
 *         more accurately than the other (in particular, in df()).
 *     (2) dbinom_raw() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *         -- but is not the case at all when called e.g., from df() or dbeta() !
 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
 *         Do this in the calling function.
 */

#include <stdio.h>
#include <math.h>
#include "mex.h"

#define ERR(s)      (mexErrMsgTxt(s))

/* These macros are not widely portable at the moment  */
#define ISNAN(x)    (isnan(x))
#define ISFINITE(x) (isfinite(x))
#define NEGINF      (-HUGE_VAL)
//#define INFINITY    (HUGE_VAL)
//#define NAN       (nan())      

/* Macros excerpted from R source code */
#define log_p   give_log
#define EXP(x) (log_p ? (x) : exp(x))
#define ZERO   (log_p ? NEGINF : 0.)
#define ONE    (log_p ? 0. : 1.)

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi)) == log(2*pi)/2 */
#endif

#define NONINT(x) (fabs((x) - floor((x)+0.5)) > 1e-7)
#define FORCEINT(x) (floor((x) + 0.5))

/* gammaln by Tom Minka from lightspeed toolbox

   Logarithm of the gamma function.
   Returns NaN for negative arguments, even though log(gamma(x)) may
   actually be defined.
   
   Written by Tom Minka (lightspeed toolbox)
*/
double gammaln(double x)
{
  static double gamma_series[] = {
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5
  };
  int i;
  double denom, x1, series;
  if(x < 0) return NAN;
  if(x == 0) return INFINITY;
  if(!finite(x)) return x;
  /* Lanczos method */
  denom = x+1;
  x1 = x + 5.5;
  series = 1.000000000190015;
  for(i = 0; i < 6; i++) {
    series += gamma_series[i] / denom;
    denom += 1.0;
  }
  return( M_LN_SQRT_2PI + (x+0.5)*log(x1) - x1 + log(series/x) );
}

/* The following functions are from the R math library */

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 *
 * see also lgammacor() in ./lgammacor.c  which computes almost the same!
 */

double stirlerr(double n)
{

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
    const static double sferr_halves[31] = {
	0.0, /* n=0 - wrong, place holder only */
	0.1534264097200273452913848,  /* 0.5 */
	0.0810614667953272582196702,  /* 1.0 */
	0.0548141210519176538961390,  /* 1.5 */
	0.0413406959554092940938221,  /* 2.0 */
	0.03316287351993628748511048, /* 2.5 */
	0.02767792568499833914878929, /* 3.0 */
	0.02374616365629749597132920, /* 3.5 */
	0.02079067210376509311152277, /* 4.0 */
	0.01848845053267318523077934, /* 4.5 */
	0.01664469118982119216319487, /* 5.0 */
	0.01513497322191737887351255, /* 5.5 */
	0.01387612882307074799874573, /* 6.0 */
	0.01281046524292022692424986, /* 6.5 */
	0.01189670994589177009505572, /* 7.0 */
	0.01110455975820691732662991, /* 7.5 */
	0.010411265261972096497478567, /* 8.0 */
	0.009799416126158803298389475, /* 8.5 */
	0.009255462182712732917728637, /* 9.0 */
	0.008768700134139385462952823, /* 9.5 */
	0.008330563433362871256469318, /* 10.0 */
	0.007934114564314020547248100, /* 10.5 */
	0.007573675487951840794972024, /* 11.0 */
	0.007244554301320383179543912, /* 11.5 */
	0.006942840107209529865664152, /* 12.0 */
	0.006665247032707682442354394, /* 12.5 */
	0.006408994188004207068439631, /* 13.0 */
	0.006171712263039457647532867, /* 13.5 */
	0.005951370112758847735624416, /* 14.0 */
	0.005746216513010115682023589, /* 14.5 */
	0.005554733551962801371038690  /* 15.0 */
    };
    double nn;

    if (n <= 15.0) {
	nn = n + n;
	if (nn == (int)nn) return(sferr_halves[(int)nn]);
	return(lgamma(n + 1.) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI);
    }

    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}

double bd0(double x, double np)
{
    double ej, s, s1, v;
    int j;

    if(!ISFINITE(x) || !ISFINITE(np) || np == 0.0) 
    {
	ERR("Invalid arguments");
    }

    if (fabs(x-np) < 0.1*(x+np)) {
	v = (x-np)/(x+np);
	s = (x-np)*v;/* s using v -- change by MM */
	ej = 2*x*v;
	v = v*v;
	for (j=1; ; j++) { /* Taylor series */
	    ej *= v;
	    s1 = s+ej/((j<<1)+1);
	    if (s1==s) /* last term was effectively 0 */
		return(s1);
	    s = s1;
	}
    }
    /* else:  | x - np |  is not too small */
    return(x*log(x/np)+np-x);
}

double
dbinom_raw(double x, double n, double p, double q, int give_log)
{
    double lf, lc;

    if (p == 0) return((x == 0) ? ONE : ZERO);
    if (q == 0) return((x == n) ? ONE : ZERO);

    if (x == 0) {
	if(n == 0) return ONE;
	lc = (p < 0.1) ? -bd0(n,n*q) - n*p : n*log(q);
	return( EXP(lc) );
    }
    if (x == n) {
	lc = (q < 0.1) ? -bd0(n,n*p) - n*q : n*log(p);
	return( EXP(lc) );
    }
    if (x < 0 || x > n) return( ZERO );

    /* n*p or n*q can underflow to zero if n and p or q are small.  This
       used to occur in dbeta, and gives NaN as from R 2.3.0.  */
    lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x,n*p) - bd0(n-x,n*q);

    /* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
    /* Upto R 2.7.1:
     * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
     * -- following is much better for  x << n : */
    lf = log(M_2PI) + log(x) + log1p(- x/n);

    return EXP(lc - 0.5*lf);
}

double dbinom(double x, double n, double p, int give_log)
{
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(n) || ISNAN(p)) return x + n + p;
#endif

    if (p < 0 || p > 1 || n < 0 || NONINT(n))
	ERR("Invalid parameters");
    
    if (NONINT(x) || x < 0 || !ISFINITE(x))
	ERR("Invalid argument");

    n = FORCEINT(n);
    x = FORCEINT(x);

    return dbinom_raw(x, n, p, 1-p, give_log);
}


/*
 * y = dbinom(x, n, p, as_log, raw)
 *     x is vector, n and p are scalar
 */ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int sizeX, sizeN, sizeP;

    if (nrhs < 3)
	mexErrMsgTxt("Three inputs required");

    // First input: x
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
	mexErrMsgTxt("Input x must be of type real-valued double");

    sizeX = mxGetNumberOfElements(prhs[0]);

    // Second input: n
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
	mexErrMsgTxt("Input n must be of type real-valued double");
    }

    sizeN = mxGetNumberOfElements(prhs[1]);
    
    if ( !(sizeN == sizeX || sizeN == 1) ) {
	mexErrMsgTxt("Input n must be a scalar or have the same size as x");
    }

    // Third input: p
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) {
	mexErrMsgTxt("Input p must be of type real-valued double");
    }

    sizeP = mxGetNumberOfElements(prhs[2]);
    
    if ( !(sizeP == sizeX || sizeP == 1) ) {
	mexErrMsgTxt("Input p must be a scalar or have the same size as x");
    }
    
    double *x = mxGetPr(prhs[0]);
    double *n = mxGetPr(prhs[1]);
    double *p = mxGetPr(prhs[2]);
    bool useLog = false;
    bool raw = false;

    if (nrhs >= 4) useLog = (bool) mxGetScalar(prhs[3]);
    if (nrhs >= 5) raw    = (bool) mxGetScalar(prhs[4]);

    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
    double *y = mxGetPr(plhs[0]);
    
    int i;
    int j; // index for n
    int k; // index for p
    for (i = 0; i < mxGetNumberOfElements(prhs[0]); i++)
    {
        j = (sizeN == sizeX) ? i : 0;  // use same index as x if n is vector; otherwise use scalar value
        k = (sizeP == sizeX) ? i : 0;  // use same index as x if p is vector; otherwise use scalar value

	if (raw)
	    y[i] = dbinom_raw(x[i], n[j], p[k], 1-p[k], useLog);
	else
	    y[i] = dbinom(x[i], n[j], p[k], useLog);
    }
}
