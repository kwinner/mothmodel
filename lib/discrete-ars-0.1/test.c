#include <math.h>
#include <strings.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define MATHLIB_STANDALONE
#include "Rmath.h"

#include "ars.h"

/****************************************/
/* Binomial distribution */
/****************************************/
typedef struct {
    int n;
    double p;
} BinomialParams;

double binomial_logpmf(double k, void *data)
{
    BinomialParams *params = (BinomialParams *) data;
    return dbinom(k, params->n, params->p, 1);
}

/****************************************/
/* Poisson   */
/****************************************/
double poisson_logpmf(double k, void *data)
{
    double *lambda = (double *) data;
    return dpois(k, *lambda, 1);
}

/****************************************/
/* gsl_urand   */
/****************************************/
gsl_rng *rng;

void gsl_urand_setup()
{
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    rng = gsl_rng_alloc (T);
}

double gsl_urand()
{
    return gsl_rng_uniform( rng );
}

void gsl_urand_destroy()
{
    gsl_rng_free(rng);
}

/****************************************/
/* test_binomial
/****************************************/
void test_binomial() 
{
    // Create Binomial distribution
    BinomialParams params;
    params.n = 1000;
    params.p = 0.87;
    

    int M = 1e5;
    double *samples = (double *) malloc(M*sizeof(double));

    double startpoints[] = {2, 3};
    int nstart = 2;

    discrete_ars(samples, M, &gsl_urand, &binomial_logpmf, &params, 0, params.n, startpoints, nstart);

    double mean = 0;
    int i;
    for (i = 0; i < M; i++)
    {
	mean += ((double) samples[i])/M;
    }
    printf("mean = %.4f, sample mean = %.4f\n", params.n*params.p, mean);
    free(samples);
}

/*************************************************************/
/* test_poisson
/*************************************************************/
void test_poisson() 
{
    double lambda = 900;

    int M = 1e5;

    double *samples = (double *) malloc(M*sizeof(double));

    double mode = find_mode(&poisson_logpmf, &lambda, 950, 960);
    printf("Mode is %.0f\n", mode);

    double startpoints[] = {2, 901};
    int nstart = 2;
    
    discrete_ars(samples, M, &gsl_urand, &poisson_logpmf, &lambda, 0, HUGE_VAL, startpoints, nstart);

    double mean = 0;
    for (int i = 0; i < M; i++)
    {
	mean += ((double) samples[i])/M;
    }
    printf("mean = %.4f, sample mean = %.4f\n", lambda, mean);
    free(samples);
}


/*************************************************************/
/* main
/*************************************************************/
int main (void)
{
    gsl_urand_setup();
    
    printf("Binomial test\n");
    test_binomial();

    printf("Poisson test\n");
    test_poisson();

    gsl_urand_destroy();
    
    return 0;    
}
