#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpz_poly_q.h"

// #include "mex.h"

/**
 * Compute the Lah number L(n,k), which corresponds to the number of ways
 * a set of n elements can be partitioned into k sets and is given by:
 *     L(n,k) = B(n,k)B(n-1,k-1)(n-k)!
 * where B(n,k) are the binomial coefficients.
 */
void arith_lah(fmpz_t l, slong n, slong k) {
	//first binomial coeff
	fmpz_bin_uiui(l, n, k);
	//fold in the second binomial coeff
	fmpz_t temp; fmpz_init(temp);
	fmpz_bin_uiui(temp, n-1, k-1);
	fmpz_mul(l, l, temp);
	//finally the third factorial
	fmpz_fac_ui(temp, n - k);
	fmpz_mul(l, l, temp);
	//release memory for the temp variable
	fmpz_clear(temp);
}

/**
 * Scales the indeterminate of poly by a
 * Does so by converting the num, den of poly to polynomials with rational coefficients,
 * rescaling them independently, then recomposing both into a single rational function.
 */
void fmpz_poly_q_rescale(fmpz_poly_q_t poly_scaled, const fmpz_poly_q_t poly, const fmpq_t a) {
	fmpz_poly_q_t result; fmpz_poly_q_init(result);

	//split the poly into two polynomials over Q
	fmpq_poly_t num; fmpq_poly_init(num);
	fmpq_poly_t den; fmpq_poly_init(den);
	fmpq_poly_set_fmpz_poly(num, poly->num);
	fmpq_poly_set_fmpz_poly(den, poly->den);

	//rescale each polynomial by a
	fmpq_poly_rescale(num, num, a);
	fmpq_poly_rescale(den, den, a);

	//extract the two components of each of the resulting polynomials
	//num and den are both implemented as a polynomial over Z with an integer denominator,
	//so first we take the two polynomials and directly combine them into the result;
	fmpq_poly_get_numerator(result->num, num);
	fmpq_poly_get_numerator(result->den, den);
	fmpz_poly_q_canonicalise(result);

	//then we multiply by the denominators of num, den
	fmpq_t scalar; fmpq_init(scalar);
	fmpq_set_fmpz_frac(scalar, den->den, num->den);
	mpq_t scalar_mpq; mpq_init(scalar_mpq);
	fmpq_get_mpq(scalar_mpq, scalar);
	fmpz_poly_q_scalar_mul_mpq(result, result, scalar_mpq);

	//copy over into the output structure
	fmpz_poly_q_set(poly_scaled, result); fmpz_poly_q_clear(result);

	fmpq_clear(scalar);
	mpq_clear(scalar_mpq);
}

/**
 * Computes the composition h(s) = (f o g)(s) = f(g(s))
 * Where f and g are both rational functions
 */
void fmpz_poly_q_compose(fmpz_poly_q_t h, const fmpz_poly_q_t f, const fmpz_poly_q_t g) {
	fmpz_poly_q_t result; fmpz_poly_q_init(result);

	//compute the numerator of result (the return value)
	slong D_n = fmpz_poly_degree(fmpz_poly_q_numref(f));
	for(slong i = 0; i <= D_n; i++) {
		//innersum will collect this sum, starting with the ith coefficient of f_num
		fmpz_poly_t innersum; fmpz_poly_init(innersum);
		fmpz_poly_set_si(innersum, fmpz_poly_get_coeff_si(fmpz_poly_q_numref(f), i));

		//multiply by (g->num)^i
		fmpz_poly_t gpower; fmpz_poly_init(gpower);
		fmpz_poly_set(gpower, fmpz_poly_q_numref(g));
		fmpz_poly_pow(gpower, gpower, i);
		fmpz_poly_mul(innersum, innersum, gpower);

		//multiply by (g->den)^(D_n-i)
		fmpz_poly_set(gpower, fmpz_poly_q_denref(g));
		fmpz_poly_pow(gpower, gpower, D_n - i);
		fmpz_poly_mul(innersum, innersum, gpower);

		//add to the numerator of result
		fmpz_poly_add(fmpz_poly_q_numref(result), fmpz_poly_q_numref(result), innersum);

		fmpz_poly_clear(innersum);
		fmpz_poly_clear(gpower);
	}

	//compute the denominator of result
	slong D_d = fmpz_poly_degree(fmpz_poly_q_denref(f));
	for(slong j = 0; j <= D_d; j++) {
		//innersum will collect this sum, starting with the jth coefficient of f_den
		fmpz_poly_t innersum; fmpz_poly_init(innersum);
		fmpz_poly_set_si(innersum, fmpz_poly_get_coeff_si(fmpz_poly_q_denref(f), j));

		//multiply by (g->num)^j
		fmpz_poly_t gpower; fmpz_poly_init(gpower);
		fmpz_poly_set(gpower, fmpz_poly_q_numref(g));
		fmpz_poly_pow(gpower, gpower, j);
		fmpz_poly_mul(innersum, innersum, gpower);

		//multiply by (g->den)^(D_d-j)
		fmpz_poly_clear(gpower); fmpz_poly_init(gpower);
		fmpz_poly_set(gpower, fmpz_poly_q_denref(g));
		fmpz_poly_pow(gpower, gpower, D_d - j);
		fmpz_poly_mul(innersum, innersum, gpower);

		//add to the denominator of result
		fmpz_poly_add(fmpz_poly_q_denref(result), fmpz_poly_q_denref(result), innersum);

		fmpz_poly_clear(innersum);
		fmpz_poly_clear(gpower);
	}

	//multiply by (g->den)^(D_d-D_n)
	if(D_d >= 0 && D_n >= 0) {
		if(D_d > D_n) {
			fmpz_poly_t gpower; fmpz_poly_init(gpower);
			fmpz_poly_set(gpower, fmpz_poly_q_denref(g));
			fmpz_poly_pow(gpower, gpower, D_d - D_n);
			fmpz_poly_mul(fmpz_poly_q_numref(result), fmpz_poly_q_numref(result), gpower);

			fmpz_poly_clear(gpower);
		}
		if(D_d < D_n) {
			fmpz_poly_t gpower; fmpz_poly_init(gpower);
			fmpz_poly_set(gpower, fmpz_poly_q_denref(g));
			fmpz_poly_pow(gpower, gpower, D_n - D_d); //invert, then multiply into the denominator
			fmpz_poly_mul(fmpz_poly_q_denref(result), fmpz_poly_q_denref(result), gpower);

			fmpz_poly_clear(gpower);
		}
	}

	//finally, canonicalise to correct the numerator and denominator
	fmpz_poly_q_canonicalise(result);

	fmpz_poly_q_set(h, result);
	fmpz_poly_q_clear(result);
}

// void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
void evidence_pgfba(fmpz_poly_q_t gamma, fmpq_t aprime, const fmpz_poly_q_t f, const fmpq_t a, const ulong y, const fmpq_t rho) {
	fmpz_poly_q_t gamma_result; fmpz_poly_q_init(gamma_result);
	fmpq_t aprime_result; fmpq_init(aprime_result);

	if(fmpq_is_zero(a)) {
		//if a is zero, then the exponential part doesn't actually exist, and the evidence function is simplified accordingly
		//first, take the yth derivative of f
		fmpz_poly_q_set(gamma_result, f);
		for(int ideriv = 0; ideriv < y; ideriv++) {
			fmpz_poly_q_derivative(gamma_result, gamma_result);
		}

		//then compose it with s(1-rho), which is just rescaling by (1-rho)
		fmpq_t minusrho; fmpq_init(minusrho);
		fmpq_one(minusrho);
		fmpq_sub(minusrho, minusrho, rho);
		fmpz_poly_q_rescale(gamma_result, gamma_result, minusrho);

		//multiply by (sp)^y/y!
		fmpz_poly_q_t multiplicand; fmpz_poly_q_init(multiplicand);
		fmpz_poly_set_coeff_si(multiplicand->num, 1, rho->num);
		fmpz_poly_set_coeff_si(multiplicand->den, 0, rho->den);
		fmpz_poly_q_pow(multiplicand, multiplicand, y);
		mpz_t yfactorial; mpz_init(yfactorial);
		mpz_fac_ui(yfactorial, y);
		fmpz_poly_q_scalar_div_mpz(multiplicand, multiplicand, yfactorial);
		fmpz_poly_q_mul(gamma_result, multiplicand, gamma_result);

		//copy over to outputs
		fmpz_poly_q_set(gamma, gamma_result);
		fmpq_set(aprime, a);
	} else {
		for(ulong k = 0; k <= y; k++) {
			//use innersum to collect the sum of all terms of the innermost sum
			fmpz_poly_q_t innersum; fmpz_poly_q_init(innersum);
			for(ulong z = 1; z <= k; z++) {
				//compute the numerator of this term
				//start with lah number
				fmpz_t lah; fmpz_init(lah);
				arith_lah(lah, k, z);

				//and a raised to the z power
				fmpz_t a2z; fmpz_init(a2z);
				fmpz_pow_ui(a2z, &a->num, z);

				//compute the coefficient for the numerator
				fmpz_t num_coeff; fmpz_init(num_coeff);
				fmpz_mul(num_coeff, lah, a2z);

				fmpz_poly_t numerator; fmpz_poly_init(numerator);
				fmpz_poly_set_coeff_fmpz(numerator, 0, num_coeff);

				//compute the denominator of this term
				fmpz_t den_coeff; fmpz_init(den_coeff);
				fmpz_pow_ui(den_coeff, &a->den, z);

				fmpz_poly_t denominator; fmpz_poly_init(denominator);
				fmpz_poly_set_coeff_fmpz(denominator, k+z, den_coeff);

				//combine numerator and denominator, canonicalise and then fold into innersum
				fmpz_poly_q_t term; fmpz_poly_q_init(term);
				fmpz_poly_set(fmpz_poly_q_numref(term), numerator);
				fmpz_poly_set(fmpz_poly_q_denref(term), denominator);
				fmpz_poly_q_canonicalise(term);
				fmpz_poly_q_add(innersum, innersum, term);

				//clean up
				fmpz_clear(lah);
				fmpz_clear(a2z);
				fmpz_clear(num_coeff);
				fmpz_poly_clear(numerator);
				fmpz_clear(den_coeff);
				fmpz_poly_clear(denominator);
				fmpz_poly_q_clear(term);
			}

			//take the y-kth derivative of f
			fmpz_poly_q_t fderiv; fmpz_poly_q_init(fderiv);
			fmpz_poly_q_set(fderiv, f);

			for(int ideriv = 0; ideriv < y-k; ideriv++) {
				fmpz_poly_q_derivative(fderiv, fderiv);
			}

			//multiply by (-1)^k, is -1 if k mod 2 == 1
			if(k % 2 == 1)
				fmpz_poly_q_scalar_mul_si(fderiv, fderiv, -1);

			//multiply by y choose k
			mpz_t ychoosek; mpz_init(ychoosek);
			mpz_bin_uiui(ychoosek, y, k);
			fmpz_poly_q_scalar_mul_mpz(fderiv, fderiv, ychoosek);

			//multiply fderiv with innersum
			fmpz_poly_q_mul(fderiv, fderiv, innersum);

			//fold fderiv into gamma_result
			fmpz_poly_q_add(gamma_result, gamma_result, fderiv);

			fmpz_poly_q_clear(fderiv);
			fmpz_poly_q_clear(innersum);
		}

		//take the composition F^(y)(s(1-p))
		//first update g(s) = as^{-1}+b by updating aprime_result = a(1-p)
		fmpq_one(aprime_result);
		fmpq_sub(aprime_result, aprime_result, rho);
		fmpq_mul(aprime_result, a, aprime_result);

		//compose f(s) (which is currently in gamma_result)
		//with s(1-p), which is basically rescaling by 1-p
		fmpq_t minusrho; fmpq_init(minusrho);
		fmpq_one(minusrho);
		fmpq_sub(minusrho, minusrho, rho);
		fmpz_poly_q_rescale(gamma_result, gamma_result, minusrho);

		//multiply by (sp)^y/y!
		fmpz_poly_q_t multiplicand; fmpz_poly_q_init(multiplicand);
		fmpz_poly_set_coeff_si(multiplicand->num, 1, rho->num);
		fmpz_poly_set_coeff_si(multiplicand->den, 0, rho->den);
		fmpz_poly_q_pow(multiplicand, multiplicand, y);
		mpz_t yfactorial; mpz_init(yfactorial);
		mpz_fac_ui(yfactorial, y);
		fmpz_poly_q_scalar_div_mpz(multiplicand, multiplicand, yfactorial);
		fmpz_poly_q_mul(gamma_result, multiplicand, gamma_result);

		//copy out the results
		fmpz_poly_q_set(gamma, gamma_result); fmpz_poly_q_clear(gamma_result);
		fmpq_set(aprime, aprime_result);      fmpq_clear(aprime_result);

		fmpq_clear(minusrho);
		fmpz_poly_q_clear(multiplicand);
		mpz_clear(yfactorial);
	}
}

void transition_pgfba(fmpq_t aprime, fmpq_t bprime, const fmpq_t a, const fmpq_t b, const fmpq_t lambda) {
	fmpq_add(aprime, a, lambda);
	fmpq_sub(bprime, b, lambda);
}

void condition_pgfba(fmpz_poly_q_t fprime, fmpq_t aprime, fmpq_t bprime, const fmpz_poly_q_t f, const fmpq_t a, const fmpq_t b, const fmpq_t rho) {
	fmpz_poly_q_t fprime_result; fmpz_poly_q_init(fprime_result);
	fmpq_t        aprime_result; fmpq_init(       aprime_result);
	fmpq_t        bprime_result; fmpq_init(       bprime_result);

	// we need the difference between the numerator and den of rho to construct g
	fmpz_t rho_diff; fmpz_init(rho_diff);
	fmpz_sub(rho_diff, fmpq_numref(rho), fmpq_denref(rho));

	fmpz_poly_q_t g; fmpz_poly_q_init(g);
	fmpz_poly_set_coeff_fmpz(fmpz_poly_q_numref(g), 1, fmpq_numref(rho));
	fmpz_poly_set_coeff_fmpz(fmpz_poly_q_denref(g), 1, rho_diff);
	fmpz_poly_set_coeff_fmpz(fmpz_poly_q_denref(g), 0, fmpq_denref(rho));
	fmpz_poly_q_canonicalise(g);

	fmpz_poly_q_compose(fprime_result, f, g);

	//multiply by 1/((rho-1)*s + 1)
	//but since rho is rational, this actually becomes beta/(alpha-beta)s+beta
	fmpz_poly_q_clear(g); fmpz_poly_q_init(g);
	fmpz_poly_set_coeff_fmpz(fmpz_poly_q_numref(g), 0, fmpq_denref(rho));
	fmpz_poly_set_coeff_fmpz(fmpz_poly_q_denref(g), 0, fmpq_denref(rho));
	fmpz_poly_set_coeff_fmpz(fmpz_poly_q_denref(g), 1, rho_diff);

	fmpz_poly_q_mul(fprime_result, fprime_result, g);

	//update the terms of the exponential, aprime and bprime
	//aprime = a/rho
	fmpq_div(aprime_result, a, rho);
	//bprime = a-a/rho+b
	fmpq_sub(bprime_result, a, aprime_result);
	fmpq_add(bprime_result, bprime_result, b);

	//copy over to output
	fmpz_poly_q_set(fprime, fprime_result); fmpz_poly_q_clear(fprime_result);
	fmpq_set(aprime, aprime_result); fmpq_clear(aprime_result);
	fmpq_set(bprime, bprime_result); fmpq_clear(bprime_result);

	fmpz_clear(rho_diff);
	fmpz_poly_q_clear(g);
}

void pgfba(fmpz_poly_q_t beta, fmpq_t a, fmpq_t b, const ulong* y, const fmpq* rho, const fmpq* lambda, const ulong K) {
	//1. beta(s) := (1-s)^-1
	fmpz_poly_q_set_str(beta, "1  1/2  1 -1");
	fmpq_zero(a);
	fmpq_zero(b);

	//2. gamma(s) := evidence(beta(s), y_K, rho_K)
	evidence_pgfba(beta, a, beta, a, y[K-1], rho + (K - 1));

	//3. for k = k-1 down to 1
	for(int k = K-2; k >= 0; k--) {
		transition_pgfba(a, b, a, b, lambda + k);
		condition_pgfba(beta, a, b, beta, a, b, rho + k);
		evidence_pgfba(beta, a, beta, b, y[k], rho + k);
	}
}

int main(int argc, const char * argv[]) {
    //configure parameters
    // ulong K    = 5;
    ulong K = 3;
    ulong y[3] = {1, 2, 3};

    fmpq* rho = _fmpq_vec_init(K);
    fmpq_set_si(rho + 0, 1, 2);
    fmpq_set_si(rho + 1, 1, 3);
    fmpq_set_si(rho + 2, 1, 4);
    // fmpq_set_si(rho + 3, 1, 5);
    // fmpq_set_si(rho + 4, 1, 6);

    fmpq* lambda = _fmpq_vec_init(K-1);
    fmpq_set_si(lambda + 0, 1, 2);
    fmpq_set_si(lambda + 1, 2, 3);
    // fmpq_set_si(lambda + 2, 3, 4);
    // fmpq_set_si(lambda + 3, 4, 5);

    //prepare outputs
    fmpq_t a;           fmpq_init(a);
    fmpq_t b;           fmpq_init(b);
    fmpz_poly_q_t beta; fmpz_poly_q_init(beta);

    pgfba(beta, a, b, y, rho, lambda, K);

    flint_printf("beta:\t"); fmpz_poly_q_print_pretty(beta, "s"); flint_printf("\n");
    flint_printf("a:\t"); fmpq_print(a); flint_printf("\n");
    flint_printf("b:\t"); fmpq_print(b); flint_printf("\n");

    _fmpq_vec_clear(rho, K);
    _fmpq_vec_clear(lambda, K-1);
    fmpq_clear(a);
    fmpq_clear(b);
    fmpz_poly_q_clear(beta);

    return 0;
}