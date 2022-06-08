#pragma once
/*
    Original FORTRAN77 version by Jorge More, Danny Sorenson, Burton Garbow, Kenneth Hillstrom
        Jorge More, Burton Garbow, Kenneth Hillstrom,
        User Guide for MINPACK-1,
        Technical Report ANL-80-74,
        Argonne National Laboratory, 1980.
        Available at: https://www.osti.gov/biblio/6997568
	The FORTRAN77 version of the code was converted to C++ using fable fortran to c++ converter
		https://cci.lbl.gov/fable/
		and was than manually edited to eliminate superfluous libraries, and to conform to C++ array indexing.
		I have only converted the functions needed to be able to use hybrd1 solver.
	I have used the C++ version of the MINPACK library by John Burkardt for checking correctness.
		https://people.sc.fsu.edu/~jburkardt/cpp_src/minpack/minpack.html
	Note that care should be taken if a parallelization of the algorithm will be pursued due to the use of a common
		work array wa, which is crated once outside the hybrd1 function and reused in all of the subroutines being
		executed inside hybrd1.  This current version is more efficient as it does not require subsequent/repeated
		memory allocations within each subroutine for working memory.  However, it may lead to conflicts if things
		are parallelized within subroutines.  (Note that, the version by Burkhardt allocates some work arrays inside
		the subroutines.)

		------
    Minpack Copyright Notice(1999) University of Chicago.All rights reserved

    Redistributionand use in sourceand binary forms, with or
    without modification, are permitted provided that the
    following conditions are met :

    1. Redistributions of source code must retain the above
    copyright notice, this list of conditionsand the following
    disclaimer.

    2. Redistributions in binary form must reproduce the above
    copyright notice, this list of conditionsand the following
    disclaimer in the documentationand /or other materials
    provided with the distribution.

    3. The end - user documentation included with the
    redistribution, if any, must include the following
    acknowledgment :

    "This product includes software developed by the
    University of Chicago, as Operator of Argonne National
    Laboratory.

    Alternately, this acknowledgment may appear in the software
    itself, ifand wherever such third - party acknowledgments
    normally appear.

    4. WARRANTY DISCLAIMER.THE SOFTWARE IS SUPPLIED "AS IS"
    WITHOUT WARRANTY OF ANY KIND.THE COPYRIGHT HOLDER, THE
    UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
    THEIR EMPLOYEES : (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
    OR NON - INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
    OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
    USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
    THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
    DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
    UNINTERRUPTED, THAT IT IS ERROR - FREE OR THAT ANY ERRORS WILL
    BE CORRECTED.

    5. LIMITATION OF LIABILITY.IN NO EVENT WILL THE COPYRIGHT
    HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
    ENERGY, OR THEIR EMPLOYEES : BE LIABLE FOR ANY INDIRECT,
    INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
    ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
    PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
    SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
    (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
    EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
    POSSIBILITY OF SUCH LOSS OR DAMAGES.
    */

#include <cmath>

//#define MAX(a,b) 	( ((b)>(a))?(b):(a) )
//#define MIN(a,b) 	( ((b)<(a))?(b):(a) )
#define SQR(X) 		( (X)*(X) )
#define FABS(x)		( ((x)>=0)?(x):(-x) )

//C     Function dpmpar
//C
//C     Double precision machine parameters
//C         If the machine has
//C         t base b digits and its smallest and largest exponents are
//C         emin and emax, respectively, then these parameters are
//C
//C         dpmpar(1) = b**(1 - t), the machine precision,
//C
//C         dpmpar(2) = b**(emin - 1), the smallest magnitude,
//C
//C         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
//C
//C     Argonne National Laboratory. MINPACK Project. November 1996.
//C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
#define DPMPAR1 2.220446049250313e-16
//#define DPMPAR2 2.22507385852e-308
//#define DPMPAR3 1.79769313485e+308
#define DPMPAR2 0.4450147717014e-307
#define DPMPAR3 1.0e+30

template <typename adType>
void TrustRegionDogleg(	int n, adType x[], adType fvec[], double tol, int& info, adType wa[], int lwa, NLEqnParams<adType> Params) {
	//C     **********
	//C
	//C     subroutine hybrd1
	//C
	//C     the purpose of hybrd1 is to find a zero of a system of
	//C     n nonlinear functions in n variables by a modification
	//C     of the powell hybrid method. this is done by using the
	//C     more general nonlinear equation solver hybrd. the user
	//C     must provide a subroutine which calculates the functions.
	//C     the jacobian is then calculated by a forward-difference
	//C     approximation.
	//C
	//C     the subroutine statement is
	//C
	//C       subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
	//C
	//C     where
	//C
	//C       fcn is the name of the user-supplied subroutine which    /// fcn() is replaced with hard coded NLEquation()
	//C         calculates the functions. fcn must be declared
	//C         in an external statement in the user calling
	//C         program, and should be written as follows.
	//C
	//C         subroutine fcn(n,x,fvec,iflag)
	//C         integer n,iflag
	//C         double precision x(n),fvec(n)
	//C         ----------
	//C         calculate the functions at x and
	//C         return this vector in fvec.
	//C         ---------
	//C         return
	//C         end
	//C
	//C         the value of iflag should not be changed by fcn unless
	//C         the user wants to terminate execution of hybrd1.
	//C         in this case set iflag to a negative integer.
	//C
	//C       n is a positive integer input variable set to the number
	//C         of functions and variables.
	//C
	//C       x is an array of length n. on input x must contain
	//C         an initial estimate of the solution vector. on output x
	//C         contains the final estimate of the solution vector.
	//C
	//C       fvec is an output array of length n which contains
	//C         the functions evaluated at the output x.
	//C
	//C       tol is a nonnegative input variable. termination occurs
	//C         when the algorithm estimates that the relative error
	//C         between x and the solution is at most tol.
	//C
	//C       info is an integer output variable. if the user has
	//C         terminated execution, info is set to the (negative)
	//C         value of iflag. see description of fcn. otherwise,
	//C         info is set as follows.
	//C
	//C         info = 0   improper input parameters.
	//C
	//C         info = 1   algorithm estimates that the relative error
	//C                    between x and the solution is at most tol.
	//C
	//C         info = 2   number of calls to fcn has reached or exceeded
	//C                    200*(n+1).
	//C
	//C         info = 3   tol is too small. no further improvement in
	//C                    the approximate solution x is possible.
	//C
	//C         info = 4   iteration is not making good progress.
	//C
	//C       wa is a work array of length lwa.
	//C
	//C       lwa is a positive integer input variable not less than
	//C         (n*(3*n+13))/2.
	//C
	//C     subprograms called
	//C
	//C       user-supplied ...... fcn
	//C
	//C       minpack-supplied ... hybrd
	//C
	//C     argonne national laboratory. minpack project. march 1980.
	//C     burton s. garbow, kenneth e. hillstrom, jorge j. more
	//C
	//C     **********
	double factor = 100.0;
	int maxfev = 200 * (n + 1);
	double xtol = tol;
	int ml = n - 1;
	int mu = n - 1;
	double epsfcn = 0.0;
	int mode = 2;
	int j;
	int nprint = 0;
	int lr = (n * (n + 1)) / 2;
	int index = 6 * n + lr;
	int nfev;	
	info = 0;
	//C
	//C     check the input parameters for errors.
	//C
	if (n <= 0 || tol < 0.0 || lwa < (n * (3 * n + 13)) / 2) {
		return;
	}
	//C
	//C     call hybrd.
	//C
	for (j = 0; j < n; j++) {
		wa[j] = 1.0;
	}
	hybrd(n, x, fvec, xtol, maxfev, ml, mu, epsfcn, wa,
		mode, factor, nprint, info, nfev, wa+index, n, wa+6*n,
		lr, wa+n, wa+2*n, wa+3*n, wa+4*n, wa+5*n, Params);
	if (info == 5) {
		info = 4;
	}
}


template <typename adType>
void hybrd(int n, adType x[], adType fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn,
	adType diag[], int mode, double factor, int nprint, int& info, int& nfev,
	adType fjac[], int ldfjac, adType r[], int lr, adType qtf[],
	adType wa1[], adType wa2[], adType wa3[], adType wa4[], NLEqnParams<adType> Params) {
	//C     **********
	//C
	//C     subroutine hybrd
	//C
	//C     the purpose of hybrd is to find a zero of a system of
	//C     n nonlinear functions in n variables by a modification
	//C     of the powell hybrid method. the user must provide a
	//C     subroutine which calculates the functions. the jacobian is
	//C     then calculated by a forward-difference approximation.
	//C
	//C     the subroutine statement is
	//C
	//C       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,
	//C                        diag,mode,factor,nprint,info,nfev,fjac,
	//C                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
	//C
	//C     where
	//C
	//C       fcn is the name of the user-supplied subroutine which
	//C         calculates the functions. fcn must be declared
	//C         in an external statement in the user calling
	//C         program, and should be written as follows.
	//C
	//C         subroutine fcn(n,x,fvec,iflag)
	//C         integer n,iflag
	//C         double precision x(n),fvec(n)
	//C         ----------
	//C         calculate the functions at x and
	//C         return this vector in fvec.
	//C         ---------
	//C         return
	//C         end
	//C
	//C         the value of iflag should not be changed by fcn unless
	//C         the user wants to terminate execution of hybrd.
	//C         in this case set iflag to a negative integer.
	//C
	//C       n is a positive integer input variable set to the number
	//C         of functions and variables.
	//C
	//C       x is an array of length n. on input x must contain
	//C         an initial estimate of the solution vector. on output x
	//C         contains the final estimate of the solution vector.
	//C
	//C       fvec is an output array of length n which contains
	//C         the functions evaluated at the output x.
	//C
	//C       xtol is a nonnegative input variable. termination
	//C         occurs when the relative error between two consecutive
	//C         iterates is at most xtol.
	//C
	//C       maxfev is a positive integer input variable. termination
	//C         occurs when the number of calls to fcn is at least maxfev
	//C         by the end of an iteration.
	//C
	//C       ml is a nonnegative integer input variable which specifies
	//C         the number of subdiagonals within the band of the
	//C         jacobian matrix. if the jacobian is not banded, set
	//C         ml to at least n - 1.
	//C
	//C       mu is a nonnegative integer input variable which specifies
	//C         the number of superdiagonals within the band of the
	//C         jacobian matrix. if the jacobian is not banded, set
	//C         mu to at least n - 1.
	//C
	//C       epsfcn is an input variable used in determining a suitable
	//C         step length for the forward-difference approximation. this
	//C         approximation assumes that the relative errors in the
	//C         functions are of the order of epsfcn. if epsfcn is less
	//C         than the machine precision, it is assumed that the relative
	//C         errors in the functions are of the order of the machine
	//C         precision.
	//C
	//C       diag is an array of length n. if mode = 1 (see
	//C         below), diag is internally set. if mode = 2, diag
	//C         must contain positive entries that serve as
	//C         multiplicative scale factors for the variables.
	//C
	//C       mode is an integer input variable. if mode = 1, the
	//C         variables will be scaled internally. if mode = 2,
	//C         the scaling is specified by the input diag. other
	//C         values of mode are equivalent to mode = 1.
	//C
	//C       factor is a positive input variable used in determining the
	//C         initial step bound. this bound is set to the product of
	//C         factor and the euclidean norm of diag*x if nonzero, or else
	//C         to factor itself. in most cases factor should lie in the
	//C         interval (.1,100.). 100. is a generally recommended value.
	//C
	//C       nprint is an integer input variable that enables controlled
	//C         printing of iterates if it is positive. in this case,
	//C         fcn is called with iflag = 0 at the beginning of the first
	//C         iteration and every nprint iterations thereafter and
	//C         immediately prior to return, with x and fvec available
	//C         for printing. if nprint is not positive, no special calls
	//C         of fcn with iflag = 0 are made.
	//C
	//C       info is an integer output variable. if the user has
	//C         terminated execution, info is set to the (negative)
	//C         value of iflag. see description of fcn. otherwise,
	//C         info is set as follows.
	//C
	//C         info = 0   improper input parameters.
	//C
	//C         info = 1   relative error between two consecutive iterates
	//C                    is at most xtol.
	//C
	//C         info = 2   number of calls to fcn has reached or exceeded
	//C                    maxfev.
	//C
	//C         info = 3   xtol is too small. no further improvement in
	//C                    the approximate solution x is possible.
	//C
	//C         info = 4   iteration is not making good progress, as
	//C                    measured by the improvement from the last
	//C                    five jacobian evaluations.
	//C
	//C         info = 5   iteration is not making good progress, as
	//C                    measured by the improvement from the last
	//C                    ten iterations.
	//C
	//C       nfev is an integer output variable set to the number of
	//C         calls to fcn.
	//C
	//C       fjac is an output n by n array which contains the
	//C         orthogonal matrix q produced by the qr factorization
	//C         of the final approximate jacobian.
	//C
	//C       ldfjac is a positive integer input variable not less than n
	//C         which specifies the leading dimension of the array fjac.
	//C
	//C       r is an output array of length lr which contains the
	//C         upper triangular matrix produced by the qr factorization
	//C         of the final approximate jacobian, stored rowwise.
	//C
	//C       lr is a positive integer input variable not less than
	//C         (n*(n+1))/2.
	//C
	//C       qtf is an output array of length n which contains
	//C         the vector (q transpose)*fvec.
	//C
	//C       wa1, wa2, wa3, and wa4 are work arrays of length n.
	//C
	//C     subprograms called
	//C
	//C       user-supplied ...... fcn
	//C
	//C       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1,
	//C                            qform,qrfac,r1mpyq,r1updt
	//C
	//C       fortran-supplied ... dabs,dmax1,dmin1,min0,mod
	//C
	//C     argonne national laboratory. minpack project. march 1980.
	//C     burton s. garbow, kenneth e. hillstrom, jorge j. more
	//C
	//C     **********
	const double p0001 = 0.0001;
	const double p001 = 0.001;
	const double p1 = 0.1;
	const double p5 = 0.5;
	int iflag = 0;
	int j;
	adType fnorm;
	int msum;
	int iter;
	int ncsuc;
	int ncfail;
	int nslow1;
	int nslow2;
	bool jeval;
	int iwa[1];
	adType xnorm;
	adType delta;
	int i;
	adType sum;
	adType temp;
	bool sing;
	int l;
	adType pnorm;
	adType fnorm1;
	adType actred;
	adType prered;
	adType ratio;
	//C
	//C     epsmch is the machine precision.
	//C
	const double epsmch = DPMPAR1;
	//C
	info = 0;
	nfev = 0;
	//C
	//C     check the input parameters for errors.
	//C
	if (n <= 0 || xtol < 0.0 || maxfev <= 0 || ml < 0 || mu < 0 ||
		factor <= 0.0 || ldfjac < n || lr < (n * (n + 1)) / 2) {
		return;
	}
	if (mode == 2) {
		for (j = 0; j < n; j++) {
			if (diag[j] <= 0.0) {
				info = 0;
				return;
			}
		}
	}
	//C
	//C     evaluate the function at the starting point
	//C     and calculate its norm.
	//C
	iflag = 1;
	NLEquation(x, fvec, Params);
	nfev = 1;
	if (iflag < 0) {
		info = iflag;
		return;
	}
	fnorm = enorm(n, fvec);
	//C
	//C     determine the number of calls to fcn needed to compute
	//C     the jacobian matrix.
	//C
	msum = MIN(ml + mu + 1, n);
	//C
	//C     initialize iteration counter and monitors.
	//C
	iter = 1;
	ncsuc = 0;
	ncfail = 0;
	nslow1 = 0;
	nslow2 = 0;
	//C
	//C     beginning of the outer loop.
	//C
	while (1) {
		jeval = true;
		//C
		//C        calculate the jacobian matrix.
		//C
		iflag = 2;
		fdjac1(n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn, wa1, wa2, Params);
		nfev += msum;
		if (iflag < 0) {
			info = iflag;
			return;
		}
		//C
		//C        compute the qr factorization of the jacobian.
		//C
		qrfac(n, n, fjac, ldfjac, false, iwa, 1, wa1, wa2, wa3);
		//C
		//C        on the first iteration and if mode is 1, scale according
		//C        to the norms of the columns of the initial jacobian.
		//C
		if (iter == 1) {
			if (mode != 2) {
				for (j = 0; j < n; j++) {
					if (wa2[j] == 0.0) {
						diag[j] = 1.0;
					}
					else {
						diag[j] = wa2[j];
					}
				}
			}			//C
			//C        on the first iteration, calculate the norm of the scaled x
			//C        and initialize the step bound delta.
			//C
			for (j = 0; j < n; j++) {
				wa3[j] = diag[j] * x[j];
			}
			xnorm = enorm(n, wa3);
			delta = factor * xnorm;
			if (delta == 0.0) {
				delta = factor;
			}
		}
		//C
		//C        form (q transpose)*fvec and store in qtf.
		//C
		for (i = 0; i < n; i++) {
			qtf[i] = fvec[i];
		}
		for (j = 0; j < n; j++) {
			if (fjac[j + j * ldfjac] != 0.0) {
				sum = 0.0;
				for (i = j; i < n; i++) {
					sum += fjac[i + j * ldfjac] * qtf[i];
				}
				temp = -sum / fjac[j + j * ldfjac];
				for (i = j; i < n; i++) {
					qtf[i] += fjac[i + j * ldfjac] * temp;
				}
			}
		}
		//C
		//C        copy the triangular factor of the qr factorization into r.
		//C
		sing = false;
		for (j = 1; j <= n; j++) {
			l = j;
			for (i = 1; i <= j - 1; i++) {
				r[l - 1] = fjac[(i - 1) + (j - 1) * ldfjac];
				l += n - i;
			}
			r[l - 1] = wa1[j - 1];
			if (wa1[j - 1] == 0.0) {
				sing = true;
			}
		}
		//C
		//C        accumulate the orthogonal factor in fjac.
		//C
		qform(n, n, fjac, ldfjac, wa1);
		//C
		//C        rescale if necessary.
		//C
		if (mode != 2) {
			for (j = 0; j < n; j++) {
				diag[j] = MAX(diag[j], wa2[j]);
			}
		}
		//C
		//C        beginning of the inner loop.
		//C
		while (1) {
			//C
			//C           if requested, call fcn to enable printing of iterates.
			//C
			if (nprint > 0) {
				iflag = 0;
				if (((iter - 1) % nprint) == 0) {
					NLEquation(x, fvec, Params);
				}
				if (iflag < 0) {
					info = iflag;
					return;
				}
			}
			//C
			//C           determine the direction p.
			//C
			dogleg(n, r, lr, diag, qtf, delta, wa1, wa2, wa3);
			//C
			//C           store the direction p and x + p. calculate the norm of p.
			//C
			for (j = 0; j < n; j++) {
				wa1[j] = -wa1[j];
				wa2[j] = x[j] + wa1[j];
				wa3[j] = diag[j] * wa1[j];
			}
			pnorm = enorm(n, wa3);
			//C
			//C           on the first iteration, adjust the initial step bound.
			//C
			if (iter == 1) {
				delta = MIN(delta, pnorm);
			}
			//C
			//C           evaluate the function at x + p and calculate its norm.
			//C
			iflag = 1;
			NLEquation(wa2, wa4, Params);
			nfev++;
			if (iflag < 0) {
				info = iflag;
				return;
			}
			fnorm1 = enorm(n, wa4);
			//C
			//C           compute the scaled actual reduction.
			//C
			actred = -1.0;
			if (fnorm1 < fnorm) {
				actred = 1.0 - SQR(fnorm1 / fnorm);
			}
			//C
			//C           compute the scaled predicted reduction.
			//C
			l = 1;
			for (i = 1; i <= n; i++) {
				sum = 0.0;
				for (j = i; j <= n; j++) {
					sum += r[l - 1] * wa1[j - 1];
					l++;
				}
				wa3[i - 1] = qtf[i - 1] + sum;
			}
			temp = enorm(n, wa3);
			prered = 0.0;
			if (temp < fnorm) {
				prered = 1.0 - SQR(temp / fnorm);
			}
			//C
			//C           compute the ratio of the actual to the predicted
			//C           reduction.
			//C
			ratio = 0.0;
			if (prered > 0.0) {
				ratio = actred / prered;
			}
			//C
			//C           update the step bound.
			//C
			if (ratio < p1) {
				ncsuc = 0;
				ncfail++;
				delta = p5 * delta;
			}
			else {
				ncfail = 0;
				ncsuc++;
				if (ratio >= p5 || ncsuc > 1) {
					delta = MAX(delta, pnorm / p5);
				}
				if (FABS(ratio - 1.0) <= p1) {
					delta = pnorm / p5;
				}
			}
			//C
			//C           test for successful iteration.
			//C
			if (ratio >= p0001) {
				//C
				//C           successful iteration. update x, fvec, and their norms.
				//C
				for (j = 0; j < n; j++) {
					x[j] = wa2[j];
					wa2[j] = diag[j] * x[j];
					fvec[j] = wa4[j];
				}
				xnorm = enorm(n, wa2);
				fnorm = fnorm1;
				iter++;
			}
			//C
			//C           determine the progress of the iteration.
			//C
			nslow1++;
			if (actred >= p001) {
				nslow1 = 0;
			}
			if (jeval) {
				nslow2++;
			}
			if (actred >= p1) {
				nslow2 = 0;
			}
			//C
			//C           test for convergence.
			//C
			if (delta <= xtol * xnorm || fnorm == 0.0) {
				info = 1;
				return;
			}
			//C
			//C           tests for termination and stringent tolerances.
			//C
			if (nfev >= maxfev) {
				info = 2;
				return;
			}
			if (p1 * MAX(p1 * delta, pnorm) <= epsmch * xnorm) {
				info = 3;
				return;
			}
			if (nslow2 == 5) {
				info = 4;
				return;
			}
			if (nslow1 == 10) {
				info = 5;
				return;
			}
			//C
			//C           criterion for recalculating jacobian approximation
			//C           by forward differences.
			//C
			if (ncfail == 2) {
				break;
			}
			//C
			//C           calculate the rank one modification to the jacobian
			//C           and update qtf if necessary.
			//C
			for (j = 0; j < n; j++) {
				sum = 0.0;
				for (i = 0; i < n; i++) {
					sum += fjac[i + j * ldfjac] * wa4[i];
				}
				wa2[j] = (sum - wa3[j]) / pnorm;
				wa1[j] = diag[j] * ((diag[j] * wa1[j]) / pnorm);
				if (ratio >= p0001) {
					qtf[j] = sum;
				}
			}
			//C
			//C           compute the qr factorization of the updated jacobian.
			//C
			r1updt(n, n, r, lr, wa1, wa2, wa3, sing);
			r1mpyq(n, n, fjac, ldfjac, wa2, wa3);
			r1mpyq(1, n, qtf, 1, wa2, wa3);
			//C
			//C           end of the inner loop.
			//C
			jeval = false;
		}
		//C
		//C        end of the outer loop.
		//C
	}
	//C
	//C     termination, either normal or user imposed.
	//C
	if (iflag < 0) {
		info = iflag;
	}
	iflag = 0;
	if (nprint > 0) {
		NLEquation(x, fvec, Params);
	}
	//C
	//C     last card of subroutine hybrd.
	//C
}


template <typename adType>
void dogleg(int n, adType r[], int lr, adType diag[], adType qtb[], adType delta, adType x[],
	adType wa1[], adType wa2[]) {
	//C     **********
	//C
	//C     subroutine dogleg
	//C
	//C     given an m by n matrix a, an n by n nonsingular diagonal
	//C     matrix d, an m-vector b, and a positive number delta, the
	//C     problem is to determine the convex combination x of the
	//C     gauss-newton and scaled gradient directions that minimizes
	//C     (a*x - b) in the least squares sense, subject to the
	//C     restriction that the euclidean norm of d*x be at most delta.
	//C
	//C     this subroutine completes the solution of the problem
	//C     if it is provided with the necessary information from the
	//C     qr factorization of a. that is, if a = q*r, where q has
	//C     orthogonal columns and r is an upper triangular matrix,
	//C     then dogleg expects the full upper triangle of r and
	//C     the first n components of (q transpose)*b.
	//C
	//C     the subroutine statement is
	//C
	//C       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
	//C
	//C     where
	//C
	//C       n is a positive integer input variable set to the order of r.
	//C
	//C       r is an input array of length lr which must contain the upper
	//C         triangular matrix r stored by rows.
	//C
	//C       lr is a positive integer input variable not less than
	//C         (n*(n+1))/2.
	//C
	//C       diag is an input array of length n which must contain the
	//C         diagonal elements of the matrix d.
	//C
	//C       qtb is an input array of length n which must contain the first
	//C         n elements of the vector (q transpose)*b.
	//C
	//C       delta is a positive input variable which specifies an upper
	//C         bound on the euclidean norm of d*x.
	//C
	//C       x is an output array of length n which contains the desired
	//C         convex combination of the gauss-newton direction and the
	//C         scaled gradient direction.
	//C
	//C       wa1 and wa2 are work arrays of length n.
	//C
	//C     subprograms called
	//C
	//C       minpack-supplied ... dpmpar,enorm
	//C
	//C       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
	//C
	//C     argonne national laboratory. minpack project. march 1980.
	//C     burton s. garbow, kenneth e. hillstrom, jorge j. more
	//C
	//C     **********
	//C
	//C     epsmch is the machine precision.
	//C
	const double epsmch = DPMPAR1;
	int i;
	int j;
	int jj;
	int jp1;
	int k;
	int l;
	adType sum;
	adType temp;
	adType qnorm;
	adType gnorm ;
	adType sgnorm;
	adType alpha;
	adType bnorm;
	//C
	//C     first, calculate the gauss-newton direction.
	//C
	jj = (n * (n + 1)) / 2 + 1;
	for (k = 1; k <= n; k++) {
		j = n - k + 1;
		jp1 = j + 1;
		jj = jj - k;
		l = jj + 1;
		sum = 0.0;
		for (i = jp1; i <= n; i++) {
			sum += r[l-1] * x[i-1];
			l++;
		}
		temp = r[jj-1];
		if (temp == 0.0) {
			l = j;
			for (i = 1; i <= j; i++) {
				temp = MAX(temp, FABS(r[l-1]));
				l += n - i;
			}
			temp = epsmch * temp;
			if (temp == 0.0) {
				temp = epsmch;
			}
		}
		x[j-1] = (qtb[j-1] - sum) / temp;
	}
	//C
	//C     test whether the gauss-newton direction is acceptable.
	//C
	for (j = 0; j < n; j++) {
		wa1[j] = 0.0;
		wa2[j] = diag[j] * x[j];
	}
	qnorm = enorm(n, wa2);
	if (qnorm <= delta) {
		return;
	}
	//C
	//C     the gauss-newton direction is not acceptable.
	//C     next, calculate the scaled gradient direction.
	//C
	l = 0;
	for (j = 0; j < n; j++) {
		temp = qtb[j];
		for (i = j; i < n; i++) {
			wa1[i] += r[l] * temp;
			l++;
		}
		wa1[j] = wa1[j] / diag[j];
	}
	//C
	//C     calculate the norm of the scaled gradient and test for
	//C     the special case in which the scaled gradient is zero.
	//C
	gnorm = enorm(n, wa1);
	sgnorm = 0.0;
	alpha = delta / qnorm;
	if (gnorm != 0.0) {
		//C
		//C     calculate the point along the scaled gradient
		//C     at which the quadratic is minimized.
		//C
		for (j = 0; j < n; j++) {
			wa1[j] = (wa1[j] / gnorm) / diag[j];
		}
		l = 1;
		for (j = 0; j < n; j++) {
			sum = 0.0;
			for (i = j; i < n; i++) {
				sum += r[l] * wa1[i];
				l++;
			}
			wa2[j] = sum;
		}
		temp = enorm(n, wa2);
		sgnorm = (gnorm / temp) / temp;
		//C
		//C     test whether the scaled gradient direction is acceptable.
		//C
		alpha = 0.0;
		if (sgnorm < delta) {
			//C
			//C     the scaled gradient direction is not acceptable.
			//C     finally, calculate the point along the dogleg
			//C     at which the quadratic is minimized.
			//C
			bnorm = enorm(n, qtb);
			temp = (bnorm / gnorm) * (bnorm / qnorm) * (sgnorm / delta);
			temp = temp - (delta / qnorm) * SQR(sgnorm / delta) +
				sqrt(SQR(temp - (delta / qnorm)) + (1.0 - SQR(delta / qnorm)) * (1.0 - SQR(sgnorm / delta)));
			alpha = ((delta / qnorm) * (1.0 - SQR(sgnorm / delta))) / temp;
		}
	}
	//C
	//C     form appropriate convex combination of the gauss-newton
	//C     direction and the scaled gradient direction.
	//C
	temp = (1.0 - alpha) * MIN(sgnorm, delta);
	for (j = 0; j < n; j++) {
		x[j] = temp * wa1[j] + alpha * x[j];
	}
	//C
	//C     last card of subroutine dogleg.
	//C
}


template <typename adType>
void fdjac1(int n, adType x[], adType fvec[], adType fjac[], int ldfjac, int& iflag,
	int ml, int mu, double epsfcn, adType wa1[], adType wa2[], NLEqnParams<adType> Params) {
	//C     **********
	//C
	//C     subroutine fdjac1
	//C
	//C     this subroutine computes a forward-difference approximation
	//C     to the n by n jacobian matrix associated with a specified
	//C     problem of n functions in n variables. if the jacobian has
	//C     a banded form, then function evaluations are saved by only
	//C     approximating the nonzero terms.
	//C
	//C     the subroutine statement is
	//C
	//C       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
	//C                         wa1,wa2)
	//C
	//C     where
	//C
	//C       fcn is the name of the user-supplied subroutine which
	//C         calculates the functions. fcn must be declared
	//C         in an external statement in the user calling
	//C         program, and should be written as follows.
	//C
	//C         subroutine fcn(n,x,fvec,iflag)
	//C         integer n,iflag
	//C         double precision x(n),fvec(n)
	//C         ----------
	//C         calculate the functions at x and
	//C         return this vector in fvec.
	//C         ----------
	//C         return
	//C         end
	//C
	//C         the value of iflag should not be changed by fcn unless
	//C         the user wants to terminate execution of fdjac1.
	//C         in this case set iflag to a negative integer.
	//C
	//C       n is a positive integer input variable set to the number
	//C         of functions and variables.
	//C
	//C       x is an input array of length n.
	//C
	//C       fvec is an input array of length n which must contain the
	//C         functions evaluated at x.
	//C
	//C       fjac is an output n by n array which contains the
	//C         approximation to the jacobian matrix evaluated at x.
	//C
	//C       ldfjac is a positive integer input variable not less than n
	//C         which specifies the leading dimension of the array fjac.
	//C
	//C       iflag is an integer variable which can be used to terminate
	//C         the execution of fdjac1. see description of fcn.
	//C
	//C       ml is a nonnegative integer input variable which specifies
	//C         the number of subdiagonals within the band of the
	//C         jacobian matrix. if the jacobian is not banded, set
	//C         ml to at least n - 1.
	//C
	//C       epsfcn is an input variable used in determining a suitable
	//C         step length for the forward-difference approximation. this
	//C         approximation assumes that the relative errors in the
	//C         functions are of the order of epsfcn. if epsfcn is less
	//C         than the machine precision, it is assumed that the relative
	//C         errors in the functions are of the order of the machine
	//C         precision.
	//C
	//C       mu is a nonnegative integer input variable which specifies
	//C         the number of superdiagonals within the band of the
	//C         jacobian matrix. if the jacobian is not banded, set
	//C         mu to at least n - 1.
	//C
	//C       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
	//C         least n, then the jacobian is considered dense, and wa2 is
	//C         not referenced.
	//C
	//C     subprograms called
	//C
	//C       minpack-supplied ... dpmpar
	//C
	//C       fortran-supplied ... dabs,dmax1,dsqrt
	//C
	//C     argonne national laboratory. minpack project. march 1980.
	//C     burton s. garbow, kenneth e. hillstrom, jorge j. more
	//C
	//C     **********
	//C
	//C     epsmch is the machine precision.
	//C
	const double epsmch = DPMPAR1;
	double eps;
	adType temp;
	adType h;
	int msum;
	int j;
	int i;
	int k;
	//C
	eps = sqrt(MAX(epsfcn, epsmch));
	msum = ml + mu + 1;
	if (msum >= n) {
		//C
		//C        computation of dense approximate jacobian.
		//C
		for (j = 0; j < n; j++) {
			temp = x[j];
			h = eps * FABS(temp);
			if (h == 0.0) {
				h = eps;
			}
			x[j] = temp + h;
			NLEquation(x, wa1, Params);
			if (iflag < 0) {
				return;
			}
			x[j] = temp;
			for (i = 0; i < n; i++) {
				fjac[i + j * ldfjac] = (wa1[i] - fvec[i]) / h;
			}
		}
	}
	else {
		//C
		//C        computation of banded approximate jacobian.
		//C
		for (k = 0; k < msum; k++) {
			for (j = k; j < n; j+=msum) {
				wa2[j] = x[j];
				h = eps * FABS(wa2[j]);
				if (h == 0.0) {
					h = eps;
				}
				x[j] = wa2[j] + h;
			}
			NLEquation(x, wa1, Params);
			if (iflag < 0) {
				return;
			}
			for (j = k; j < n; j+=msum) {
				x[j] = wa2[j];
				h = eps * FABS(wa2[j]);
				if (h == 0.0) {
					h = eps;
				}
				for (i = 0; i < n; i++) {
					if (i >= j - mu && i <= j + ml) {
						fjac[i + j * ldfjac] = (wa1[i] - fvec[i]) / h;
					}
					else {
						fjac[i + j * ldfjac] = 0.0;
					}
				}
			}
		}
	}
	//C
	//C     last card of subroutine fdjac1.
	//C
}


template <typename adType>
adType enorm(int n, adType x[])
{
	int i;
	adType return_value = 0.0;

	for (i = 0; i < n; i++)
	{
		return_value += SQR(x[i]);
	}
	return_value = sqrt(return_value);
	return return_value;
}
/*
double enorm(int n, double x[]) {
	//C     **********
	//C
	//C     function enorm
	//C
	//C     given an n-vector x, this function calculates the
	//C     euclidean norm of x.
	//C
	//C     the euclidean norm is computed by accumulating the sum of
	//C     squares in three different sums. the sums of squares for the
	//C     small and large components are scaled so that no overflows
	//C     occur. non-destructive underflows are permitted. underflows
	//C     and overflows do not occur in the computation of the unscaled
	//C     sum of squares for the intermediate components.
	//C     the definitions of small, intermediate and large components
	//C     depend on two constants, rdwarf and rgiant. the main
	//C     restrictions on these constants are that rdwarf**2 not
	//C     underflow and rgiant**2 not overflow. the constants
	//C     given here are suitable for every known computer.
	//C
	//C     the function statement is
	//C
	//C       double precision function enorm(n,x)
	//C
	//C     where
	//C
	//C       n is a positive integer input variable.
	//C
	//C       x is an input array of length n.
	//C
	//C     subprograms called
	//C
	//C       fortran-supplied ... dabs,dsqrt
	//C
	//C     argonne national laboratory. minpack project. march 1980.
	//C     burton s. garbow, kenneth e. hillstrom, jorge j. more
	//C
	//C     **********
	double return_value;
	const double rdwarf = 3.834e-20;
	const double rgiant = 1.304e19;
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 0.0;
	double x1max = 0.0;
	double x3max = 0.0;
	double floatn = n;
	double agiant = rgiant / floatn;
	int i;
	double xabs;
	for (i = 0; i < n; i++) {
		xabs = fabs(x[i]);
		if (xabs > rdwarf && xabs < agiant) {
			goto statement_70;
		}
		if (xabs <= rdwarf) {
			goto statement_30;
		}
		//C
		//C              sum for large components.
		//C
		if (xabs <= x1max) {
			goto statement_10;
		}
		s1 = 1.0 + s1 * SQR(x1max / xabs);
		x1max = xabs;
		goto statement_20;
	statement_10:
		s1 += SQR(xabs / x1max);
	statement_20:
		goto statement_60;
	statement_30:
		//C
		//C              sum for small components.
		//C
		if (xabs <= x3max) {
			goto statement_40;
		}
		s3 = 1.0 + s3 * SQR(x3max / xabs);
		x3max = xabs;
		goto statement_50;
	statement_40:
		if (xabs != 0.0) {
			s3 += SQR(xabs / x3max);
		}
	statement_50:
	statement_60:
		goto statement_80;
	statement_70:
		//C
		//C           sum for intermediate components.
		//C
		s2 += SQR(xabs);
	statement_80:;
	}
	//C
	//C     calculation of norm.
	//C
	if (s1 == 0.0) {
		goto statement_100;
	}
	return_value = x1max * sqrt(s1 + (s2 / x1max) / x1max);
	goto statement_130;
statement_100:
	if (s2 == 0.0) {
		goto statement_110;
	}
	if (s2 >= x3max) {
		return_value = sqrt(s2 * (one + (x3max / s2) * (x3max * s3)));
	}
	if (s2 < x3max) {
		return_value = sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
	}
	goto statement_120;
statement_110:
	return_value = x3max * sqrt(s3);
statement_120:
statement_130:
	return return_value;
	//C
	//C     last card of function enorm.
	//C
}
*/


template <typename adType>
void qform(int m, int n, adType q[], int ldq, adType wa[]) {
	//C     **********
	//C
	//C     subroutine qform
	//C
	//C     this subroutine proceeds from the computed qr factorization of
	//C     an m by n matrix a to accumulate the m by m orthogonal matrix
	//C     q from its factored form.
	//C
	//C     the subroutine statement is
	//C
	//C       subroutine qform(m,n,q,ldq,wa)
	//C
	//C     where
	//C
	//C       m is a positive integer input variable set to the number
	//C         of rows of a and the order of q.
	//C
	//C       n is a positive integer input variable set to the number
	//C         of columns of a.
	//C
	//C       q is an m by m array. on input the full lower trapezoid in
	//C         the first min(m,n) columns of q contains the factored form.
	//C         on output q has been accumulated into a square matrix.
	//C
	//C       ldq is a positive integer input variable not less than m
	//C         which specifies the leading dimension of the array q.
	//C
	//C       wa is a work array of length m.
	//C
	//C     subprograms called
	//C
	//C       fortran-supplied ... min0
	//C
	//C     argonne national laboratory. minpack project. march 1980.
	//C     burton s. garbow, kenneth e. hillstrom, jorge j. more
	//C
	//C     **********
	//C
	//C     zero out upper triangle of q in the first min(m,n) columns.
	//C
	int minmn = MIN(m,n);
	int i;
	int j;
	int k;
	adType sum;
	adType temp;
		for (j = 1; j < minmn; j++) {
			for (i = 0; i < j; i++) {
				q[i + j * ldq] = 0.0;
			}
		}
	//C
	//C     initialize remaining columns to those of the identity matrix.
	//C
	for (j = n; j < m; j++) {
		for (i = 0; i < m; i++) {
			q[i + j * ldq] = 0.0;
		}
		q[j + j * ldq] = 1.0;
	}
	//C
	//C     accumulate q from its factored form.
	//C
	for (k = minmn - 1; k>=0; k--) {
		for (i = k; i < m; i++) {
			wa[i] = q[i + k * ldq];
			q[i + k * ldq] = 0.0;
		}
		q[k + k * ldq] = 1.0;
		if (wa[k] != 0.0) {
			for (j = k; j < m; j++) {
				sum = 0.0;
				for (i = k; i < m; i++) {
					sum += q[i + j * ldq] * wa[i];
				}
				temp = sum / wa[k];
				for (i=k; i<m; i++) {
					q[i + j * ldq] -= temp * wa[i];
				}
			}
		}
	}
	//C
	//C     last card of subroutine qform.
	//C

}


template <typename adType>
void qrfac(int m, int n, adType a[], int lda, bool pivot, int ipvt[], int lipvt,
	adType rdiag[], adType acnorm[], adType wa[]) {
	//C     **********
	//C
	//C     subroutine qrfac
	//C
	//C     this subroutine uses householder transformations with column
	//C     pivoting (optional) to compute a qr factorization of the
	//C     m by n matrix a. that is, qrfac determines an orthogonal
	//C     matrix q, a permutation matrix p, and an upper trapezoidal
	//C     matrix r with diagonal elements of nonincreasing magnitude,
	//C     such that a*p = q*r. the householder transformation for
	//C     column k, k = 1,2,...,min(m,n), is of the form
	//C
	//C                           t
	//C           i - (1/u(k))*u*u
	//C
	//C     where u has zeros in the first k-1 positions. the form of
	//C     this transformation and the method of pivoting first
	//C     appeared in the corresponding linpack subroutine.
	//C
	//C     the subroutine statement is
	//C
	//C       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
	//C
	//C     where
	//C
	//C       m is a positive integer input variable set to the number
	//C         of rows of a.
	//C
	//C       n is a positive integer input variable set to the number
	//C         of columns of a.
	//C
	//C       a is an m by n array. on input a contains the matrix for
	//C         which the qr factorization is to be computed. on output
	//C         the strict upper trapezoidal part of a contains the strict
	//C         upper trapezoidal part of r, and the lower trapezoidal
	//C         part of a contains a factored form of q (the non-trivial
	//C         elements of the u vectors described above).
	//C
	//C       lda is a positive integer input variable not less than m
	//C         which specifies the leading dimension of the array a.
	//C
	//C       pivot is a logical input variable. if pivot is set true,
	//C         then column pivoting is enforced. if pivot is set false,
	//C         then no column pivoting is done.
	//C
	//C       ipvt is an integer output array of length lipvt. ipvt
	//C         defines the permutation matrix p such that a*p = q*r.
	//C         column j of p is column ipvt(j) of the identity matrix.
	//C         if pivot is false, ipvt is not referenced.
	//C
	//C       lipvt is a positive integer input variable. if pivot is false,
	//C         then lipvt may be as small as 1. if pivot is true, then
	//C         lipvt must be at least n.
	//C
	//C       rdiag is an output array of length n which contains the
	//C         diagonal elements of r.
	//C
	//C       acnorm is an output array of length n which contains the
	//C         norms of the corresponding columns of the input matrix a.
	//C         if this information is not needed, then acnorm can coincide
	//C         with rdiag.
	//C
	//C       wa is a work array of length n. if pivot is false, then wa
	//C         can coincide with rdiag.
	//C
	//C     subprograms called
	//C
	//C       minpack-supplied ... dpmpar,enorm
	//C
	//C       fortran-supplied ... dmax1,dsqrt,min0
	//C
	//C     argonne national laboratory. minpack project. march 1980.
	//C     burton s. garbow, kenneth e. hillstrom, jorge j. more
	//C
	//C     **********
	//C
	//C     epsmch is the machine precision.
	//C
	double epsmch = DPMPAR1;
	const double p05 = 0.05;
	adType temp;
	adType ajnorm;
	int i;
	int j;
	int k;
	int kmax;
	int minmn;
	adType sum;
	//C
	//C     compute the initial column norms and initialize several arrays.
	//C
	for (j = 0; j < n; j++) {
		acnorm[j] = enorm(m, a+j*lda);
		rdiag[j] = acnorm[j];
		wa[j] = rdiag[j];
		if (pivot) {
			ipvt[j] = j;
		}
	}
	//C
	//C     reduce a to r with householder transformations.
	//C
	minmn = MIN(m, n);
	for (j = 0; j < minmn; j++) {
		if (pivot) {
			//C
			//C        bring the column of largest norm into the pivot position.
			//C
			kmax = j;
			for (k = j; k < n; k++) {
				if (rdiag[k] > rdiag[kmax]) {
					kmax = k;
				}
			}
			if (kmax != j) {
				for (i = 0; i < m; i++) {
					temp = a[i + j * lda];
					a[i + j * lda] = a[i + kmax * lda];
					a[i + kmax * lda] = temp;
				}
				rdiag[kmax] = rdiag[j];
				wa[kmax] = wa[j];
				k = ipvt[j];
				ipvt[j] = ipvt[kmax];
				ipvt[kmax] = k;
			}
		}
		//C
		//C        compute the householder transformation to reduce the
		//C        j-th column of a to a multiple of the j-th unit vector.
		//C
		ajnorm = enorm(m - j, a + j + j * lda);
		if (ajnorm != 0.0) {
			if (a[j + j * lda] < 0.0) {
				ajnorm = -ajnorm;
			}
			for (i = j; i < m; i++) {
				a[i + j * lda] = a[i + j * lda] / ajnorm;
			}
			a[j + j * lda] += 1.0;
			//C
			//C        apply the transformation to the remaining columns
			//C        and update the norms.
			//C
			for (k = j + 1; k < n; k++) {
				sum = 0.0;
				for (i = j; i < m; i++) {
					sum += a[i + j * lda] * a[i + k * lda];
				}
				temp = sum / a[j + j * lda];
				for (i = j; i < m; i++) {
					a[i + k * lda] -= temp * a[i + j * lda];
				}
				if (pivot && rdiag[k] != 0.0) {
					temp = a[j + k * lda] / rdiag[k];
					rdiag[k] = rdiag[k] * sqrt( MAX(0.0, 1.0 - SQR(temp)) );
					if (p05 * SQR(rdiag[k] / wa[k]) <= epsmch) {
						rdiag[k] = enorm( m-1-j, a+(j+1)+k*lda );
						wa[k] = rdiag[k];
					}
				}
			}
		}
		rdiag[j] = -ajnorm;
	}
	//C
	//C     last card of subroutine qrfac.
	//C
}


template <typename adType>
void r1mpyq(int m, int n, adType a[], int lda, adType v[], adType w[]) {
	//C     **********
	//C
	//C     subroutine r1mpyq
	//C
	//C     given an m by n matrix a, this subroutine computes a*q where
	//C     q is the product of 2*(n - 1) transformations
	//C
	//C           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
	//C
	//C     and gv(i), gw(i) are givens rotations in the (i,n) plane which
	//C     eliminate elements in the i-th and n-th planes, respectively.
	//C     q itself is not given, rather the information to recover the
	//C     gv, gw rotations is supplied.
	//C
	//C     the subroutine statement is
	//C
	//C       subroutine r1mpyq(m,n,a,lda,v,w)
	//C
	//C     where
	//C
	//C       m is a positive integer input variable set to the number
	//C         of rows of a.
	//C
	//C       n is a positive integer input variable set to the number
	//C         of columns of a.
	//C
	//C       a is an m by n array. on input a must contain the matrix
	//C         to be postmultiplied by the orthogonal matrix q
	//C         described above. on output a*q has replaced a.
	//C
	//C       lda is a positive integer input variable not less than m
	//C         which specifies the leading dimension of the array a.
	//C
	//C       v is an input array of length n. v(i) must contain the
	//C         information necessary to recover the givens rotation gv(i)
	//C         described above.
	//C
	//C       w is an input array of length n. w(i) must contain the
	//C         information necessary to recover the givens rotation gw(i)
	//C         described above.
	//C
	//C     subroutines called
	//C
	//C       fortran-supplied ... dabs,dsqrt
	//C
	//C     argonne national laboratory. minpack project. march 1980.
	//C     burton s. garbow, kenneth e. hillstrom, jorge j. more
	//C
	//C     **********
	int i;
	int j;
	adType cos;
	adType sin;
	adType temp;
	//C
	//C     apply the first set of givens rotations to a.
	//C
	for (j = n - 2; j >= 0; j--) {
		if (FABS(v[j]) > 1.0) {
			cos = 1.0 / v[j]; 
			sin = sqrt(1.0 - SQR(cos));
		}
		else {
			sin = v[j];
			cos = sqrt(1.0 - SQR(sin));
		}
		for (i = 0; i < m; i++) {
			temp = cos * a[i + j * lda] - sin * a[i + (n - 1) * lda];
			a[i + (n - 1) * lda] = sin * a[i + j * lda] + cos * a[i + (n - 1) * lda];
			a[i + j * lda] = temp;
		}
	}
	//C
	//C     apply the second set of givens rotations to a.
	//C
	for (j = 0; j < n - 1; j++) {
		if (FABS(w[j]) > 1.0) {
			cos = 1.0 / w[j];
			sin = sqrt(1.0 - SQR(cos));
		}
		else {
			sin = w[j];
			cos = sqrt(1.0 - SQR(sin));
		}
		for (i = 0; i < m; i++) {
			temp = cos * a[i + j * lda] + sin * a[i + (n - 1) * lda];
			a[i + (n - 1) * lda] = -sin * a[i + j * lda] + cos * a[i + (n - 1) * lda];
			a[i + j * lda] = temp;
		}
	}
	//C
	//C     last card of subroutine r1mpyq.
	//C
}


template <typename adType>
void r1updt(int m, int n, adType s[], int ls, adType u[], adType v[], adType w[], bool& sing) {
	//C     **********
	//C
	//C     subroutine r1updt
	//C
	//C     given an m by n lower trapezoidal matrix s, an m-vector u,
	//C     and an n-vector v, the problem is to determine an
	//C     orthogonal matrix q such that
	//C
	//C                   t
	//C           (s + u*v )*q
	//C
	//C     is again lower trapezoidal.
	//C
	//C     this subroutine determines q as the product of 2*(n - 1)
	//C     transformations
	//C
	//C           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
	//C
	//C     where gv(i), gw(i) are givens rotations in the (i,n) plane
	//C     which eliminate elements in the i-th and n-th planes,
	//C     respectively. q itself is not accumulated, rather the
	//C     information to recover the gv, gw rotations is returned.
	//C
	//C     the subroutine statement is
	//C
	//C       subroutine r1updt(m,n,s,ls,u,v,w,sing)
	//C
	//C     where
	//C
	//C       m is a positive integer input variable set to the number
	//C         of rows of s.
	//C
	//C       n is a positive integer input variable set to the number
	//C         of columns of s. n must not exceed m.
	//C
	//C       s is an array of length ls. on input s must contain the lower
	//C         trapezoidal matrix s stored by columns. on output s contains
	//C         the lower trapezoidal matrix produced as described above.
	//C
	//C       ls is a positive integer input variable not less than
	//C         (n*(2*m-n+1))/2.
	//C
	//C       u is an input array of length m which must contain the
	//C         vector u.
	//C
	//C       v is an array of length n. on input v must contain the vector
	//C         v. on output v(i) contains the information necessary to
	//C         recover the givens rotation gv(i) described above.
	//C
	//C       w is an output array of length m. w(i) contains information
	//C         necessary to recover the givens rotation gw(i) described
	//C         above.
	//C
	//C       sing is a logical output variable. sing is set true if any
	//C         of the diagonal elements of the output s are zero. otherwise
	//C         sing is set false.
	//C
	//C     subprograms called
	//C
	//C       minpack-supplied ... dpmpar
	//C
	//C       fortran-supplied ... dabs,dsqrt
	//C
	//C     argonne national laboratory. minpack project. march 1980.
	//C     burton s. garbow, kenneth e. hillstrom, jorge j. more,
	//C     john l. nazareth
	//C
	//C     **********
	//C
	//C     giant is the largest magnitude.
	//C
	const double giant = DPMPAR3;
	const double p25 = 0.25;
	const double p5 = 0.5;
	int jj;
	int l;
	int i;
	int nm1;
	int j;
	adType cotan;
	adType sin;
	adType cos;
	adType tau;
	adType tan;
	adType temp;
	//C
	//C     initialize the diagonal element pointer.
	//C
	jj = (n * (2 * m - n + 1)) / 2 - (m - n);
	//C
	//C     move the nontrivial part of the last column of s into w.
	//C
	l = jj;
	for (i = n; i <= m; i++) {
		w[i-1] = s[l-1];
		l++;
	}
	//C
	//C     rotate the vector v into a multiple of the n-th unit vector
	//C     in such a way that a spike is introduced into w.
	//C
	nm1 = n - 1;
	for (j = n - 1; j >= 1; j--) {
		jj = jj - (m - j + 1);
		w[j-1] = 0.0;
		if (v[j - 1] != 0.0) {
			//C
			//C        determine a givens rotation which eliminates the
			//C        j-th element of v.
			//C
			if (FABS(v[n - 1]) < FABS(v[j - 1])) {
				cotan = v[n - 1] / v[j - 1];
				sin = p5 / sqrt(p25 + p25 * SQR(cotan));
				cos = sin * cotan;
				tau = 1.0;
				if (FABS(cos) * giant > 1.0) {
					tau = 1.0 / cos;
				}
			}
			else {
				tan = v[j-1] / v[n-1];
				cos = p5 / sqrt(p25 + p25 * SQR(tan));
				sin = cos * tan;
				tau = sin;
			}
			//C
			//C        apply the transformation to v and store the information
			//C        necessary to recover the givens rotation.
			//C
			v[n-1] = sin * v[j-1] + cos * v[n-1];
			v[j-1] = tau;
			//C
			//C        apply the transformation to s and extend the spike in w.
			//C
			l = jj;
			for (i = j; i <= m; i++) {
				temp = cos * s[l-1] - sin * w[i-1];
				w[i-1] = sin * s[l-1] + cos * w[i-1];
				s[l-1] = temp;
				l++;
			}
		}
	}
	//C
	//C     add the spike from the rank 1 update to w.
	//C
	for (i = 1; i <= m; i++) {
		w[i-1] += v[n-1] * u[i-1];
	}
	//C
	//C     eliminate the spike.
	//C
	sing = false;
	for (j = 1; j <= nm1; j++) {
		if (w[j - 1] != 0.0) {
			//C
			//C        determine a givens rotation which eliminates the
			//C        j-th element of the spike.
			//C
			if (FABS(s[jj - 1]) < FABS(w[j - 1])) {
				cotan = s[jj-1] / w[j-1];
				sin = p5 / sqrt(p25 + p25 * SQR(cotan));
				cos = sin * cotan;
				tau = 1.0;
				if (FABS(cos) * giant > 1.0) {
					tau = 1.0 / cos;
				}
			}
			else {
				tan = w[j-1] / s[jj-1];
				cos = p5 / sqrt(p25 + p25 * SQR(tan));
				sin = cos * tan;
				tau = sin;
			}
			//C
			//C        apply the transformation to s and reduce the spike in w.
			//C
			l = jj;
			for (i = j; i <= m; i++) {
				temp = cos * s[l-1] + sin * w[i-1];
				w[i-1] = -sin * s[l-1] + cos * w[i-1];
				s[l-1] = temp;
				l++;
			}
			//C
			//C        store the information necessary to recover the
			//C        givens rotation.
			//C
			w[j-1] = tau;
		}
		//C
		//C        test for zero diagonal elements in the output s.
		//C
		if (s[jj-1] == 0.0) {
			sing = true;
		}
		jj += (m - j + 1);
	}
	//C
	//C     move w back into the last column of the output s.
	//C
	l = jj;
	for (i = n; i <= m; i++) {
		s[l-1] = w[i-1];
		l++;
	}
	if (s[jj-1] == 0.0) {
		sing = true;
	}
	//C
	//C     last card of subroutine r1updt.
	//C

}

