#pragma once
/*
	Modified from C++ version of the MINPACK library by John Burkardt.
	https://people.sc.fsu.edu/~jburkardt/cpp_src/minpack/minpack.html
	Original FORTRAN77 version by Jorge More, Danny Sorenson, Burton Garbow, Kenneth Hillstrom
		Jorge More, Burton Garbow, Kenneth Hillstrom,
		User Guide for MINPACK-1,
		Technical Report ANL-80-74,
		Argonne National Laboratory, 1980.
		Available at: https://www.osti.gov/biblio/6997568
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


//#include "CRMDYN.hpp"

template <typename adType>
void TrustRegionDogleg(		int n, adType x[], adType fvec[], double tol, int& info, adType wa[], int lwa, NLEqnParams<adType> Params);

template <typename adType>
void hybrd(		int n, adType x[], adType fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn,
						adType diag[], int mode, double factor, int nprint, int& info, int& nfev,
						adType fjac[], int ldfjac, adType r[], int lr, adType qtf[],
						adType wa1[], adType wa2[], adType wa3[], adType wa4[], NLEqnParams<adType> Params);

template <typename adType>
void dogleg(	int n, adType r[], int lr, adType diag[], adType qtb[], adType delta, adType x[],
						adType wa1[], adType wa2[]);

template <typename adType>
void fdjac1(	int n, adType x[], adType fvec[], adType fjac[], int ldfjac, int& iflag, 
						int ml, int mu, double epsfcn, adType wa1[], adType wa2[], NLEqnParams<adType> Params);

template <typename adType>
adType enorm(	int n, adType x[]);

template <typename adType>
void qform(		int m, int n, adType q[], int ldq, adType wa[]);

template <typename adType>
void qrfac(		int m, int n, adType a[], int lda, bool pivot, int ipvt[], int lipvt,
							adType rdiag[], adType acnorm[], adType wa[]);

template <typename adType>
void r1mpyq(	int m, int n, adType a[], int lda, adType v[], adType w[]);

template <typename adType>
void r1updt(	int m, int n, adType s[], int ls, adType u[], adType v[], adType w[], bool& sing);


#include "minpack_Defs.hpp"