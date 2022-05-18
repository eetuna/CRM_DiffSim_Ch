#ifndef __CRMMATRIXOPERATIONS_H__
#define __CRMMATRIXOPERATIONS_H__

//matrix multiplication C=A*B  ( A:D1xD2 B:D2xD3 gives C:D1xD3 )
template <int D1, int D2, int D3, typename Atype, typename Btype>
void mMult_AB (const Atype in_A[D1*D2], const Btype in_B[D2*D3], decltype((*(Atype*)0) * (*(Btype*)0)) out_C[D1*D3]) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D3; j++) {
			decltype((*(Atype*)0) * (*(Btype*)0)) sum=0.0;
			for (int k=0; k<D2; k++) {
				sum +=  in_A[i*D2+k] * in_B[k*D3+j];
			}
			out_C[i*D3+j]=sum;
		}
	}

}


//matrix multiplication C=A^T*B  (A transposed times B)  ( A:D1xD2, A^T:D2xD1, B:D1xD3 gives C:D2xD3 )
template <int D1, int D2, int D3, typename Atype, typename Btype>
void mMult_ATB (const Atype in_A[D1*D2], const Btype in_B[D1*D3], decltype((*(Atype*)0) * (*(Btype*)0)) out_C[D2*D3]) {

	for (int i=0; i<D2; i++) {
		for (int j=0; j<D3; j++) {
			decltype((*(Atype*)0) * (*(Btype*)0)) sum=0.0;
			for (int k=0; k<D1; k++) {
				sum +=  in_A[k*D2+i] * in_B[k*D3+j];
			}
			out_C[i*D3+j]=sum;
		}
	}

}


// matrix scalar multiplication C=s*A
template <int D1, int D2, typename stype, typename Atype>
void mMult_sA (const stype s, const Atype in_A[D1*D2], decltype((*(stype*)0) * (*(Atype*)0)) out_C[D1*D2]) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D2; j++) {
			out_C[i*D2+j]= s *  in_A[i*D2+j];
		}
	}

}


// matrix addition X=A+B, X=A+sB,  or X=A+B+C
template <int D1, int D2, typename Atype, typename Btype>
void mAdd_AB (const Atype in_A[D1*D2], const Btype in_B[D1*D2], decltype((*(Atype*)0) + (*(Btype*)0)) out_X[D1*D2]) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D2; j++) {
			out_X[i*D2+j]= in_A[i*D2+j] + in_B[i*D2+j];
		}
	}

}

template <int D1, int D2, typename Atype, typename stype, typename Btype>
void mAdd_AsB(const Atype in_A[D1 * D2], const stype in_s, const Btype in_B[D1 * D2], decltype((*(Atype*)0) + (*(stype*)0) * (*(Btype*)0)) out_X[D1 * D2]) {

	for (int i = 0; i < D1; i++) {
		for (int j = 0; j < D2; j++) {
			out_X[i * D2 + j] = in_A[i * D2 + j] + in_s * in_B[i * D2 + j];
		}
	}

}

template <int D1, int D2, typename Atype, typename Btype, typename Ctype>
void mAdd_ABC (const Atype in_A[D1*D2], const Btype in_B[D1*D2], const Ctype in_C[D1*D2], decltype((*(Atype*)0) + (*(Btype*)0) + (*(Ctype*)0)) out_X[D1*D2]) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D2; j++) {
			out_X[i*D2+j]= in_A[i*D2+j] + in_B[i*D2+j] + in_C[i*D2+j];
		}
	}

}


// matrix subtraction X=A-B
template <int D1, int D2, typename Atype, typename Btype>
void mSub_AB (const Atype in_A[D1*D2], const Btype in_B[D1*D2], decltype((*(Atype*)0) - (*(Btype*)0)) out_X[D1*D2]) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D2; j++) {
			out_X[i*D2+j]= in_A[i*D2+j] - in_B[i*D2+j];
		}
	}

}


// 2-norm of a vector
template <int D1, typename vtype>
vtype vNormSq (const vtype v[D1]) {

	vtype vn=0.0;

	for (int i=0; i<D1; i++) {
		vn += v[i]*v[i];
	}
	return vn;

}


//copy matrix: B=A  (D1xD2 matrices)
template <int D1, int D2, typename Atype, typename Btype >
void mCopy_AB (const Atype in_A[D1][D2], Btype out_B[D1][D2]) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D2; j++) {
			out_B[i][j]=in_A[i][j];
		}
	}

}


//copy vector: B=A  (D1x1 vector)
template <int D1, typename Atype, typename Btype>
void mCopy_AB (const Atype in_A[D1], Btype out_B[D1]) {

	for (int i=0; i<D1; i++) {
		out_B[i]=in_A[i];
	}

}



#endif // __CRMMATRIXOPERATIONS_H__ not defined
