#pragma once

//matrix multiplication C=A*B  ( A:D1xD2 B:D2xD3 gives C:D1xD3 )
template <int D1, int D2, int D3, typename Atype, typename Btype, typename Ctype>
void mMult_AB (const Atype& in_A, const Btype& in_B, Ctype out_C) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D3; j++) {
			decltype(in_A[0]*in_B[0]) sum=0.0;
			for (int k=0; k<D2; k++) {
				sum +=  in_A[i*D2+k] * in_B[k*D3+j];
			}
			out_C[i*D3+j]=sum;
		}
	}

}


//matrix multiplication C=A^T*B  (A transposed times B)  ( A:D1xD2, A^T:D2xD1, B:D1xD3 gives C:D2xD3 )
template <int D1, int D2, int D3, typename Atype, typename Btype, typename Ctype>
void mMult_ATB (const Atype& in_A, const Btype& in_B, Ctype& out_C) {

	for (int i=0; i<D2; i++) {
		for (int j=0; j<D3; j++) {
			decltype(in_A[0] * in_B[0]) sum=0.0;
			for (int k=0; k<D1; k++) {
				sum +=  in_A[k*D2+i] * in_B[k*D3+j];
			}
			out_C[i*D3+j]=sum;
		}
	}

}


//matrix multiplication C=A*B^T  (A times B transposed)  ( A:D1xD2, B:D3xD2, B^T:D2xD3,  gives C:D1xD3 )
template <int D1, int D2, int D3, typename Atype, typename Btype, typename Ctype>
void mMult_ABT(const Atype& in_A, const Btype& in_B, Ctype& out_C) {

	for (int i = 0; i < D1; i++) {
		for (int j = 0; j < D3; j++) {
			decltype(in_A[0] * in_B[0]) sum = 0.0;
			for (int k = 0; k < D2; k++) {
				sum += in_A[i * D2 + k] * in_B[j * D2 + k];
			}
			out_C[i * D3 + j] = sum;
		}
	}

}


// matrix scalar multiplication C=s*A  ( s:scalar, A,C:D1xD2 )
template <int D1, int D2, typename stype, typename Atype, typename Ctype>
void mMult_sA (const stype s, const Atype& in_A, Ctype& out_C) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D2; j++) {
			out_C[i*D2+j]= s *  in_A[i*D2+j];
		}
	}

}


//matrix multiply and add C=C+A*B  ( A:D1xD2 B:D2xD3 gives C:D1xD3 )
template <int D1, int D2, int D3, typename Atype, typename Btype, typename Ctype>
void mMultAdd_AB(const Atype& in_A, const Btype& in_B, Ctype& out_C) {

	for (int i = 0; i < D1; i++) {
		for (int j = 0; j < D3; j++) {
			decltype(in_A[0] * in_B[0]) sum = 0.0;
			for (int k = 0; k < D2; k++) {
				sum += in_A[i * D2 + k] * in_B[k * D3 + j];
			}
			out_C[i * D3 + j] += sum;
		}
	}

}


//matrix multiply and subtract C=C-A*B  ( A:D1xD2 B:D2xD3 gives C:D1xD3 )
template <int D1, int D2, int D3, typename Atype, typename Btype, typename Ctype>
void mMultSub_AB(const Atype& in_A, const Btype& in_B, Ctype& out_C) {

	for (int i = 0; i < D1; i++) {
		for (int j = 0; j < D3; j++) {
			decltype(in_A[0] * in_B[0]) sum = 0.0;
			for (int k = 0; k < D2; k++) {
				sum += in_A[i * D2 + k] * in_B[k * D3 + j];
			}
			out_C[i * D3 + j] -= sum;
		}
	}

}


// matrix addition X=A+B, A=A+B, X=A+sB,  or X=A+B+C   ( A,B,C,X:D1xD2, s:scalar )
template <int D1, int D2, typename Atype, typename Btype, typename Xtype>
void mAdd_AB (const Atype& in_A, const Btype& in_B, Xtype& out_X) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D2; j++) {
			out_X[i*D2+j]= in_A[i*D2+j] + in_B[i*D2+j];
		}
	}

}

template <int D1, int D2, typename Atype, typename Btype>
void mAdd_AB(Atype& inout_A, const Btype& in_B) {

	for (int i = 0; i < D1; i++) {
		for (int j = 0; j < D2; j++) {
			inout_A[i * D2 + j] += in_B[i * D2 + j];
		}
	}

}

template <int D1, int D2, typename Atype, typename stype, typename Btype, typename Xtype>
void mAdd_AsB(const Atype& in_A, const stype in_s, const Btype& in_B, Xtype& out_X) {

	for (int i = 0; i < D1; i++) {
		for (int j = 0; j < D2; j++) {
			out_X[i * D2 + j] = in_A[i * D2 + j] + in_s * in_B[i * D2 + j];
		}
	}

}

template <int D1, int D2, typename Atype, typename Btype, typename Ctype, typename Xtype>
void mAdd_ABC (const Atype& in_A, const Btype& in_B, const Ctype& in_C, Xtype& out_X) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D2; j++) {
			out_X[i*D2+j]= in_A[i*D2+j] + in_B[i*D2+j] + in_C[i*D2+j];
		}
	}

}


// matrix subtraction X=A-B   ( A,B,X:D1xD2 )
template <int D1, int D2, typename Atype, typename Btype, typename Xtype>
void mSub_AB (const Atype& in_A, const Btype& in_B, Xtype out_X) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D2; j++) {
			out_X[i*D2+j]= in_A[i*D2+j] - in_B[i*D2+j];
		}
	}

}


// matrix subtraction A=A-B   ( A,B:D1xD2 )
template <int D1, int D2, typename Atype, typename Btype>
void mSub_AB(Atype& inout_A, const Btype& in_B) {

	for (int i = 0; i < D1; i++) {
		for (int j = 0; j < D2; j++) {
			inout_A[i * D2 + j] -= in_B[i * D2 + j];
		}
	}

}


// 2-norm of a vector  ( D1x1 vector )
template <int D1, typename vtype>
auto vNormSq (const vtype& v) {

	decltype(v[0]*v[0]) vn=0.0;

	for (int i=0; i<D1; i++) {
		vn += v[i]*v[i];
	}
	return vn;

}


//copy matrix: B=A  (D1xD2 matrices)
template <int D1, int D2, typename Atype, typename Btype >
void mCopy_ABm (const Atype in_A[D1][D2], Btype out_B[D1][D2]) {

	for (int i=0; i<D1; i++) {
		for (int j=0; j<D2; j++) {
			out_B[i][j]=in_A[i][j];
		}
	}

}



//copy vector: B=A  (D1x1 vector)
template <int D1, typename Atype, typename Btype>
void mCopy_AB (const Atype& in_A, Btype out_B) {

	for (int i=0; i<D1; i++) {
		out_B[i]=in_A[i];
	}

}

