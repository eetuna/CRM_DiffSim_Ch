#ifndef __CRMMATRIXOPERATIONS_H__
#define __CRMMATRIXOPERATIONS_H__

#include<bits/stdc++.h>

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

// Function to get cofactor of A[p][q] in temp[][]. n is current
// dimension of A[][]
template <int N>
void getCofactor(double A[N][N], double temp[N][N], int p, int q, int n)
{
	int i = 0, j = 0;

	// Looping for each element of the matrix
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			//  Copying into temporary matrix only those element
			//  which are not in given row and column
			if (row != p && col != q)
			{
				temp[i][j++] = A[row][col];

				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
template <int N>
double determinant(double A[N][N], int n)
{
	double D = 0; // Initialize result

	//  Base case : if matrix contains single element
	if (n == 1)
		return A[0][0];

	double temp[N][N]; // To store cofactors

	int sign = 1;  // To store sign multiplier

	// Iterate for each element of first row
	for (int f = 0; f < n; f++)
	{
		// Getting Cofactor of A[0][f]
		getCofactor<N>(A, temp, 0, f, n);
		D += sign * A[0][f] * determinant<N>(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
template <int N>
void adjoint(double A[N][N], double adj[N][N])
{
	if (N == 1)
	{
		adj[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][]
	int sign = 1;
	double temp[N][N];

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			// Get cofactor of A[i][j]
			getCofactor<N>(A, temp, i, j, N);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i+j)%2==0)? 1: -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = (sign)*(determinant<N>(temp, N-1));
		}
	}
}

// Function to calculate and store inverse, returns false if
// matrix is singular
template <int N>
bool inverse(double A[N][N], double inverse_[N][N])
{
	// Find determinant of A[][]
	double det = determinant<N>(A, N);
	if (det == 0)
	{
		std::cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint
	double adj[N][N];
	adjoint<N>(A, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			inverse_[i][j] = adj[i][j]/double(det);

	return true;
}


#endif // __CRMMATRIXOPERATIONS_H__ not defined
