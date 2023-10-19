#include "strassen_mm.hpp"
#include <stdexcept>


// function to add matrix that is called by strassen_mm
template <typename T>

std::vector<std::vector<T>> add_matrix(std::vector<std::vector<T>> A, std::vector<std::vector<T>> B, int split_index, int multiplier = 1)
{
	for (int i = 0; i < split_index; i++)
		for (int j = 0; j < split_index; j++)
			A[i][j] = A[i][j] + (multiplier * B[i][j]);
	return A;
}


// core function that computes the multiplication of the matrices
template <typename T>

std::vector<std::vector<T>> strassen_mm(const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B)
{
	int col_1 = A[0].size();
	int row_1 = A.size();
	int col_2 = B[0].size();
	int row_2 = B.size();

	if (col_1 != row_2) {
		std::cout << "\nError: The number of columns in Matrix "
				"A must be equal to the number of rows in "
				"Matrix B\n";
		return {};
	}

	std::vector<T> result_matrix_row(col_2, 0);
	std::vector<std::vector<T> > result_matrix(row_1, result_matrix_row);

	
    // this is the base case for the recursive function
    if (col_1 == 1)
		result_matrix[0][0]	= A[0][0] * B[0][0];
	else {
		int split_index = col_1 / 2; // split the matrix and determine size of splits

		std::vector<T> row_vector(split_index, 0);

		// initialize the matrix elements for the operation
		std::vector<std::vector<T>> a00(split_index, row_vector);
		std::vector<std::vector<T>> a01(split_index, row_vector);
		std::vector<std::vector<T>> a10(split_index, row_vector);
		std::vector<std::vector<T>> a11(split_index, row_vector);
		std::vector<std::vector<T>> b00(split_index, row_vector);
		std::vector<std::vector<T>> b01(split_index, row_vector);
		std::vector<std::vector<T>> b10(split_index, row_vector);
		std::vector<std::vector<T>> b11(split_index, row_vector);


		// assign values to a11.. a2, b11...b22 from matrix A and B
		for (int i = 0; i < split_index; i++)
			for (int j = 0; j < split_index; j++) {
				a00[i][j] = A[i][j];
				a01[i][j] = A[i][j + split_index];
				a10[i][j] = A[split_index + i][j];
				a11[i][j] = A[i + split_index][j + split_index];
				b00[i][j] = B[i][j];
				b01[i][j] = B[i][j + split_index];
				b10[i][j] = B[split_index + i][j];
				b11[i][j] = B[i + split_index][j + split_index];
			}

		// calculate P1...P7 as q, r, s, t...
		std::vector<std::vector<T>> p(strassen_mm(a00, add_matrix(b01, b11, split_index, -1)));
		std::vector<std::vector<T>> q(strassen_mm(add_matrix(a00, a01, split_index), b11));
		std::vector<std::vector<T>> r(strassen_mm(add_matrix(a10, a11, split_index), b00));
		std::vector<std::vector<T>> s(strassen_mm(a11, add_matrix(b10, b00, split_index, -1)));
		std::vector<std::vector<T>> t(strassen_mm(add_matrix(a00, a11, split_index), add_matrix(b00, b11, split_index)));
		std::vector<std::vector<T>> u(strassen_mm(add_matrix(a01, a11, split_index, -1), add_matrix(b10, b11, split_index)));
		std::vector<std::vector<T>> v(strassen_mm( add_matrix(a00, a10, split_index, -1), add_matrix(b00, b01, split_index)));

		// result matrix "C"
		std::vector<std::vector<T>> result_matrix_00(add_matrix(add_matrix(add_matrix(t, s, split_index), u, split_index),	q, split_index, -1));
		std::vector<std::vector<T>> result_matrix_01(add_matrix(p, q, split_index));
		std::vector<std::vector<T>> result_matrix_10(add_matrix(r, s, split_index));
		std::vector<std::vector<T>> result_matrix_11(add_matrix(add_matrix(add_matrix(t, p, split_index), r, split_index, -1), v, split_index, -1));

		for (int i = 0; i < split_index; i++)
			for (int j = 0; j < split_index; j++) {
				result_matrix[i][j]	= result_matrix_00[i][j];
				result_matrix[i][j + split_index]= result_matrix_01[i][j];
				result_matrix[split_index + i][j]= result_matrix_10[i][j];
				result_matrix[i + split_index][j + split_index]	= result_matrix_11[i][j];
			}

		a00.clear();
		a01.clear();
		a10.clear();
		a11.clear();
		b00.clear();
		b01.clear();
		b10.clear();
		b11.clear();
		p.clear();
		q.clear();
		r.clear();
		s.clear();
		t.clear();
		u.clear();
		v.clear();
		result_matrix_00.clear();
		result_matrix_01.clear();
		result_matrix_10.clear();
		result_matrix_11.clear();
	}


	return result_matrix;
}



template std::vector<std::vector<double>> strassen_mm <double> (const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B);

