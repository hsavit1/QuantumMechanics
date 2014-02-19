/*
 * EigenSystem.cpp
 *
 *  Created on: 05/02/2014
 *      Author: gregersen
 */

#include "EigenSystem.h"
#include <mkl.h>

using namespace Eigen;

namespace QuantumMechanics {

HermitianEigenSolver::HermitianEigenSolver() :
		points(0),
		order(0),
		M_array(nullptr),
		M_function(nullptr)
{ }

HermitianEigenSolver::HermitianEigenSolver(size_t n, MatrixXcd * const M) :
		points(n),
		order((M->rows() == M->cols() && M->rows() > 0) ? M->rows() : 0),
		M_array(M),
		M_function(nullptr)
{ }

HermitianEigenSolver::HermitianEigenSolver(size_t n, MatrixXcd * const M, const size_t &order) :
		points(n),
		order((M->rows() == M->cols() && M->rows() >= order && order > 0) ? order : 0),
		M_array(M),
		M_function(nullptr)
{ }

HermitianEigenSolver::HermitianEigenSolver(size_t n, const std::function<MatrixXcd(int)> &M, const size_t &order) :
		points(n),
		order((order > 0) ? order : 0),
		M_array(nullptr),
		M_function(M)
{ }

HermitianEigenSolver::~HermitianEigenSolver()
{ }

ArrayXXd HermitianEigenSolver::eigenvalues(EigenSystem::range range) const {
	if (points <= 0 || order <= 0 || (!M_array && !M_function))
		return ArrayXXd();

	range.fit_indices_to_order(order);

	const char range_token = (range.value == EigenSystem::range::full_range) ? 'A' :
								((range.value == EigenSystem::range::value_range) ?	'V' : 'I');
	const char job_token = 'N'; // Normal (only eigenvalues).
	const char upper_lower_token = 'U';

	ArrayXXd values(order, points);

	for (size_t p_index = 0; p_index < points; p_index++) {
		MatrixXcd M = (M_array) ? M_array[p_index] :
					((M_function) ? M_function(p_index) : MatrixXcd());

		if (M.rows() < order || M.cols() < order)
		{
			values.col(p_index).fill(nan(0));
			continue;
		}

		// Problem dimensions:
		const lapack_int major_dim_order =
				(M.IsRowMajor) ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
		const lapack_int major_dim_length =
				(M.IsRowMajor) ? M.rows() : M.cols();
		lapack_int value_count = 0;

		// Additional eigenvector dimensions:
		MKL_Complex16* vectors = nullptr;
		lapack_int lead_dim_vectors = 1; // must be >= 1, or an error occurs. Only used for vector calculation.
		lapack_int* vectors_suppliements = nullptr;

		lapack_int info = LAPACKE_zheevr(
											major_dim_order,				// param. 0
											job_token,						// param. 1
											range_token,					// param. 2
											upper_lower_token,				// param. 3
											order,							// param. 4
											M.data(),						// param. 5
											major_dim_length,				// param. 6
											range.lowest_value,				// param. 7
											range.highest_value,			// param. 8
											range.begin_index,				// param. 9
											range.end_index,				// param. 10
											0.,								// param. 11
											&value_count,					// param. 12
											values.data() + p_index * order,// param. 13
											vectors,						// param. 14
											lead_dim_vectors,				// param. 15
											vectors_suppliements			// param. 16
										);

		if (info != 0) {
			values.col(p_index).fill(nan(0));
		}
	}

	return values;
}

ArrayXd HermitianEigenSolver::eigenvalues(MatrixXcd M, EigenSystem::range range)
{

	if (M.rows() != M.cols())
	{
		return ArrayXd();
	}

	const size_t order = M.rows();

	ArrayXd values(order);

	range.fit_indices_to_order(order);

	const char range_token = (range.value == EigenSystem::range::full_range) ? 'A' :
					((range.value == EigenSystem::range::value_range) ? 'V' : 'I');
	const char job_token = 'N'; // Normal (only eigenvalues).
	const char upper_lower_token = 'U';

	// Problem dimensions:
	const lapack_int major_dim_order = (M.IsRowMajor) ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
	const lapack_int major_dim_length = (M.IsRowMajor) ? M.rows() : M.cols();
	lapack_int value_count = 0;

	// Additional eigenvector dimensions:
	MKL_Complex16* vectors = nullptr;
	lapack_int lead_dim_vectors = 1; // must be >= 1, or an error occurs. Only used for vector calculation.
	lapack_int* vectors_suppliements = nullptr;

	lapack_int info = LAPACKE_zheevr(
										major_dim_order,	// param. 0
										job_token,			// param. 1
										range_token,		// param. 2
										upper_lower_token,	// param. 3
										order,				// param. 4
										M.data(),			// param. 5
										major_dim_length,	// param. 6
										range.lowest_value,	// param. 7
										range.highest_value,// param. 8
										range.begin_index,	// param. 9
										range.end_index,	// param. 10
										0.,					// param. 11
										&value_count,		// param. 12
										values.data(),		// param. 13
										vectors,			// param. 14
										lead_dim_vectors,	// param. 15
										vectors_suppliements// param. 16
									);

	if (info != 0) {
		values.fill(nan(0));
	}

	return values;
}

} /* namespace QuantumMechanics */
