/*
Header file for QuantumMechanics::Eigensystem::HermitianSolver: 

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _EIGENSYSTEM_HERMITIANSOLVER_H_
#define _EIGENSYSTEM_HERMITIANSOLVER_H_

#include <iostream>
#include <functional>
#include <vector>
#include <Eigen/Dense>

namespace QuantumMechanics{

namespace Eigensystem{

class HermitianSolver {
  
private:
	const size_t matrices_count;
	const size_t matrices_size;
	const MatrixXcd * const matrices_array;
	const std::function<MatrixXcd(int)> matrices_function;
	const std::vector<MatrixXcd> * const matrices_vector;

	range computed_range;
	MatrixXd computed_eigenvalues;
	std::vector<MatrixXcd> computed_eigenvectors;

public:
	HermitianSolver() :
		matrices_count(0),
		matrices_size(0),
		matrices_array(nullptr),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		computed_range(range::full()) 
		{ }

	HermitianSolver(const MatrixXcd &M) :
		matrices_count(1),
		matrices_size((M.rows() == M.cols() && M.rows() > 0) ? M.rows() : 0),
		matrices_array(&M),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		computed_range(range::full())
		{ }

	HermitianSolver(size_t n, const MatrixXcd * const M) :
		matrices_count(n),
		matrices_size((M->rows() == M->cols() && M->rows() > 0) ? M->rows() : 0),
		matrices_array(M),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		computed_range(range::full())
		{ }

	HermitianSolver(size_t n, const MatrixXcd * const M, const size_t &size) :
		matrices_count(n),
		matrices_size((M->rows() == M->cols() && M->rows() >= size && size > 0) ? size : 0),
		matrices_array(M),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		computed_range(range::full())
		{ }

	HermitianSolver(size_t n, const std::vector<MatrixXcd> &M) :
		matrices_count(n),
		matrices_size((M[0].rows() == M[0].cols()) ? M[0].rows() : 0),
		matrices_array(nullptr),
		matrices_vector(&M),
		matrices_function(nullptr),
		computed_range(range::full())
		{ }

	HermitianSolver(size_t n, const std::vector<MatrixXcd> &M, const size_t &size) :
		matrices_count(n),
		matrices_size((M[0].rows() == M[0].cols() && M[0].rows() == size) ? size : 0),
		matrices_array(nullptr),
		matrices_vector(&M),
		matrices_function(nullptr),
		computed_range(range::full())
		{ }

	HermitianSolver(size_t n, const std::function<MatrixXcd(int)> &M, const size_t &size) :
		matrices_count(n),
		matrices_size((size > 0) ? size : 0),
		matrices_array(nullptr),
		matrices_vector(nullptr),
		matrices_function(M),
		computed_range(range::full())
		{ }

	virtual ~HermitianSolver() { }

	enum compute_action {
		EigenvaluesOnly,
		EigenvaluesAndVectors
	};

	void compute(compute_action action, range compute_range = range::full())
	{
		if(action == EigenvaluesOnly && (computed_eigenvalues.size() == 0 || computed_range != compute_range))
		{
			computed_range = compute_range;
			compute_eigenvalues();
		}
		else if(action == EigenvaluesAndVectors && (computed_eigenvalues.size() == 0 || computed_range != compute_range))
		{
			computed_range = compute_range;
			compute_eigenvectors();
		}
	}

	MatrixXd eigenvalues();
	MatrixXd eigenvalues(range compute_range);

	static VectorXd eigenvalues(const MatrixXcd &M, range compute_range = range::full());

	std::vector<MatrixXcd> eigenvectors();
	std::vector<MatrixXcd> eigenvectors(range compute_range);

protected:
	void compute_eigenvalues() 
	{
		if (matrices_count <= 0 || matrices_size <= 0 || (!matrices_array && !matrices_vector && !matrices_function))
			return;

		computed_range.fit_indices_to_size(matrices_size);

		const char range_token = (computed_range.value == range::full_range) ? 'A' : ((computed_range.value == range::value_range) ? 'V' : 'I');
		const char job_token = 'N'; // Normal (only eigenvalues).
		const char upper_lower_token = 'U';

		MatrixXd &values = computed_eigenvalues;
		values.resize(matrices_size, matrices_count);

		for (size_t m_index = 0; m_index < matrices_count; m_index++) {

			MatrixXcd M = (matrices_array) ? ((matrices_count == 1) ? *matrices_array : matrices_array[m_index]) :
						((matrices_vector) ? matrices_vector[m_index] : 
							((matrices_function) ? matrices_function(m_index) : MatrixXcd())
						);

			if (M.rows() < matrices_size || M.cols() < matrices_size)
			{
				values.col(m_index).fill(nan(0));
				continue;
			}

			// Problem dimensions:
			const lapack_int major_dim_order = (M.IsRowMajor) ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
			const lapack_int major_dim_length =	(M.IsRowMajor) ? M.rows() : M.cols();
			lapack_int value_count;

			// Additional eigenvector dimensions:
			MKL_Complex16* vectors = nullptr;
			lapack_int lead_dim_vectors = 1; // must be >= 1, or an error occurs. Only used for vector calculation.
			lapack_int* vectors_suppliements = new lapack_int[2*matrices_size];

			lapack_int info = LAPACKE_zheevr(
												major_dim_order,								// param. 0
												job_token,										// param. 1
												range_token,									// param. 2
												upper_lower_token,								// param. 3
												matrices_size,									// param. 4
												M.data(),										// param. 5
												major_dim_length,								// param. 6
												computed_range.lowest_value,					// param. 7
												computed_range.highest_value,					// param. 8
												computed_range.begin_index,						// param. 9
												computed_range.end_index,						// param. 10
												0.,												// param. 11
												&value_count,									// param. 12
												&(values.data()[ m_index * major_dim_length]),	// param. 13
												vectors,										// param. 14
												lead_dim_vectors,								// param. 15
												vectors_suppliements							// param. 16
											);

			delete[] vectors_suppliements;

			if (info != 0) {
				values.col(m_index).fill(nan(0));
			}
		}
	}

	void compute_eigenvectors()
	{
		if (matrices_count <= 0 || matrices_size <= 0 || (!matrices_array && !matrices_vector && !matrices_function))
			return;

		computed_range.fit_indices_to_size(matrices_size);

		const char range_token = (computed_range.value == range::full_range) ? 'A' : ((computed_range.value == range::value_range) ? 'V' : 'I');
		const char job_token = 'V';
		const char upper_lower_token = 'U';


		MatrixXd &values = computed_eigenvalues;
		std::vector<MatrixXcd> &eigenvectors = computed_eigenvectors;

		lapack_int lead_dim_vectors = (range_token == 'I') ? computed_range.end_index - computed_range.begin_index + 1 : matrices_size;

		values.resize(matrices_size, matrices_count);
		eigenvectors.resize(matrices_count, MatrixXcd(matrices_size, lead_dim_vectors));

		for (size_t m_index = 0; m_index < matrices_count; m_index++) {

			MatrixXcd M = (matrices_array) ? ((matrices_count == 1) ? *matrices_array : matrices_array[m_index]) :
						((matrices_vector) ? matrices_vector[m_index] : 
							((matrices_function) ? matrices_function(m_index) : MatrixXcd())
						);

			if (M.rows() < matrices_size || M.cols() < matrices_size)
			{
				values.col(m_index).fill(nan(0));
				continue;
			}

			// Problem dimensions:
			const lapack_int major_dim_order = (M.IsRowMajor) ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
			const lapack_int major_dim_length = (M.IsRowMajor) ? M.rows() : M.cols();
			lapack_int value_count = 0;

			// Additional eigenvector dimensions:
			lapack_int* vectors_suppliements = new lapack_int[2*matrices_size];

			lapack_int info = LAPACKE_zheevr(
												major_dim_order,							// param. 0
												job_token,									// param. 1
												range_token,								// param. 2
												upper_lower_token,							// param. 3
												matrices_size,								// param. 4
												M.data(),									// param. 5
												major_dim_length,							// param. 6
												computed_range.lowest_value,				// param. 7
												computed_range.highest_value,				// param. 8
												computed_range.begin_index,					// param. 9
												computed_range.end_index,					// param. 10
												0.,											// param. 11
												&value_count,								// param. 12
												values.data() + m_index * matrices_size,	// param. 13
												eigenvectors[m_index].data(),				// param. 14
												lead_dim_vectors,							// param. 15
												vectors_suppliements						// param. 16
											);

			delete[] vectors_suppliements;

			if (info != 0) {
				values.col(m_index).fill(nan(0));
			}
		}
	}
};

}

}

#endif                                                      // _EIGENSYSTEM_HERMITIANSOLVER_H_