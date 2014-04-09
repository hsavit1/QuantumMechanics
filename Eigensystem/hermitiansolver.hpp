/*
Header file for QuantumMechanics::Eigensystem::HermitianSolver: 

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the 

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _EIGENSYSTEM_HERMITIANSOLVER_H_
#define _EIGENSYSTEM_HERMITIANSOLVER_H_

#include <iostream> // For printing useful messages.
#include <fstream> // For empty streams.

// Including two different storage-forms of matrices.
#include <functional>
#include <vector>

#include <Math/Eigen/Dense> 	// Includes math in namespace Eigen.
using namespace Eigen; 	// Escaping namespace Eigen!

#include <mkl.h>	// Includes lapack.
#include <tbb/mutex.h>	// Includes lapack.

#include "range.hpp" 	// Includes the custom struct range.

namespace QuantumMechanics {

namespace Eigensystem {


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
	
	std::function<void(double)> progress_function;

public:
	HermitianSolver() :
		matrices_count(0),
		matrices_size(0),
		matrices_array(nullptr),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		computed_range(range::full()),
		progress_function(nullptr)
		{ }

	HermitianSolver(const MatrixXcd &M) :
		matrices_count(1),
		matrices_size((M.rows() == M.cols() && M.rows() > 0) ? M.rows() : 0),
		matrices_array(&M),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		computed_range(range::full()),
		progress_function(nullptr)
		{ }

	HermitianSolver(size_t n, const MatrixXcd * const M) :
		matrices_count(n),
		matrices_size((M->rows() == M->cols() && M->rows() > 0) ? M->rows() : 0),
		matrices_array(M),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		computed_range(range::full()),
		progress_function(nullptr)
		{ }

	HermitianSolver(size_t n, const MatrixXcd * const M, const size_t &size) :
		matrices_count(n),
		matrices_size((M->rows() >= size && M->cols() >= size && size > 0) ? size : 0),
		matrices_array(M),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		computed_range(range::full()),
		progress_function(nullptr)
		{ }

	HermitianSolver(size_t n, const std::vector<MatrixXcd> &M) :
		matrices_count(n),
		matrices_size((M[0].rows() == M[0].cols()) ? M[0].rows() : 0),
		matrices_array(nullptr),
		matrices_vector(&M),
		matrices_function(nullptr),
		computed_range(range::full()),
		progress_function(nullptr)
		{ }

	HermitianSolver(size_t n, const std::vector<MatrixXcd> &M, const size_t &size) :
		matrices_count(n),
		matrices_size((M[0].rows() == M[0].cols() && M[0].rows() == size) ? size : 0),
		matrices_array(nullptr),
		matrices_vector(&M),
		matrices_function(nullptr),
		computed_range(range::full()),
		progress_function(nullptr)
		{ }

	HermitianSolver(size_t n, const std::function<MatrixXcd(int)> &M, const size_t &size) :
		matrices_count(n),
		matrices_size((size > 0) ? size : 0),
		matrices_array(nullptr),
		matrices_vector(nullptr),
		matrices_function(M),
		computed_range(range::full()),
		progress_function(nullptr)
		{ }

	virtual ~HermitianSolver() { }
	
	// TODO: to disable all log from these elements set this to false.	
	static bool logging_enabled;
	
	void enableProgressFeedback(std::function<void(double)> feedback_function) {
		progress_function = feedback_function;
	}

protected:
	static std::ostream null_stream;
	
	std::ostream & log() 
	{ 
		if(logging_enabled)
			return (std::clog << "Eigensystem::HermitianSolver message: ");
		else
			return null_stream;
	}
	
	std::ostream & logAppend() 
	{ 
		if(logging_enabled)
			return std::clog;
		else
			return null_stream;
	}
	
	void compute_eigenvalues() 
	{
		log() << "Initiating attempt to compute eigenvalues (only) in";
		if(computed_range.type_value == range::full_range)
			logAppend() << " full range." << std::endl;
		else if(computed_range.type_value == range::index_range)
			logAppend() << " index range " << computed_range.begin_index << " to " << computed_range.end_index << "." << std::endl;
		else if(computed_range.type_value == range::value_range)
			logAppend() << " value range " << computed_range.lowest_value << " to " << computed_range.highest_value << "." << std::endl;
		else if(computed_range.type_value == range::mid_index_range)
			logAppend() << " middle range " << computed_range.begin_index << " to " << computed_range.end_index << "." << std::endl;
		
		if (matrices_count <= 0 || matrices_size <= 0 || (!matrices_array && !matrices_vector && !matrices_function))
		{
			log() << "Failed to compute eigenvalues due to ";
			if(matrices_count <= 0)
				logAppend() << "matrices_count = " << matrices_count << " ";
			if(matrices_size <= 0)
				logAppend() << "matrices_size = " << matrices_size << " ";
			if((!matrices_array && !matrices_vector && !matrices_function))
				logAppend() << "matrices_array = matrices_vector = matrices_function = null";
			logAppend() << std::endl;
			return;
		}

		computed_range.fit_indices_to_size(matrices_size);

		const char range_token = (computed_range.type_value == range::full_range) ? 'A' : ((computed_range.type_value == range::value_range) ? 'V' : 'I');
		const char job_token = 'N'; // Normal (only eigenvalues).
		const char upper_lower_token = 'U';

		MatrixXd &values = computed_eigenvalues;
		values.resize(matrices_size, matrices_count);
		
		log() << "Computing eigenvalues (only) for " << matrices_count << " matrices in";
		if(computed_range.type_value == range::full_range)
			logAppend() << " full range." << std::endl;
		else if(computed_range.type_value == range::index_range)
			logAppend() << " index range " << computed_range.begin_index << " to " << computed_range.end_index << "." << std::endl;
		else if(computed_range.type_value == range::value_range)
			logAppend() << " value range " << computed_range.lowest_value << " to " << computed_range.highest_value << "." << std::endl;
		
		lapack_int value_count_max = 0;
		
		const double progress_step((progress_function) ? 1. / double(matrices_count) : 0.);
		double current_progress(0.);
		tbb::mutex progress_mutex; 
		
		if(progress_function)
			progress_function(current_progress);

		_Cilk_for (size_t m_index = 0; m_index < matrices_count; m_index++) {

			MatrixXcd M = (matrices_array) ? ((matrices_count == 1) ? *matrices_array : matrices_array[m_index]) :
						((matrices_vector) ? (*matrices_vector)[m_index] : 
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
												computed_range.begin_index+1,					// param. 9
												computed_range.end_index+1,						// param. 10
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
			else if (value_count_max < value_count)
				value_count_max = value_count;
			
			if(progress_function)
			{
				progress_mutex.lock();
				current_progress += progress_step;
				progress_function(current_progress);
				progress_mutex.unlock();
			}
		}
		
		if(progress_function && current_progress < 1.)
			progress_function(1.);
			
		log() << "Done computing eigenvalues." << std::endl;
		
		
		if(value_count_max < matrices_size)
		{
			log() << "The results are trimmed to a maximum eigenvalue count of " << value_count_max << "." << std::endl;
			values.conservativeResize(value_count_max, matrices_count);
		}
	}

	void compute_eigenvectors()
	{
		log() << "Initiating attempt to compute eigenvalues and eigenvectors in";
		if(computed_range.type_value == range::full_range)
			logAppend() << " full range." << std::endl;
		else if(computed_range.type_value == range::index_range)
			logAppend() << " index range " << computed_range.begin_index << " to " << computed_range.end_index << "." << std::endl;
		else if(computed_range.type_value == range::value_range)
			logAppend() << " value range " << computed_range.lowest_value << " to " << computed_range.highest_value << "." << std::endl;
		else if(computed_range.type_value == range::mid_index_range)
			logAppend() << " middle range " << computed_range.begin_index << " to " << computed_range.end_index << "." << std::endl;
		
		if (matrices_count <= 0 || matrices_size <= 0 || (!matrices_array && !matrices_vector && !matrices_function))
		{
			log() << "Failed to compute eigenvalues and eigenvectors due to ";
			if(matrices_count <= 0)
				logAppend() << "matrices_count = " << matrices_count << " ";
			if(matrices_size <= 0)
				logAppend() << "matrices_size = " << matrices_size << " ";
			if((!matrices_array && !matrices_vector && !matrices_function))
				logAppend() << "matrices_array = matrices_vector = matrices_function = null";
			logAppend() << std::endl;
			return;
		}

		computed_range.fit_indices_to_size(matrices_size);

		const char range_token = (computed_range.type_value == range::full_range) ? 'A' : ((computed_range.type_value == range::value_range) ? 'V' : 'I');
		const char job_token = 'V';
		const char upper_lower_token = 'U';


		MatrixXd &values = computed_eigenvalues;
		std::vector<MatrixXcd> &eigenvectors = computed_eigenvectors;

		const lapack_int lead_dim_vectors = matrices_size;
		const lapack_int size_vectors = (range_token == 'I') ? computed_range.end_index - computed_range.begin_index + 1 : matrices_size;
		
		values = MatrixXd(matrices_size, matrices_count);
		eigenvectors = std::vector<MatrixXcd>(matrices_count, MatrixXcd(lead_dim_vectors, size_vectors));
		
		log() << "Computing eigenvalues and eigenvectors for " << matrices_count << " matrices in";
		if(computed_range.type_value == range::full_range)
			logAppend() << " full range ("<< range_token <<")." << std::endl;
		else if(computed_range.type_value == range::index_range)
			logAppend() << " index range " << computed_range.begin_index << " to " << computed_range.end_index << " ("<< range_token <<")." << std::endl;
		else if(computed_range.type_value == range::value_range)
			logAppend() << " value range " << computed_range.lowest_value << " to " << computed_range.highest_value << " ("<< range_token <<")." << std::endl;
		
		lapack_int value_count_max = 0;
		
		_Cilk_for (size_t m_index = 0; m_index < matrices_count; m_index++) {

			MatrixXcd M = (matrices_array) ? ((matrices_count == 1) ? *matrices_array : matrices_array[m_index]) :
						((matrices_vector) ? (*matrices_vector)[m_index] : 
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
												computed_range.begin_index+1,				// param. 9
												computed_range.end_index+1,					// param. 10
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
			else if (value_count_max < value_count)
				value_count_max = value_count;
			
			if(value_count < matrices_size)
				eigenvectors[m_index].conservativeResize(matrices_size, value_count);
		}
		
		log() << "Done computing eigenvalues and eigenvectors." << std::endl;
		
		if(value_count_max < matrices_size)
		{
			log() << "The results are trimmed to a maximum eigenvalue count of " << value_count_max << "." << std::endl;
			values.conservativeResize(value_count_max, matrices_count);
		}
			
	}
	
public:
	
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
		else if(action == EigenvaluesAndVectors && (computed_eigenvectors.size() == 0 || computed_range != compute_range))
		{
			computed_range = compute_range;
			compute_eigenvectors();
		}
	}

	MatrixXd eigenvalues()
	{
		if(computed_eigenvalues.size() == 0)
		{
			compute_eigenvalues();
		}

		return computed_eigenvalues;
	}
	
	MatrixXd eigenvalues(range compute_range)
	{
		if(computed_eigenvalues.size() == 0 || computed_range != compute_range)
		{
			computed_range = compute_range;
			compute_eigenvalues();
		}

		return computed_eigenvalues;
	}

	static VectorXd eigenvalues(const MatrixXcd &M, range compute_range = range::full())
	{
		return HermitianSolver(M).eigenvalues(compute_range);
	}

	std::vector<MatrixXcd> eigenvectors()
	{
		if(computed_eigenvectors.size() == 0)
		{
			compute_eigenvectors();
		}

		return computed_eigenvectors;
	}
	
	std::vector<MatrixXcd> eigenvectors(range compute_range)
	{
		if(computed_eigenvectors.size() == 0 || computed_range != compute_range)
		{
			computed_range = compute_range;
			compute_eigenvectors();
		}

		return computed_eigenvectors;
	}
};


std::ostream HermitianSolver::null_stream(0);
bool HermitianSolver::logging_enabled(false);

}

}

#endif                                                      // _EIGENSYSTEM_HERMITIANSOLVER_H_