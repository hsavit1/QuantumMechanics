/*
Header file for QuantumMechanics::GreensFormalism::SemiInfiniteChain: 

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the 

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _GREENSFORMALISM_SEMIINFINITECHAIN_H_
#define _GREENSFORMALISM_SEMIINFINITECHAIN_H_

#include <iostream> // For printing useful messages.
#include <fstream> // For empty streams.

// Including two different storage-forms of matrices.
#include <functional>
#include <vector>

#include <Math/Dense>

#include <tbb/mutex.h>	// Includes a thread lock mutex.

namespace QuantumMechanics {

namespace GreensFormalism {
	
class SemiInfiniteChain {
	
private:
	const size_t matrices_count;
	const size_t matrices_size;
	const MatrixXcd * const inverse_matrices_array;
	const std::function<MatrixXcd(int)> inverse_matrices_function;
	const std::vector<MatrixXcd> * const inverse_matrices_vector;
	
	Array3i block_sizes;

	std::vector<MatrixXcd> computed_matrices;
	
	std::function<void(double)> progress_function;
	
public:
	SemiInfiniteChain() :
		matrices_count(0),
		matrices_size(0),
		inverse_matrices_array(nullptr),
		inverse_matrices_vector(nullptr),
		inverse_matrices_function(nullptr),
		block_sizes(Array::Zero<3>()),
		progress_function(nullptr)
		{ }

	SemiInfiniteChain(const MatrixXcd &M) :
		matrices_count(1),
		matrices_size((M.rows() == M.cols() && M.rows() > 0) ? M.rows() : 0),
		inverse_matrices_array(&M),
		inverse_matrices_vector(nullptr),
		inverse_matrices_function(nullptr),
		block_sizes(Array::Zero<3>()),
		computed_matrices(std::vector<MatrixXcd>(1, MatrixXcd())),
		progress_function(nullptr)
		{ }

	SemiInfiniteChain(size_t n, const MatrixXcd * const M) :
		matrices_count(n),
		matrices_size((M->rows() == M->cols() && M->rows() > 0) ? M->rows() : 0),
		inverse_matrices_array(M),
		inverse_matrices_vector(nullptr),
		inverse_matrices_function(nullptr),
		block_sizes(Array::Zero<3>()),
		computed_matrices(std::vector<MatrixXcd>(n, MatrixXcd())),
		progress_function(nullptr)
		{ }

	SemiInfiniteChain(size_t n, const MatrixXcd * const M, const size_t &size) :
		matrices_count(n),
		matrices_size((M->rows() >= size && M->cols() >= size && size > 0) ? size : 0),
		inverse_matrices_array(M),
		inverse_matrices_vector(nullptr),
		inverse_matrices_function(nullptr),
		block_sizes(Array::Zero<3>()),
		computed_matrices(std::vector<MatrixXcd>(n, MatrixXcd())),
		progress_function(nullptr)
		{ }

	SemiInfiniteChain(size_t n, const std::vector<MatrixXcd> &M) :
		matrices_count(n),
		matrices_size((M[0].rows() == M[0].cols()) ? M[0].rows() : 0),
		inverse_matrices_array(nullptr),
		inverse_matrices_vector(&M),
		inverse_matrices_function(nullptr),
		block_sizes(Array::Zero<3>()),
		computed_matrices(std::vector<MatrixXcd>(n, MatrixXcd())),
		progress_function(nullptr)
		{ }

	SemiInfiniteChain(size_t n, const std::vector<MatrixXcd> &M, const size_t &size) :
		matrices_count(n),
		matrices_size((M[0].rows() == M[0].cols() && M[0].rows() == size) ? size : 0),
		inverse_matrices_array(nullptr),
		inverse_matrices_vector(&M),
		inverse_matrices_function(nullptr),
		block_sizes(Array::Zero<3>()),
		computed_matrices(std::vector<MatrixXcd>(n, MatrixXcd())),
		progress_function(nullptr)
		{ }

	SemiInfiniteChain(size_t n, const std::function<MatrixXcd(int)> &M, const size_t &size) :
		matrices_count(n),
		matrices_size((size > 0) ? size : 0),
		inverse_matrices_array(nullptr),
		inverse_matrices_vector(nullptr),
		inverse_matrices_function(M),
		block_sizes(Array::Zero<3>()),
		computed_matrices(std::vector<MatrixXcd>(n, MatrixXcd())),
		progress_function(nullptr)
		{ }

	virtual ~SemiInfiniteChain() { }
	
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
			return (std::clog << "GreensFormalism::SemiInfiniteChain message: ");
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
	
	void compute_matrix_from_left() 
	{
		computed_matrices = std::vector<MatrixXcd>(matrices_count, MatrixXcd(matrices_size, matrices_size));
		
		_Cilk_for (long m = 0; m < matrices_count; m++) {
			compute_matrix_from_left(m);
		}
	}
	
	void compute_matrix_from_right() 
	{
		computed_matrices = std::vector<MatrixXcd>(matrices_count, MatrixXcd(matrices_size, matrices_size));
		
		_Cilk_for (long m = 0; m < matrices_count; m++) {
			compute_matrix_from_right(m);
		}
	}
	
	void compute_matrix_from_left(const long &matrix_index) 
	{
		if((block_sizes == 0).any())
		{
			log() << "The sizes given do not fit in the matrices. Cannot guess the value."
			return;
		}
		
		ArrayXi offsets(4);
		offsets[0] = 0;
		// A cumulative sum (first,last,destination_start).
		std::partial_sum(&block_sizes[0], &block_sizes[2], &offsets[1]);
		
		MatrixXcd M = (matrices_array) ? 
						(
							(matrices_count == 1) ? 
								*matrices_array : 
								matrices_array[matrix_index]
						) :
						(
							(matrices_vector) ? 
								(*matrices_vector)[matrix_index] : 
								((matrices_function) ? matrices_function(matrix_index) : MatrixXcd())
						);
		
		std::function<Block<MatrixXcd>(long,long)> matrix_block = [&](long I, long J) { return M.block(offsets[I], offsets[J], block_sizes[I], block_sizes[J]); };


		MatrixXcd epsilon = matrix_block(0,0);
		MatrixXcd g = epsilon.inverse();
		MatrixXcd epsilonsurf = epsilon;
		MatrixXcd alpha = matrix_block(1,0);
		MatrixXcd beta = matrix_block(0,1);

		auto valid = [&]() {

			if ((alpha > 1.0e-10).any() + (beta > 1.0e-10).any())
				return true;

			return false;
		};
		
		for (int iter = 0; valid(); iter++)
		{
			epsilon += (beta * g * alpha + alpha * g * beta);
			epsilonsurf += (alpha * g * beta);

			alpha = (alpha * g * alpha);
			beta = (beta * g * beta);
			
			g = epsilon.inverse();
		}

		epsilonsurf += (alpha * g * beta);
		
		g = epsilonsurf.inverse();

		computed_matrices[matrix_index] = matrix_block(2,1) * g * matrix_block(1,2);
	}
	
	void compute_matrix_from_right(const long &matrix_index) 
	{
		if((block_sizes == 0).any())
		{
			log() << "The sizes given do not fit in the matrices. Cannot guess the value."
			return;
		}
		
		ArrayXi offsets(4);
		offsets[0] = 0;
		// A cumulative sum (first,last,destination_start).
		std::partial_sum(&block_sizes[0], &block_sizes[2], &offsets[1]);
		
		MatrixXcd M = (matrices_array) ? 
						(
							(matrices_count == 1) ? 
								*matrices_array : 
								matrices_array[matrix_index]
						) :
						(
							(matrices_vector) ? 
								(*matrices_vector)[matrix_index] : 
								((matrices_function) ? matrices_function(matrix_index) : MatrixXcd())
						);
		
		std::function<Block<MatrixXcd>(long,long)> matrix_block = [&](long I, long J) { return M.block(offsets[I], offsets[J], block_sizes[I], block_sizes[J]); };


		MatrixXcd epsilon = matrix_block(2,2);
		MatrixXcd g = epsilon.inverse();
		MatrixXcd epsilonsurf = epsilon;
		MatrixXcd alpha = matrix_block(1,2);
		MatrixXcd beta = matrix_block(2,1);

		auto valid = [&]() {

			if ((alpha > 1.0e-10).any() + (beta > 1.0e-10).any())
				return true;

			return false;
		};
		
		for (int iter = 0; valid(); iter++)
		{
			epsilon += (beta * g * alpha + alpha * g * beta);
			epsilonsurf += (alpha * g * beta);

			alpha = (alpha * g * alpha);
			beta = (beta * g * beta);
			
			g = epsilon.inverse();
		}

		epsilonsurf += (alpha * g * beta);
		
		g = epsilonsurf.inverse();

		computed_matrices[matrix_index] = matrix_block(0,1) * g * matrix_block(1,0);
	}
	
public:
	
	enum compute_action {
		FromLeft,
		FromRight
		
	};

	void setBlockSizes(const Array3i &sizes) {
		if((sizes == 0).any())
		{
			log() << "The sizes given do not fit in the matrices. Cannot guess the value."
			return;
		}
	}
	
	void compute(const compute_action &action)
	{
		if(action == FromLeft)
			compute_matrix_from_left();
		else if(action == FromRight)
			compute_matrix_from_right();
	}

	void compute(const compute_action &action, const Array3i &sizes)
	{
		if((sizes == 0).any())
		{
			log() << "The sizes given do not fit in the matrices. Cannot guess the value."
			return;
		}
		
		if(action == FromLeft)
			compute_matrix_from_left();
		else if(action == FromRight)
			compute_matrix_from_right();
	}
	
	void computeIndex(const compute_action &action, const long &index)
	{
		if(action == FromLeft)
			compute_matrix_from_left(index);
		else if(action == FromRight)
			compute_matrix_from_right(index);
	}
	
	std::vector<MatrixXcd> &maitrices() const
	{
		return computed_matrices;
	}
		
};

std::ostream SemiInfiniteChain::null_stream(0);
bool SemiInfiniteChain::logging_enabled(false);

}

}

#endif  
