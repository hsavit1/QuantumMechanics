/*
Header file for QuantumMechanics::GreensFormalism::GreensMatrixSolver: 

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the 

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _GREENSFORMALISM_GREENSMATRIXSOLVER_H_
#define _GREENSFORMALISM_GREENSMATRIXSOLVER_H_

#include <iostream> // For printing useful messages.
#include <fstream> // For empty streams.

// Including two different storage-forms of matrices.
#include <functional>
#include <vector>

#include <Math/Dense>

#include <tbb/mutex.h>	// Includes a thread lock mutex.

namespace QuantumMechanics {

namespace GreensFormalism {
	
class GreensMatrixSolver {
	
private:
	const size_t matrices_count;
	const size_t matrices_size;
	const MatrixXcd * const inverse_matrices_array;
	const std::function<MatrixXcd(int)> inverse_matrices_function;
	const std::vector<MatrixXcd> * const inverse_matrices_vector;
	
	ArrayXi block_sizes;

	std::vector<MatrixXcd> computed_matrices;
	
	std::function<void(double)> progress_function;
	
public:
	GreensMatrixSolver() :
		matrices_count(0),
		matrices_size(0),
		inverse_matrices_array(nullptr),
		inverse_matrices_vector(nullptr),
		inverse_matrices_function(nullptr),
		block_sizes(ArrayXi(1) << 0;),
		progress_function(nullptr)
		{ }

	GreensMatrixSolver(const MatrixXcd &M) :
		matrices_count(1),
		matrices_size((M.rows() == M.cols() && M.rows() > 0) ? M.rows() : 0),
		inverse_matrices_array(&M),
		inverse_matrices_vector(nullptr),
		inverse_matrices_function(nullptr),
		block_sizes(ArrayXi(1) << ((M.rows() == M.cols() && M.rows() > 0) ? M.rows() : 0)),
		computed_matrices(std::vector<MatrixXcd>(1, MatrixXcd())),
		progress_function(nullptr)
		{ }

	GreensMatrixSolver(size_t n, const MatrixXcd * const M) :
		matrices_count(n),
		matrices_size((M->rows() == M->cols() && M->rows() > 0) ? M->rows() : 0),
		inverse_matrices_array(M),
		inverse_matrices_vector(nullptr),
		inverse_matrices_function(nullptr),
		block_sizes(ArrayXi(1) << ((M->rows() == M->cols() && M->rows() > 0) ? M->rows() : 0)),
		computed_matrices(std::vector<MatrixXcd>(n, MatrixXcd())),
		progress_function(nullptr)
		{ }

	GreensMatrixSolver(size_t n, const MatrixXcd * const M, const size_t &size) :
		matrices_count(n),
		matrices_size((M->rows() >= size && M->cols() >= size && size > 0) ? size : 0),
		inverse_matrices_array(M),
		inverse_matrices_vector(nullptr),
		inverse_matrices_function(nullptr),
		block_sizes(ArrayXi(1) << (M->rows() >= size && M->cols() >= size && size > 0) ? size : 0),
		computed_matrices(std::vector<MatrixXcd>(n, MatrixXcd())),
		progress_function(nullptr)
		{ }

	GreensMatrixSolver(size_t n, const std::vector<MatrixXcd> &M) :
		matrices_count(n),
		matrices_size((M[0].rows() == M[0].cols()) ? M[0].rows() : 0),
		inverse_matrices_array(nullptr),
		inverse_matrices_vector(&M),
		inverse_matrices_function(nullptr),
		block_sizes(ArrayXi(1) << ((M[0].rows() == M[0].cols()) ? M[0].rows() : 0)),
		computed_matrices(std::vector<MatrixXcd>(n, MatrixXcd())),
		progress_function(nullptr)
		{ }

	GreensMatrixSolver(size_t n, const std::vector<MatrixXcd> &M, const size_t &size) :
		matrices_count(n),
		matrices_size((M[0].rows() == M[0].cols() && M[0].rows() == size) ? size : 0),
		inverse_matrices_array(nullptr),
		inverse_matrices_vector(&M),
		inverse_matrices_function(nullptr),
		block_sizes(ArrayXi(1) << ((M[0].rows() == M[0].cols() && M[0].rows() == size) ? size : 0)),
		computed_matrices(std::vector<MatrixXcd>(n, MatrixXcd())),
		progress_function(nullptr)
		{ }

	GreensMatrixSolver(size_t n, const std::function<MatrixXcd(int)> &M, const size_t &size) :
		matrices_count(n),
		matrices_size((size > 0) ? size : 0),
		inverse_matrices_array(nullptr),
		inverse_matrices_vector(nullptr),
		inverse_matrices_function(M),
		block_sizes(ArrayXi(1) << ((size > 0) ? size : 0)),
		computed_matrices(std::vector<MatrixXcd>(n, MatrixXcd())),
		progress_function(nullptr)
		{ }

	virtual ~GreensMatrixSolver() { }
	
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
			return (std::clog << "GreensFormalism::GreensMatrixSolver message: ");
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
	
	void compute_full_matrix() 
	{
		computed_matrices = std::vector<MatrixXcd>(matrices_count, MatrixXcd(matrices_size, matrices_size));
		
		_Cilk_for (long m = 0; m < matrices_count; m++) {
			compute_full_matrix(m);
		}
	}
	
	void compute_first_block() 
	{
		computed_matrices = std::vector<MatrixXcd>(matrices_count, MatrixXcd(block_sizes[block_count - 1], block_sizes[block_count - 1]));
		
		_Cilk_for (long m = 0; m < matrices_count; m++) {
			compute_first_block(m);
		}
	}
	
	void compute_last_block() 
	{
		computed_matrices = std::vector<MatrixXcd>(matrices_count, MatrixXcd(block_sizes[block_count - 1], block_sizes[block_count - 1]));
		
		_Cilk_for (long m = 0; m < matrices_count; m++) {
			compute_last_block(m);
		}
	}
	
	void compute_first_block_column() 
	{
		computed_matrices = std::vector<MatrixXcd>(matrices_count, MatrixXcd(matrices_size, block_sizes[block_count - 1]));
		
		_Cilk_for (long m = 0; m < matrices_count; m++) {
			compute_first_block_column(m);
		}
	}
	
	void compute_last_block_column() 
	{
		computed_matrices = std::vector<MatrixXcd>(matrices_count, MatrixXcd(matrices_size, block_sizes[block_count - 1]));
		
		_Cilk_for (long m = 0; m < matrices_count; m++) {
			compute_last_block_column(m);
		}
	}
	
	void compute_full_matrix(const long &matrix_index) 
	{		
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
					
		computed_matrices[matrix_index] = M.inverse();
		
	}
	
	void compute_first_block(const long &matrix_index) 
	{		
		const long block_count = block_sizes.size();
		
		ArrayXi offsets(block_count+1);
		offsets[0] = 0;
		// A cumulative sum (first,last,destination_start).
		std::partial_sum(&block_sizes[0], &block_sizes[block_count - 1], &offsets[1]);
		
		
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
						).reverse();
						
		
		std::function< Block<MatrixXcd>(long,long) > block = [&](const long &i, const long &j)
		{
			return M.block(offsets[i],offsets[j],block_sizes[i],block_sizes[j]);
		};
						
		MatrixXcd g;
		MatrixXcd sigma = MatrixXcd::Zero(block_sizes(0),block_sizes(0));
		
		for (long b = 0; b < block_count; b++)
		{
			g = (block(b, b) - sigma).inverse();
			
			if(b != 0)
				sigma = block(b + 1, b) * g * block(b, b + 1);
		}
		
		computed_matrices[matrix_index] = g.reverse();
	}
	
	void compute_last_block(const long &matrix_index) 
	{
		const long block_count = block_sizes.size();
		
		ArrayXi offsets(block_count+1);
		offsets[0] = 0;
		// A cumulative sum (first,last,destination_start).
		std::partial_sum(&block_sizes[0], &block_sizes[block_count - 1], &offsets[1]);
				
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
						
		
		std::function< Block<MatrixXcd>(long,long) > block = [&](const long &i, const long &j)
		{
			return M.block(offsets[i],offsets[j],block_sizes[i],block_sizes[j]);
		};
						
		MatrixXcd g;
		MatrixXcd sigma = MatrixXcd::Zero(block_sizes(0),block_sizes(0));
		
		for (long b = 0; b < block_count; b++)
		{
			g = (block(b, b) - sigma).inverse();
			
			if(b != 0)
				sigma = block(b + 1, b) * g * block(b, b + 1);
		}
		
		computed_matrices[matrix_index] = g;
	}
	
	void compute_first_block_column(const long &matrix_index) 
	{		
		const long block_count = block_sizes.size();
		
		ArrayXi offsets(block_count+1);
		offsets[0] = 0;
		// A cumulative sum (first,last,destination_start).
		std::partial_sum(&block_sizes[0], &block_sizes[block_count - 1], &offsets[1]);
		
		
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
						).reverse();
						
		
		std::function< Block<MatrixXcd>(long,long) > block = [&](const long &i, const long &j)
		{
			return M.block(offsets[i],offsets[j],block_sizes[i],block_sizes[j]);
		};
		
		std::vector<MatrixXcd> g(block_count, MatrixXcd());
		MatrixXcd sigma = MatrixXcd::Zero(block_sizes(0),block_sizes(0));
		
		for (long b = 0; b < block_count; b++)
		{
			g[b] = (block(b, b) - sigma).inverse();
			
			if(b != 0)
				sigma = block(b + 1, b) * g[b] * block(b, b + 1);
		}
		
		computed_matrices[matrix_index].resize(matrix_size,block_sizes[block_count - 1]);
		std::function< Block< Reverse<MatrixXcd,BothDirections> >(long,long) > col_block = [&](const long &i)
		{
			return computed_matrices[matrix_index].reverse().block(offsets[i],0,block_sizes[i],block_sizes[block_count-1]);
		};
		
		col_block(block_count - 1) = g[block_count - 1];
		
		for (long b = 2; b <= block_count; b++)
		{
			col_block(block_count - b) = g[block_count - b] * block(b, b + 1) * g[block_count - b + 1];
		}
	}
	
	void compute_last_block_column(const long &matrix_index) 
	{		
		const long block_count = block_sizes.size();
		
		ArrayXi offsets(block_count+1);
		offsets[0] = 0;
		// A cumulative sum (first,last,destination_start).
		std::partial_sum(&block_sizes[0], &block_sizes[block_count - 1], &offsets[1]);
		
		
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
						
		
		std::function< Block<MatrixXcd>(long,long) > block = [&](const long &i, const long &j)
		{
			return M.block(offsets[i],offsets[j],block_sizes[i],block_sizes[j]);
		};
		
		std::vector<MatrixXcd> g(block_count, MatrixXcd());
		MatrixXcd sigma = MatrixXcd::Zero(block_sizes(0),block_sizes(0));
		
		for (long b = 0; b < block_count; b++)
		{
			g[b] = (block(b, b) - sigma).inverse();
			
			if(b != 0)
				sigma = block(b + 1, b) * g[b] * block(b, b + 1);
		}
		
		computed_matrices[matrix_index].resize(matrix_size,block_sizes[block_count - 1]);
		std::function< Block<MatrixXcd>(long,long) > col_block = [&](const long &i)
		{
			return computed_matrices[matrix_index].block(offsets[i],0,block_sizes[i],block_sizes[block_count-1]);
		};
		
		col_block(block_count - 1) = g[block_count - 1];
		
		for (long b = 2; b <= block_count; b++)
		{
			col_block(block_count - b) = g[block_count - b] * block(b, b + 1) * g[block_count - b + 1];
		}
	}
	
public:
	
	enum compute_action {
		FullMatrix,
		FirstBlock,
		LastBlock,
		FirstBlockColumn,
		LastBlockColumn
	};

	void setBlockSizes(const ArrayXi &sizes) {
		if(sizes.sum() > matrices_size)
		{
			log() << "The sizes given do not fit in the matrices. Instead the full matrix is used."
			return;
		}
		
		block_sizes = sizes;
	}
	
	void compute(compute_action action = FullMatrix)
	{
		switch(action)
		{
		case FullMatrix:
			compute_full_matrix();
			break;
		case FirstBlock:
			compute_first_block();
			break;
		case LastBlock:
			compute_last_block();
			break;
		case FirstBlockColumn:
			compute_first_block_column();
			break;
		case LastBlockColumn:
			compute_last_block_column();
			break;
		}
	}

	void compute(compute_action action, const VectorXi &sizes)
	{
		setBlockSizes(sizes);
		compute(action);
	}
	
	void computeIndex(const long &index, compute_action action = FullMatrix)
	{
		switch(action)
		{
		case FullMatrix:
			compute_full_matrix(index);
			break;
		case FirstBlock:
			compute_first_block(index);
			break;
		case LastBlock:
			compute_last_block(index);
			break;
		case FirstBlockColumn:
			compute_first_block_column(index);
			break;
		case LastBlockColumn:
			compute_last_block_column(index);
			break;
		}
	}
	
	std::vector<MatrixXcd> &maitrices() const
	{
		return computed_matrices;
	}
};

std::ostream GreensMatrixSolver::null_stream(0);
bool GreensMatrixSolver::logging_enabled(false);

}

}

#endif  