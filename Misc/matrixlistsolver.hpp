/*
Header file for QuantumMechanics::MatrixListSolver:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, S�ren Schou Gregersen <sorge@nanotech.dtu.dk>
*/
#ifndef _MATRIXLISTSOLVER_H_
#define _MATRIXLISTSOLVER_H_

#include <functional>
#include <vector>

namespace QuantumMechanics {

template <typename  MATRIXSOLVER>
class MatrixListSolver : public FeedbackObject {
public:
	typedef MATRIXSOLVER MatrixSolver;
	typedef MATRIXSOLVER::InputMatrixType InputMatrixType;
	typedef MATRIXSOLVER::OutputMatrixType OutputMatrixType;
	typedef MATRIXSOLVER::ComputeEnumType ComputeEnumType;

private:
	const size_t matrices_count;
	const size_t matrices_size;
	const  * const matrices_array;
	const std::function<InputMatrixType(int)> matrices_function;
	const std::vector<InputMatrixType> * const matrices_vector;

	ArrayXi block_sizes;
	ArrayXi block_offsets;

	std::vector<OutputMatrixType> solved_matrices;

public:
	MatrixListSolver() :
		matrices_count(0),
		matrices_size(0),
		matrices_array(nullptr),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		block_sizes(ArrayXi(1) << 0),
		block_offsets(ArrayXi(1) << 0)
	{ }

	MatrixListSolver(size_t n, const InputMatrixType * const M) :
		matrices_count(n),
		matrices_size((M->rows() == M->cols() && M->rows() > 0) ? M->rows() : 0),
		matrices_array(M),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		block_sizes(ArrayXi(1) << ((M->rows() == M->cols() && M->rows() > 0) ? M->rows() : 0)),
		block_offsets(ArrayXi(1) << 0),
		solved_matrices(std::vector<OutputMatrixType>(n, OutputMatrixType()))
	{ }

	MatrixListSolver(size_t n, const InputMatrixType * const M, const size_t &size) :
		matrices_count(n),
		matrices_size((M->rows() >= size && M->cols() >= size && size > 0) ? size : 0),
		matrices_array(M),
		matrices_vector(nullptr),
		matrices_function(nullptr),
		block_sizes(ArrayXi(1) << ((M->rows() >= size && M->cols() >= size && size > 0) ? size : 0)),
		block_offsets(ArrayXi(1) << 0),
		solved_matrices(std::vector<OutputMatrixType>(n, OutputMatrixType()))
	{ }

	MatrixListSolver(const std::vector<InputMatrixType> &M) :
		matrices_count(M.size()),
		matrices_size((M[0].rows() == M[0].cols()) ? M[0].rows() : 0),
		matrices_array(nullptr),
		matrices_vector(&M),
		matrices_function(nullptr),
		block_sizes(ArrayXi(1) << ((M[0].rows() == M[0].cols()) ? M[0].rows() : 0)),
		block_offsets(ArrayXi(1) << 0),
		solved_matrices(std::vector<OutputMatrixType>(n, OutputMatrixType()))
	{ }

	MatrixListSolver(const std::vector<InputMatrixType> &M, const size_t &size) :
		matrices_count(M.size()),
		matrices_size((M[0].rows() >= size && M[0].cols() >= size && size > 0) ? size : 0),
		matrices_array(nullptr),
		matrices_vector(&M),
		matrices_function(nullptr),
		block_sizes(ArrayXi(1) << ((M[0].rows() >= size && M[0].cols() >= size && size > 0) ? size : 0)),
		block_offsets(ArrayXi(1) << 0),
		solved_matrices(std::vector<OutputMatrixType>(n, OutputMatrixType()))
	{ }

	MatrixListSolver(size_t n, const std::vector<InputMatrixType> &M) :
		matrices_count(n),
		matrices_size((M[0].rows() == M[0].cols()) ? M[0].rows() : 0),
		matrices_array(nullptr),
		matrices_vector(&M),
		matrices_function(nullptr),
		block_sizes(ArrayXi(1) << ((M[0].rows() == M[0].cols()) ? M[0].rows() : 0)),
		block_offsets(ArrayXi(1) << 0),
		solved_matrices(std::vector<OutputMatrixType>(n, OutputMatrixType()))
	{ }

	MatrixListSolver(size_t n, const std::vector<InputMatrixType> &M, const size_t &size) :
		matrices_count(n),
		matrices_size((M[0].rows() >= size && M[0].cols() >= size && size > 0) ? size : 0),
		matrices_array(nullptr),
		matrices_vector(&M),
		matrices_function(nullptr),
		block_sizes(ArrayXi(1) << ((M[0].rows() >= size && M[0].cols() >= size && size > 0) ? size : 0)),
		block_offsets(ArrayXi(1) << 0),
		solved_matrices(std::vector<OutputMatrixType>(n, OutputMatrixType()))
	{ }

	MatrixListSolver(size_t n, const std::function<InputMatrixType(int)> &M, const size_t &size) :
		matrices_count(n),
		matrices_size((size > 0) ? size : 0),
		matrices_array(nullptr),
		matrices_vector(nullptr),
		matrices_function(M),
		block_sizes(ArrayXi(1) << ((size > 0) ? size : 0)),
		block_offsets(ArrayXi(1) << 0),
		solved_matrices(std::vector<OutputMatrixType>(n, OutputMatrixType()))
	{ }

	virtual ~MatrixListSolver() { }

	void setBlockSizes(const ArrayXi &sizes) {
		if (sizes.sum() > matrices_size)
		{
			block_sizes = ArrayXi(ArrayXi(1) << matrices_size);
			block_offsets = ArrayXi(ArrayXi(1) << 0);
			return;
		}

		block_sizes = sizes;
		const long block_count = block_sizes.size();
		block_offsets = ArrayXi(block_count);
		offsets[0] = 0;
		// A cumulative sum (first,last,destination_start).
		std::partial_sum(&block_sizes[0], &block_sizes[block_count - 1], &block_offsets[1]);
	}

	inline bool blockedMatrices() const {
		return (block_sizes.size() > 1);
	}

	inline InputMatrixType matrixAt(const long &index) const {
		return (matrices_array) ?
			(
			(matrices_count == 1) ?
			*matrices_array :
			matrices_array[index]
			) :
			(
			(matrices_vector) ?
			(*matrices_vector)[index] :
			((matrices_function) ? matrices_function(index) : Min())
			);
	}

	inline OutputMatrixType &solutionAt(const long &index) {
		return solved_matrices[index];
	}

	inline const OutputMatrixType &solutionAt(const long &index) const {
		return solved_matrices[index];
	}

	inline void computeAt(const long index, const ComputeEnumType &action)
	{
		MatrixSolver solver(matrixAt(i), matrices_size);
		if (blockedMatrices())
			solver.setBlockSizes(block_sizes);

		solver.compute(action);

		solved_matrices[i] = solver.soluition();
	}

	inline void compute(const ComputeEnumType &action)
	{
		resetFeedback();
		const double delta = 1. / double(matrices_count);
		_Cilk_for(long m = 0; m < matrices_count; m++) {
			computeAt(m, action);
			updateFeedback(delta);
		}
	}

protected:

	inline const Block<Eigen::Ref<const Min>> &blockAt(const Min &full, const long &i, const long &j) const {
		return full.block(block_offsets[i], block_offsets[j], block_sizes[i], block_sizes[j]);
	}

	inline Block<Eigen::Ref<Min>> &blockAt(const Min &full, const long &i, const long &j) {
		return full.block(block_offsets[i], block_offsets[j], block_sizes[i], block_sizes[j]);
	}

	inline Min zeroBlock(const long &i, const long &j) const {
		return Min::Zero(block_sizes(0), block_sizes(0));
	}

	inline long blockCount() const {
		return block_sizes.size();
	}

	inline const Mout &solutionAt(const long &index) const {
		return solved_matrices[index];
	}
};

};

#endif //namespace _MATRIXLISTSOLVER_H_