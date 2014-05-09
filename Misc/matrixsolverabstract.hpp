/*
Header file for QuantumMechanics::MatrixSolverAbstract:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/
#ifndef _MATRIXSOLVERABSTRACT_H_
#define _MATRIXSOLVERABSTRACT_H_

#include "..Math/Dense"

namespace QuantumMechanics {

template <typename  IN, typename  OUT, typename COMPUTEENUM>
class MatrixSolverAbstract {
public:
	typedef IN InputMatrixType;
	typedef OUT OutputMatrixType;
	typedef COMPUTEENUM ComputeEnumType;

private:
	const size_t matrix_size;
	const InputMatrixType &matrix;
	OutputMatrixType solution_matrix;

	ArrayXi block_sizes;
	ArrayXi block_offsets;

	static const InputMatrixType empty_matrix;

public:
	MatrixListSolver() :
		matrix_size(0),
		matrix(empty_matrix),
		block_sizes(ArrayXi(1) << 0),
		block_offsets(ArrayXi(1) << 0)
	{ }

	MatrixListSolver(const InputMatrixType &M) :
		matrix_size((M.rows() == M.cols() && M.rows() > 0) ? M.rows() : 0),
		matrix(M),
		block_sizes(ArrayXi(1) << ((M.rows() == M.cols() && M.rows() > 0) ? M.rows() : 0)),
		block_offsets(ArrayXi(1) << 0)
	{ }

	MatrixListSolver(const InputMatrixType &M, const size_t &size) :
		matrix_size((M.rows() >= size && M.cols() >= size && size > 0) ? size : 0),
		matrix(M),
		block_sizes(ArrayXi(1) << ((M.rows() >= size && M.cols() >= size && size > 0) ? size : 0)),
		block_offsets(ArrayXi(1) << 0)
	{ }

	MatrixListSolver(const InputMatrixType * const M) : 
		MatrixListSolver(*M)
	{ }

	MatrixListSolver(const InputMatrixType * const M, const size_t &size) :
		MatrixListSolver(*M, size)
	{ }

	virtual ~MatrixListSolver() { }

	virtual void setBlockSizes(const ArrayXi &sizes) {
		if (sizes.sum() > matrix_size)
		{
			block_sizes = ArrayXi(ArrayXi(1) << matrix_size);
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

	virtual inline bool blockedMatrix() const {
		return (block_sizes.size() > 1);
	}

	virtual inline void compute(const ComputeEnumType &action);

	inline const OutputMatrixType &solution() const {
		return solution_matrix;
	}

	inline OutputMatrixType &solution() {
		return solution_matrix;
	}

protected:
	virtual inline const Block<const InputMatrixType> block(const long &i, const long &j) const {
		return matrix.block(block_offsets[i], block_offsets[j], block_sizes[i], block_sizes[j]);
	}

	virtual inline const Block<const InputMatrixType> blockRange(const long &i, const long &j, const long &n, const long &m) const {
		return matrix.block(block_offsets[i], block_offsets[j], block_sizes.segment(i, i + n).sum(), block_sizes.segment(j, j + n).sum());
	}

	virtual inline const Block<const OutputMatrixType> solutionBlock(const long &i, const long &j) const {
		return solution_matrix.block(block_offsets[i], block_offsets[j], block_sizes[i], block_sizes[j]);
	}

	virtual inline const Block<const OutputMatrixType> solutionBlockRange(const long &i, const long &j, const long &n, const long &m) const {
		return solution_matrix.block(block_offsets[i], block_offsets[j], block_sizes.segment(i, i + n).sum(), block_sizes.segment(j, j + n).sum());
	}

	virtual inline const Block<Eigen::Reverse<const InputMatrixType, BothDirections>> reverseBlock(const long &i, const long &j) const {
		return matrix.reverse().block(block_offsets[i], block_offsets[j], block_sizes[i], block_sizes[j]);
	}

	virtual inline const Block<Eigen::Reverse<const InputMatrixType, BothDirections>> reverseBlockRange(const long &i, const long &j, const long &n, const long &m) const {
		return matrix.reverse().block(block_offsets[i], block_offsets[j], block_sizes.segment(i, i + n).sum(), block_sizes.segment(j, j + n).sum());
	}

	virtual inline const Block<Eigen::Reverse<const OutputMatrixType, BothDirections>> reverseSolutionBlock(const long &i, const long &j) const {
		return solution_matrix.reverse().block(block_offsets[i], block_offsets[j], block_sizes[i], block_sizes[j]);
	}

	virtual inline const Block<Eigen::Reverse<const OutputMatrixType, BothDirections>> reverseSolutionBlockRange(const long &i, const long &j, const long &n, const long &m) const {
		return solution_matrix.reverse().block(block_offsets[i], block_offsets[j], block_sizes.segment(i, i + n).sum(), block_sizes.segment(j, j + n).sum());
	}

	template <typename LOCAL>
	virtual inline LOCAL zeroBlock(const long &i, const long &j) const {
		return LOCAL::Zero(block_sizes(0), block_sizes(0));
	}

	template <typename LOCAL>
	virtual inline LOCAL zeroBlockRange(const long &i, const long &j, const long &n, const long &m) const {
		return LOCAL::Zero(block_sizes.segment(i, i + n).sum(), block_sizes.segment(j, j + n).sum());
	}

	template <typename LOCAL>
	virtual inline LOCAL reverseZeroBlock(const long &i, const long &j) const {
		return LOCAL::Zero(block_sizes.reverse()(0), block_sizes.reverse()(0));
	}

	template <typename LOCAL>
	virtual inline LOCAL reverseZeroBlockRange(const long &i, const long &j, const long &n, const long &m) const {
		return LOCAL::Zero(block_sizes.reverse().segment(i, i + n).sum(), block_sizes.reverse().segment(j, j + n).sum());
	}

	virtual inline long blockCount() const {
		return block_sizes.size();
	}

};

const InputMatrixType MatrixSolverAbstract::empty_matrix;

};

#endif //namespace _MATRIXSOLVERABSTRACT_H_