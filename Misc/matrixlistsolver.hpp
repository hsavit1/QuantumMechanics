/*
Header file for QuantumMechanics::MatrixListSolver:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/
#ifndef _MATRIXLISTSOLVER_H_
#define _MATRIXLISTSOLVER_H_

#include <functional>
#include <vector>

namespace QuantumMechanics {

template <typename MATRIXSOLVER>
class MatrixListSolver : public FeedbackObject {
public:
	typedef MATRIXSOLVER MatrixSolver;

private:
	const size_t input_count;
	const InputType * const input_array;
	const std::function<InputType(int)> input_function;
	const std::vector<InputType> * const input_vector;

	std::vector<MatrixSolver *> solvers;

public:
	template <typename InputType>
	MatrixListSolver(size_t n, const InputType * const M) :
		input_count(n),
		input_array(M),
		input_vector(nullptr),
		input_function(nullptr),
		solvers(n, nullptr)
	{ }

	template <typename InputType>
	MatrixListSolver(const std::vector<InputType> &M) :
		input_count(M.size()),
		input_array(nullptr),
		input_vector(&M),
		input_function(nullptr),
		solvers(n, nullptr)
	{ }

	template <typename InputType>
	MatrixListSolver(size_t n, const std::function<InputType(int)> &M) :
		input_count(n),
		input_array(nullptr),
		input_vector(nullptr),
		input_function(M),
		solvers(n, nullptr)
	{ }

	virtual ~MatrixListSolver() { }
	
	inline InputType input(const long &index) const {
		return (input_array) ?
			(
			(input_count == 1) ?
			*input_array :
			input_array[index]
			) :
			(
			(input_vector) ?
			(*input_vector)[index] :
			((input_function) ? input_function(index) : InputType())
			);
	}

	inline MatrixSolver *solver(const long &index) {
		return solvers[index];
	}

	inline const MatrixSolver *solver(const long &index) const {
		return solvers[index];
	}

	inline void compute(const long index, const ComputeType &action)
	{
		InputType M = input(i);
		if (block_row_sizes.size() > 0 && block_column_sizes.size() > 0)
			M.setInternalBlocks(block_row_sizes, block_column_sizes);
		
		MatrixSolver solver;
		
		solver.input(input(i));
		solver.compute(action);

		solutions[i] = solver.soluition();
	}

	inline void compute(const ComputeEnumType &action)
	{
		resetFeedback();
		const double delta = 1. / double(input_count);
		_Cilk_for(long m = 0; m < matrices_count; m++) {
			compute(m, action);
			updateFeedback(delta);
		}
	}
};

};

#endif //namespace _MATRIXLISTSOLVER_H_
