/*
Header file for QuantumMechanics::GreensFormalism::ChainSolver: 

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the 

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _GREENSFORMALISM_CHAINSOLVER_H_
#define _GREENSFORMALISM_CHAINSOLVER_H_

#include "../Misc/MatrixSolverAbstract"
#include "../Misc/LoggingOject"
#include "../Misc/FeedbackObject"

namespace QuantumMechanics {

namespace GreensFormalism {

	enum ChainType {
		LeftSemiInfinite,
		RightSemiInfinite
	};
	
class ChainSolver : public MatrixSolverAbstract<MatrixXcd, MatrixXcd, ChainType>  {

public:

	long max_iterations;

	ChainSolver() : MatrixSolverAbstract(), max_iterations(100) {}

	ChainSolver(const MatrixXcd &M) : MatrixSolverAbstract(M), max_iterations(100) { }

	ChainSolver(const MatrixXcd &M, const size_t &size) : MatrixSolverAbstract(M, size), max_iterations(100) { }

	ChainSolver(const MatrixXcd * const M) : MatrixSolverAbstract(M), max_iterations(100) { }

	ChainSolver(const MatrixXcd * const M, const size_t &size) : MatrixSolverAbstract(M, size), max_iterations(100) { }

	static LoggingObject log;
	
protected:	
	void compute_semi_infinite_matrix_from_left() 
	{
		if (blockCount() != 2)
		{
			log() << "The chain matrix is not of the type [h,v], where both h and v are block matrices. Cannot move on alone.";
			return;
		}
				
		MatrixXcd epsilon = block(0,0);
		MatrixXcd g = epsilon.inverse();
		MatrixXcd epsilonsurf = epsilon;

		if (matrix.rows() >= block_sizes[0] && matrix.cols() >= block_sizes.segments(0, 1).sum())
		{
			MatrixXcd alpha = block(0, 1).adjoint();
			MatrixXcd beta = block(0, 1);
		}
		else if (matrix.cols() >= block_sizes[0] && matrix.rows() >= block_sizes.segments(0, 1).sum())
		{
			MatrixXcd alpha = block(1, 0);
			MatrixXcd beta = block(1, 0).adjoint();
		}
		else
		{
			log() << "The [h,v] blocks cannot be found. Cannot move on alone.";
			return;
		}

		auto valid = [&]() {

			if ((alpha > 1.0e-10).any() + (beta > 1.0e-10).any())
				return true;

			return false;
		};
		
		for (int iter = 0; iter < max_iterations && valid(); iter++)
		{
			epsilon += (beta * g * alpha + alpha * g * beta);
			epsilonsurf += (alpha * g * beta);

			alpha = (alpha * g * alpha);
			beta = (beta * g * beta);
			
			g = epsilon.inverse();
		}

		epsilonsurf += (alpha * g * beta);
		
		solution_matrix = epsilonsurf.inverse();
	}
	
	void compute_semi_infinite_matrix_from_right()
	{
		if (blockCount() != 2)
		{
			log() << "The chain matrix is not of the type [h,v], where both h and v are block matrices. Cannot move on alone."
				return;
		}

		MatrixXcd epsilon = block(0, 0);
		MatrixXcd g = epsilon.inverse();
		MatrixXcd epsilonsurf = epsilon;
		if (matrix.rows() >= block_sizes[0] && matrix.cols() >= block_sizes.segments(0, 1).sum())
		{
			MatrixXcd alpha = block(0, 1);
			MatrixXcd beta = block(0, 1).adjoint();
		}
		else if (matrix.cols() >= block_sizes[0] && matrix.rows() >= block_sizes.segments(0, 1).sum())
		{
			MatrixXcd alpha = block(1, 0).adjoint();
			MatrixXcd beta = block(1, 0);
		}
		else
		{
			log() << "The [h,v] blocks cannot be found. Cannot move on alone.";
			return;
		}

		auto valid = [&]() {

			if ((alpha > 1.0e-10).any() + (beta > 1.0e-10).any())
				return true;

			return false;
		};

		for (int iter = 0; iter < max_iterations && valid(); iter++)
		{
			epsilon += (beta * g * alpha + alpha * g * beta);
			epsilonsurf += (alpha * g * beta);

			alpha = (alpha * g * alpha);
			beta = (beta * g * beta);

			g = epsilon.inverse();
		}

		epsilonsurf += (alpha * g * beta);

		solution_matrix = epsilonsurf.inverse();
	}
	
public:
	inline void compute(const ChainType &action)
	{
		if (action == LeftSemiInfinite)
			compute_semi_infinite_matrix_from_left();
		else if (action == RightSemiInfinite)
			compute_semi_infinite_matrix_from_right();
	}

	void compute(const ChainType &action, const Array3i &sizes)
	{
		setBlockSizes(sizes);
		compute(action);
	}
};

LoggingObject ChainSolver::log("GreensFormalism::ChainSolver", false);

}

}

#include "../Misc/MatrixListSolver"

typedef MatrixListSolver<ChainSolver> ChainListSolver;

#endif  
