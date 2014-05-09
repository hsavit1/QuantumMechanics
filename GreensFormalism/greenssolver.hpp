/*
Header file for QuantumMechanics::GreensFormalism::GreensSolver: 

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the 

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _GREENSFORMALISM_GREENSSOLVER_H_
#define _GREENSFORMALISM_GREENSSOLVER_H_

#include "../Misc/MatrixSolverAbstract"
#include "../Misc/LoggingOject"
#include "../Misc/FeedbackObject"

namespace QuantumMechanics {

namespace GreensFormalism {

	enum GreenMatrixSubType {
		FullMatrix,
		FirstBlock,
		LastBlock,
		FirstBlockColumn,
		LastBlockColumn
	};
	
class GreensSolver : public MatrixSolverAbstract<MatrixXcd, MatrixXcd, GreenMatrixSubType> {

public:
	GreensSolver() : MatrixSolverAbstract() { }

	GreensSolver(const MatrixXcd &M) : MatrixSolverAbstract(M) { }

	GreensSolver(const MatrixXcd &M, const size_t &size) : MatrixSolverAbstract(M, size) { }

	GreensSolver(const MatrixXcd * const M) : MatrixSolverAbstract(M) { }

	GreensSolver(const MatrixXcd * const M, const size_t &size) : MatrixSolverAbstract(M, size) { }

	static LoggingObject log;

protected:
	void compute_full_matrix() 
	{
		solution_matrix = matrix.inverse();
	}
	
	void compute_first_block() 
	{		
		if (!blockedMatrices())
		{
			solution_matrix = matrix.inverse();
			return;
		}

		const long block_count = blockCount();
						
		MatrixXcd g;
		MatrixXcd sigma = reverseZeroBlock(0, 0);
		
		for (long b = 0; b < block_count; b++)
		{
			g = (reverseBlock(b, b) - sigma).inverse();
			
			if(b != 0)
				sigma = reverseBlock(b + 1, b) * g * reverseBlock(b, b + 1);
		}
		
		solution_matrix.reverse() = g;
	}
	
	void compute_last_block() 
	{
		if (!blockedMatrices())
		{
			solution_matrix = matrix.inverse();
			return;
		}

		const long block_count = blockCount();

		MatrixXcd g;
		MatrixXcd sigma = zeroBlock(0, 0);
		
		for (long b = 0; b < block_count; b++)
		{
			g = (block(b, b) - sigma).inverse();
			
			if(b != 0)
				sigma = block(b + 1, b) * g * block(b, b + 1);
		}

		solution_matrix = g;
	}
	
	void compute_first_block_column() 
	{
		const long block_count = blockCount();
		
		std::vector<MatrixXcd> g(block_count, MatrixXcd());
		MatrixXcd sigma = reverseZeroBlock(0, 0);
		
		for (long b = 0; b < block_count; b++)
		{
			g[b] = (reverseBlock(b, b) - sigma).inverse();
			
			if(b != 0)
				sigma = reverseBlock(M, b + 1, b) * g[b] * reverseBlock(M, b, b + 1);
		}

		solution_matrix = reverseZeroBlockRange(0, 0, 0, block_count - 1);

		reverseSolutionBlock(0, block_count - 1) = g[block_count - 1];

		for (long b = 2; b <= block_count; b++)
		{
			reverseSolutionBlock(0, block_count - b) = g[block_count - b] * reverseBlock(b, b + 1) * g[block_count - b + 1];
		}
	}
	
	void compute_last_block_column() 
	{
		const long block_count = blockCount();

		std::vector<MatrixXcd> g(block_count, MatrixXcd());
		MatrixXcd sigma = zeroBlock(0, 0);

		for (long b = 0; b < block_count; b++)
		{
			g[b] = (block(b, b) - sigma).inverse();

			if (b != 0)
				sigma = block(M, b + 1, b) * g[b] * block(M, b, b + 1);
		}

		solution_matrix = zeroBlockRange(0, 0, 0, block_count - 1);

		solutionBlock(0, block_count - 1) = g[block_count - 1];

		for (long b = 2; b <= block_count; b++)
		{
			solutionBlock(0, block_count - b) = g[block_count - b] * block(b, b + 1) * g[block_count - b + 1];
		}
	}
	
public:
	void compute(const GreenMatrixSubType &action)
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

	void compute(compute_action action, const ArrayXi &sizes)
	{
		setBlockSizes(sizes);
		compute(action);
	}
};

LoggingObject GreensMatrixSolver::log("GreensFormalism::GreensSolver", false);

}

}

#include "../Misc/MatrixListSolver"

typedef MatrixListSolver<GreensSolver> GreensListSolver;

#endif  