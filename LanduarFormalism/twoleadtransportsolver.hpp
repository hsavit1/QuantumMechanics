/*
Header file for QuantumMechanics::LanduarFormalism::TwoLeadTransportSolver: 

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the 

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _GREENSFORMALISM_GREENSSOLVER_H_
#define _GREENSFORMALISM_GREENSSOLVER_H_

#include "../Misc/MatrixSolverAbstract"
#include "../Misc/LoggingOject"

#include "../GreensFormalism/GreensSolver"
#include "../GreensFormalism/ChainSolver"

namespace QuantumMechanics {

namespace LanduarFormalism{

	enum TwoLeadTransportCalculation {
		LeftToRight,
		RightToLeft,
		CurrentsLeftToRight,
		CurrentsRightToLeft
	};
	
class TwoLeadTransportSolver : public MatrixSolverAbstract<MatrixXcd, MatrixXcd, TwoLeadTransportCalculation> {

	MatrixXcd left_matrix, right_matrix;
	MatrixXcd left_coupling, right_coupling;
	ArrayXi left_block_sizes, right_block_sizes;
	ArrayXi left_block_offsets, right_block_sizes;

public:
	TwoLeadTransportSolver() : MatrixSolverAbstract() { }

	TwoLeadTransportSolver(const MatrixXcd &M) : MatrixSolverAbstract(M) { }

	TwoLeadTransportSolver(const MatrixXcd &M, const size_t &size) : MatrixSolverAbstract(M, size) { }

	TwoLeadTransportSolver(const MatrixXcd * const M) : MatrixSolverAbstract(M) { }

	TwoLeadTransportSolver(const MatrixXcd * const M, const size_t &size) : MatrixSolverAbstract(M, size) { }

	static LoggingObject log;

	void setLeads(const MatrixXcd &left, const MatrixXcd &right, const ArrayXi &left_sizes, const ArrayXi &right_sizes) {

		left_matrix = left.block(0, 0, left_sizes.segments(0, 1).sum(), left_sizes[0]);
		right_matrix = right.reverse().block(0, 0, right_sizes.segments(0, 1).sum(), right_sizes[0]);
		left_coupling = left.block(left_sizes.segments(0, 1).sum(), left_sizes.segments(0, 1).sum(), left_sizes[2], left_sizes[2]);
		right_coupling = right.reverse().block(right_sizes.segments(0, 1).sum(), right_sizes.segments(0, 1).sum(), right_sizes[2], right_sizes[2]);
		left_block_sizes = left_sizes;
		right_block_sizes = right_sizes;
	}

	

protected:
	void compute_left_to_right()
	{
		using GreensFormalism;

		ChainSolver left_chain(left_matrix);
		ChainSolver right_chain(right_matrix);

		left_chain.compute(LeftSemiInfinite);
		right_chain.compute(RightSemiInfinite);

		MatrixXcd sigma_left = left_coupling * left_chain.solution() * left_coupling.adjoint();
		MatrixXcd sigma_right = right_coupling.adjoint() * right_chain.solution() * right_coupling;

		GreensSolver solver(matrix - sigma_left - sigma_right);

		solver.compute(FirstBlock);

		sigma_left = block(sigma_left, 0, 0).eval();
		sigma_right = solver.reducedSigma();

		sigma_left = (sigma_left - sigma_left * std::complex(0, 1));
		sigma_right = (sigma_right - sigma_right * std::complex(0, 1));

		solution_matrix = Eigen::trace(sigma_left * solver.solution() * sigma_right * solver.solution().adjoint());
	}

	void compute_right_to_left()
	{
		using GreensFormalism;

		ChainSolver left_chain(left_matrix);
		ChainSolver right_chain(right_matrix);

		left_chain.compute(LeftSemiInfinite);
		right_chain.compute(RightSemiInfinite);

		MatrixXcd sigma_left = left_coupling * left_chain.solution() * left_coupling.adjoint();
		MatrixXcd sigma_right = right_coupling.adjoint() * right_chain.solution() * right_coupling;

		GreensSolver solver(matrix - sigma_left - sigma_right);

		solver.compute(LastBlock);

		sigma_left = solver.reducedSigma();
		sigma_right = block(sigma_right, 0, 0).eval();

		sigma_left = (sigma_left - sigma_left * std::complex(0, 1));
		sigma_right = (sigma_right - sigma_right * std::complex(0, 1));

		solution_matrix = Eigen::trace(sigma_right * solver.solution() * sigma_left * solver.solution().adjoint());
	}

	void compute_currents_left_to_right()
	{
		solution_matrix = matrix.inverse();
	}

	void compute_currents_right_to_left()
	{
		solution_matrix = matrix.inverse();
	}
	
public:
	void compute(const TwoLeadTransportCalculation &action)
	{
		switch(action)
		{
		case LeftToRight:
			compute_left_to_right();
			break;
		case RightToLeft:
			compute_right_to_left();
			break;
		case CurrentsLeftToRight:
			compute_currents_left_to_right();
			break;
		case CurrentsRightToLeft:
			compute_currents_right_to_left();
			break;
		}
	}

	void compute(TwoLeadTransportCalculation action, const ArrayXi &sizes)
	{
		setBlockSizes(sizes);
		compute(action);
	}
};

LoggingObject TwoLeadTransportSolver::log("LanduarFormalism::TwoLeadTransportSolver", false);

}

}

#include "../Misc/MatrixListSolver"

typedef MatrixListSolver<TwoLeadTransportSolver> TwoLeadTransportListSolver;

#endif  