/*
Header file for QuantumMechanics::LanduarFormalism::TwoLeadTransportSolver: 

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the 

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _LANDUARFORMALISM_TWOLEADTRANSPORTSOLVER_H_
#define _LANDUARFORMALISM_TWOLEADTRANSPORTSOLVER_H_

#include "Misc/LoggingObject"

#include "GreensFormalism/GreensSolver"
#include "GreensFormalism/ChainSolver"

namespace QuantumMechanics {

namespace LanduarFormalism{

	enum TwoLeadTransportCalculation {
		LeftToRight,
		RightToLeft,
		CurrentsLeftToRight,
		CurrentsRightToLeft
	};
	
class TwoLeadTransportSolver{

	/*
	We assume this type of matrix!
	
	[	h_ll	v_ll	0		0		0		]
	[											]
	[	v_ll*	h_ll	v_l		0		0		]
	[											]
	[	0		v_l*	h_d		v_r		0		]
	[											]
	[	0		0		v_r*	h_rl	v_rl	]
	[											]
	[	0		0		0		v_rl*	h_rl	]

	Note that all of these matrices are assumed a block matrix!
	
	*/

	BlockMatrixXcd full;

	// left lead chain!
	BlockMatrixXcd h_ll;
	BlockMatrixXcd v_ll;

	// device couplings and isolated hamilton!
	BlockMatrixXcd v_l;
	BlockMatrixXcd h_d;
	BlockMatrixXcd v_r;

	// right lead chain!
	BlockMatrixXcd h_rl;
	BlockMatrixXcd v_rl;

	double transport;

	MatrixXd current;

	static LoggingObject log;

public:
	// as default, the h_ii and v_ii are assumed as a single block.
	TwoLeadTransportSolver(const BlockMatrixXcd &m) : 
		full(m),

		h_ll(m.blocks(0, 0, 1, 1)),
		v_ll(m.blocks(0, 1, 1, 1)),

		v_l(m.blocks(1, 3, 1, m.blockRows() - 4)),
		h_d(m.blocks(3, 3, m.blockRows() - 4, m.blockRows() - 4)),
		v_r(m.blocks(3, -2, m.blockRows() - 4, 1)),

		h_rl(m.blocks(-1, -1, 1, 1)),
		v_rl(m.blocks(-2, -1, 1, 1))
		{}

	void setLeftLeadBlockCount(const long &left_lead_count)
	{
		// Note that the lead cell are square and have equal size!
		const long right_lead_count = h_rl.blockRows();
		const long device_block_count = full.blockRows() - 2 * left_lead_count - 2 * right_lead_count;

		h_ll.assignToReference(full.blocks(0, 0, left_lead_count, left_lead_count));
		v_ll.assignToReference(full.blocks(0, left_lead_count, left_lead_count, left_lead_count));


		v_l.assignToReference(full.blocks(left_lead_count, 2 * left_lead_count, left_lead_count, device_block_count));
		h_d.assignToReference(full.blocks(2 * left_lead_count, 2 * left_lead_count, device_block_count, device_block_count));
		v_r.assignToReference(full.blocks(2 * left_lead_count, -2 * left_lead_count, device_block_count, right_lead_count));
	}

	void setRightLeadBlockCount(const long &right_lead_count)
	{
		// Note that the lead cell are square and have equal size!
		const long left_lead_count = h_ll.blockRows();
		const long device_block_count = full.blockRows() - 2 * left_lead_count - 2 * right_lead_count;

		h_rl.assignToReference(full.blocks(-right_lead_count, -right_lead_count, right_lead_count, right_lead_count));
		v_rl.assignToReference(full.blocks(-2 * right_lead_count, -right_lead_count, right_lead_count, right_lead_count));


		v_l.assignToReference(full.blocks(left_lead_count, 2 * left_lead_count, left_lead_count, device_block_count));
		h_d.assignToReference(full.blocks(2 * left_lead_count, 2 * left_lead_count, device_block_count, device_block_count));
		v_r.assignToReference(full.blocks(2 * left_lead_count, -2 * left_lead_count, device_block_count, right_lead_count));
	}

protected:
	void compute_left_to_right()
	{
		using namespace GreensFormalism;

		ChainSolver left_chain(h_ll, v_ll);
		ChainSolver right_chain(h_rl, v_rl);

		left_chain.compute(SurfaceGreensMatrix);
		right_chain.compute(SurfaceGreensMatrix);
		
		BlockMatrixXcd sigma_left = h_d;
		BlockMatrixXcd sigma_right = h_d;

		sigma_left = v_l * left_chain.greensMatrix() * v_l.adjoint();
		sigma_right = v_r.adjoint() * right_chain.greensMatrix() * v_r;

		GreensSolver solver(BlockMatrixXcd(h_d - sigma_left - sigma_right).withBlocks(h_d));

		solver.compute(FirstBlock);

		sigma_left = sigma_left.block(0, 0);
		sigma_right = solver.reducedSigma();

		sigma_left = (sigma_left.adjoint() - sigma_left).eval() * std::complex<double>(0, 1);
		sigma_right = (sigma_right.adjoint() - sigma_right).eval() * std::complex<double>(0, 1);

		transport = (sigma_right * solver.greensMatrix() * sigma_left * solver.greensMatrix().adjoint()).trace().real();
	}

	void compute_right_to_left()
	{
		using namespace GreensFormalism;

		ChainSolver left_chain(h_ll, v_ll.adjoint());
		ChainSolver right_chain(h_rl, v_rl.adjoint());

		left_chain.compute(SurfaceGreensMatrix);
		right_chain.compute(SurfaceGreensMatrix);

		BlockMatrixXcd sigma_left = h_d;
		BlockMatrixXcd sigma_right = h_d;

		sigma_left = v_l.adjoint() * left_chain.greensMatrix() * v_l;
		sigma_right = v_r * right_chain.greensMatrix() * v_r.adjoint();

		GreensSolver solver(BlockMatrixXcd(h_d - sigma_left - sigma_right).withBlocks(h_d));

		solver.compute(LastBlock);

		sigma_left = solver.reducedSigma();
		sigma_right = sigma_right.block(-1, -1);

		sigma_left = (sigma_left.adjoint() - sigma_left).eval() * std::complex<double>(0, 1);
		sigma_right = (sigma_right.adjoint() - sigma_right).eval() * std::complex<double>(0, 1);

		transport = (sigma_left * solver.greensMatrix() * sigma_right * solver.greensMatrix().adjoint()).trace().real();
	}

	void compute_currents_left_to_right()
	{
		current = full.inverse().real();
	}

	void compute_currents_right_to_left()
	{
		current = full.inverse().real();
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
};

LoggingObject TwoLeadTransportSolver::log("LanduarFormalism::TwoLeadTransportSolver", false);

}

}

#endif  