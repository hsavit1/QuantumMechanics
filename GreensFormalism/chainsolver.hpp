/*
Header file for QuantumMechanics::GreensFormalism::ChainSolver: 

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the 

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _GREENSFORMALISM_CHAINSOLVER_H_
#define _GREENSFORMALISM_CHAINSOLVER_H_

#include "Math/Dense"
#include "Misc/LoggingObject"

namespace QuantumMechanics {

namespace GreensFormalism {

	enum ResultType {
		SurfaceGreensMatrix
	};
	
class ChainSolver{

	const BlockMatrixXcd H;
	const BlockMatrixXcd V;
	BlockMatrixXcd G;

	static LoggingObject log;

public:

	long max_iterations;

	ChainSolver(const BlockMatrixXcd &h, const BlockMatrixXcd &v) : H(h), V(v), max_iterations(1000) { }

	ChainSolver(const MatrixXcd &h, const MatrixXcd &v) : H(h), V(v), max_iterations(1000) { }

	static inline void enableLog()
	{
		log.enable();
	}

protected:
	void compute_matrix()
	{
		const long block_count = (H.isSquare() || H.blockRows() < H.blockCols() ? H.blockRows() : H.blockCols());

		log() << "Preparing to calculate the surface solution of " << block_count << "-by-" << block_count << " blocks chain parts." << std::endl;

		MatrixXcd epsilon = H;
		G = epsilon.inverse();
		MatrixXcd epsilonsurf = epsilon;

		MatrixXcd alpha = V;
		MatrixXcd beta = V.adjoint();

		auto valid = [&]() {

			if (alpha.isZero(1.0e-12) && beta.isZero(1.0e-12))
				return true;

			return false;
		};

		for (int iter = 0; iter < max_iterations && valid(); iter++)
		{
			epsilon += (beta * G * alpha + alpha * G * beta);
			epsilonsurf += (alpha * G * beta);

			alpha = (alpha * G * alpha);
			beta = (beta * G * beta);

			G = epsilon.inverse();
		}

		epsilonsurf += (alpha * G * beta);

		G = epsilonsurf.inverse();
	}
		
public:
	inline void compute(const ResultType &action)
	{
		if (action == SurfaceGreensMatrix)
			compute_matrix();
	}

	const BlockMatrixXcd &greensMatrix() const {
		return G;
	}
};

LoggingObject ChainSolver::log("GreensFormalism::ChainSolver", false);

}

}

#endif  
