/*
 * HamiltonSolver.cpp
 *
 *  Created on: 05/02/2014
 *      Author: gregersen
 */

#include "HamiltonSolver.h"

using namespace Eigen;

namespace QuantumMechanics {

HamiltonSolver::HamiltonSolver(std::function<MatrixXcd(int)> H) {
	// TODO Auto-generated constructor stub
	if(H == nullptr)
		return;
}

HamiltonSolver::~HamiltonSolver() {
	// TODO Auto-generated destructor stub
}

} /* namespace QuantumMechanics */
