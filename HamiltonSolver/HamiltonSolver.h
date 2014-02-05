/*
 * HamiltonSolver.h
 *
 *  Created on: 05/02/2014
 *      Author: gregersen
 */

#ifndef HAMILTONSOLVER_H_
#define HAMILTONSOLVER_H_

#include <Eigen/Dense>
#include <functional>

namespace QuantumMechanics {

using namespace Eigen;

class HamiltonSolver {
public:
	HamiltonSolver();
	virtual ~HamiltonSolver();
};

} /* namespace QuantumMechanics */

#endif /* HAMILTONSOLVER_H_ */
