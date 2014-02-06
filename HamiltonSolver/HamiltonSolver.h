/*
 * HamiltonSolver.h
 *
 *  Created on: 05/02/2014
 *      Author: gregersen
 */

#ifndef HAMILTONSOLVER_H_
#define HAMILTONSOLVER_H_

#include <iostream>
#include <functional>
#include <Eigen/Dense>


namespace QuantumMechanics {

using namespace Eigen;

class HamiltonSolver {
public:
	HamiltonSolver(std::function<MatrixXcd(int)> H = nullptr);
	virtual ~HamiltonSolver();

private:
	std::function<MatrixXcd(int)> hamilton;
};

} /* namespace QuantumMechanics */

#endif /* HAMILTONSOLVER_H_ */
