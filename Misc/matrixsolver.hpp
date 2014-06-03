/*
Header file for QuantumMechanics::MatrixSolver:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/
#ifndef _MATRIXSOLVER_H_
#define _MATRIXSOLVER_H_

namespace QuantumMechanics {

template <typename INPUT, typename COMPUTEENUM>
class MatrixSolver {
	typedef INPUT InputType;
	typedef COMPUTEENUM ComputeType;

public:
	virtual inline void compute(const ComputeType &action) = 0;
};

};

#endif //namespace _MATRIXSOLVERABSTRACT_H_
