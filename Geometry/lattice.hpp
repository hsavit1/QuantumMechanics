/*
Header file for QuantumMechanics::Geometry::Lattice:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/
#ifndef _EIGENSYSTEM_HERMITIANSOLVER_H_
#define _EIGENSYSTEM_HERMITIANSOLVER_H_

#include <iostream> // For printing useful messages.
#include <fstream> // For empty streams.

#include <Math/Dense>

namespace QuantumMechanics {

namespace Geometry {

class Lattice {

private:
};


std::ostream Lattice::null_stream(0);
bool Lattice::logging_enabled(false);

}

}

#endif                                                      // _EIGENSYSTEM_HERMITIANSOLVER_H_