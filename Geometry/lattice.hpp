/*
Header file for QuantumMechanics::Geometry::Lattice:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/
#ifndef _GEOMETRY_LATTICE_H_
#define _GEOMETRY_LATTICE_H_

#include <iostream> // For printing useful messages.
#include <fstream> // For empty streams.

#include <Math/Dense>

namespace QuantumMechanics {

namespace Geometry {

class Lattice {

private:
	MatrixXd lattice_matrix;

public:
	inline set(const MatrixXd &lattice)
	{
		if (lattice.rows() > 0 && lattice.rows() < 4 && lattice.cols() > 0 && lattice.cols() < 4)
			lattice_matrix = lattice;
	}

	inline set(const std::vector<VectorXd> &lattice)
	{
		const long lattice_dimension = lattice.size();

		if (lattice_dimension == 1)
		{
			long vector_size = lattice[0].size();
			if (vector_size > 0 && vector_size < 4)
				lattice_matrix = MatrixXd(vector_size, lattice_dimension) << lattice[0];
		}
		else if (lattice_dimension == 2)
		{
			long vector_size = lattice[0].size();
			long vector_size = (vector_size > lattice[1].size()) ? vector_size : lattice[1].size();
			if (vector_size > 0 && vector_size < 4)
				lattice_matrix = MatrixXd(vector_size, lattice_dimension) << lattice[0] << lattice[1];
		}
		else if (lattice_dimension == 3)
		{
			long vector_size = lattice[0].size();
			long vector_size = (vector_size > lattice[1].size()) ? vector_size : lattice[1].size();
			long vector_size = (vector_size > lattice[2].size()) ? vector_size : lattice[2].size();
			if (vector_size > 0 && vector_size < 4)
				lattice_matrix = MatrixXd(vector_size, lattice_dimension) << lattice[0] << lattice[1] << lattice[2];
		}
		else
			return;
	}

	inline set(const VectorXd &v1)
	{
		long vector_size = v1.size();
		if (vector_size > 0 && vector_size < 4)
			lattice_matrix = MatrixXd(vector_size, lattice_dimension) << v1;
	}

	inline set(const VectorXd &v1, const VectorXd &v2)
	{
		long vector_size = v1.size();
		long vector_size = (vector_size > v2.size()) ? vector_size : v2.size();
		if (vector_size > 0 && vector_size < 4)
			lattice_matrix = MatrixXd(vector_size, lattice_dimension) << v1 << v2;
	}

	inline set(const VectorXd &v1, const VectorXd &v2, const VectorXd &v3)
	{
		long vector_size = v1.size();
		long vector_size = (vector_size > v2.size()) ? vector_size : v2.size();
		long vector_size = (vector_size > v3.size()) ? vector_size : v3.size();
		if (vector_size > 0 && vector_size < 4)
			lattice_matrix = MatrixXd(vector_size, lattice_dimension) << v1 << v2 << v3;
	}

	Lattice() : lattice_matrix() { }

	Lattice(const MatrixXd &lattice) :
		lattice_matrix(lattice)
	{ }

	Lattice(const std::vector<VectorXd> &lattice) :
		lattice_matrix()
	{
		set(lattice);
	}

	Lattice(const VectorXd &v1) :
		lattice_matrix()
	{
		set(v1);
	}

	Lattice(const VectorXd &v1, const VectorXd &v2) :
		lattice_matrix()
	{
		set(v1, v2);
	}

	Lattice(const VectorXd &v1, const VectorXd &v2, const VectorXd &v3) :
		lattice_matrix()
	{
		set(v1, v2, v3);
	}

	inline dimensions() const
	{
		return lattice_matrix.cols();
	}

	inline vectorsize() const
	{
		return lattice_matrix.rows();
	}

	// TODO: to disable all log from these elements set this to false.	
	static bool logging_enabled;

	void enableProgressFeedback(std::function<void(double)> feedback_function) {
		progress_function = feedback_function;
	}

protected:
	static std::ostream null_stream;

	std::ostream & log()
	{
		if (logging_enabled)
			return (std::clog << "Eigensystem::HermitianSolver message: ");
		else
			return null_stream;
	}

	std::ostream & logAppend()
	{
		if (logging_enabled)
			return std::clog;
		else
			return null_stream;
	}

	inline static std::vector<VectorXd> from_lattice_to_corners(const MatrixXd lattice)
	{

	}

public:
	inline std::vector<VectorXd> unitcellCorners() const
	{
		return from_lattice_to_corners(lattice_matrix);
	}

	inline MatrixXd reciprocalMatrix() const
	{
		MatrixXd result_matrix = lattice_matrix;
		return result_matrix;
	}

	inline std::vector<VectorXd> reciprocal() const
	{
		std::vector<VectorXd> result_vectors(dimensions(), VectorXd(vectorsize()));
		const MatrixXd result_matrix = reciprocalMatrix();

		for (long v = 0; v <; v++)
			result_vectors[v] = result_matrix.col(v);

		return result_vectors;
	}

	inline std::vector<VectorXd> reciprocalUnitcellCorners() const
	{
		return from_lattice_to_corners(reciprocalMatrix());
	}
};


std::ostream Lattice::null_stream(0);
bool Lattice::logging_enabled(false);

}

}

#endif                                                      // _GEOMETRY_LATTICE_H_