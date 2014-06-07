/*
Header file for QuantumMechanics::Geometry::Lattice:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/
#ifndef _GEOMETRY_LATTICE_H_
#define _GEOMETRY_LATTICE_H_

#include <Math/Dense>

namespace QuantumMechanics {

namespace Geometry {

class Lattice {

private:
	MatrixXd lattice_matrix;
	MatrixXd reciprocal_lattice_matrix;

protected:
	inline static MatrixXd create_lattice_matrix(const Vector2d &v1, const Vector2d &v2)
	{
		return MatrixXd(2, 2) << v1 << v2;
	}

	inline static MatrixXd create_lattice_matrix(const Vector3d &v1, const Vector3d &v2)
	{
		return MatrixXd(3, 2) << v1 << v2;
	}

	inline static MatrixXd create_lattice_matrix(const Vector3d &v1, const Vector3d &v2, const Vector3d &v3)
	{
		return MatrixXd(3, 3) << v1 << v2 << v3;
	}

	inline static MatrixXd calculate_reciprocal_matrix(const MatrixXd lattice)
	{
		const long dimensions = lattice.cols();
		const long v_length = lattice.rows();

		MatrixXd result;

		if (dimensions == 1)
		{
			result = 2 * M_PI * lattice.col(0) / lattice.col(0).squaredNorm();
		}
		else if (dimensions == v_length)
		{
			result = 2 * M_PI * lattice.inverse().transpose();
		}
		else if (dimensions == 2 && v_length == 3)
		{
			MatrixXd tmp_lattice = MatrixXd(3, 3) << lattice << lattice.col(0).cross(lattice.col(1)).normalized();

			result = 2 * M_PI * tmp_lattice.inverse().transpose();
		}
		else
		{
			//Error
		}

		return result;
	}

public:
	inline void set(const Matrix<double, 1, 1> &v)
	{
		lattice_matrix = v;
	}

	inline void set(const Vector2d &v)
	{
		lattice_matrix = v;
	}

	inline void set(const Vector3d &v)
	{
		lattice_matrix = v;
	}

	inline void set(const Vector2d &v1, const Vector2d &v2)
	{
		lattice_matrix = create_lattice_matrix(v1, v2);
	}

	inline void set(const Vector3d &v1, const Vector3d &v2)
	{
		lattice_matrix = create_lattice_matrix(v1, v2);
	}

	inline void set(const Vector3d &v1, const Vector3d &v2, const Vector3d &v3)
	{
		lattice_matrix = create_lattice_matrix(v1, v2, v3);
	}

	inline void set(const Matrix2d &m)
	{
		lattice_matrix = m;
	}

	inline void set(const Matrix<double, 3, 2> &m)
	{
		lattice_matrix = m;
	}

	inline void set(const Matrix3d &m)
	{
		lattice_matrix = m;
	}

	Lattice() { }

	Lattice(const Matrix<double, 1, 1> &v) :
		lattice_matrix(v),
		reciprocal_lattice_matrix(calculate_reciprocal_matrix(v))
	{ }

	Lattice(const Vector2d &v) :
		lattice_matrix(v),
		reciprocal_lattice_matrix(calculate_reciprocal_matrix(v))
	{ }

	Lattice(const Vector3d &v) :
		lattice_matrix(v),
		reciprocal_lattice_matrix(calculate_reciprocal_matrix(v))
	{ }

	Lattice(const Vector2d &v1, const Vector2d &v2) :
		lattice_matrix(create_lattice_matrix(v1, v2)),
		reciprocal_lattice_matrix(calculate_reciprocal_matrix(create_lattice_matrix(v1, v2)))
	{ }

	Lattice(const Vector3d &v1, const Vector3d &v2) :
		lattice_matrix(create_lattice_matrix(v1, v2)),
		reciprocal_lattice_matrix(calculate_reciprocal_matrix(create_lattice_matrix(v1, v2)))
	{ }

	Lattice(const Vector3d &v1, const Vector3d &v2, const Vector3d &v3) :
		lattice_matrix(create_lattice_matrix(v1, v2, v3)),
		reciprocal_lattice_matrix(calculate_reciprocal_matrix(create_lattice_matrix(v1, v2, v3)))
	{ }

	Lattice(const Matrix2d &m) :
		lattice_matrix(m),
		reciprocal_lattice_matrix(calculate_reciprocal_matrix(m))
	{ }

	Lattice(const Matrix<double, 3, 2> &m) :
		lattice_matrix(m),
		reciprocal_lattice_matrix(calculate_reciprocal_matrix(m))
	{ }

	Lattice(const Matrix3d &m) :
		lattice_matrix(m),
		reciprocal_lattice_matrix(calculate_reciprocal_matrix(m))
	{ }

	inline long dimensions() const
	{
		return lattice_matrix.cols();
	}

	inline long vectorsize() const
	{
		return lattice_matrix.rows();
	}

	inline const MatrixXd &latticeMatrix() const
	{
		return lattice_matrix;
	}
	
	inline const MatrixXd &reciprocalMatrix() const
	{
		return reciprocal_lattice_matrix;
	}
};

}

}

#endif                                                      // _GEOMETRY_LATTICE_H_