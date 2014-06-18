/*
Header file for QuantumMechanics::Geometry::Field:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/
#ifndef _GEOMETRY_FIELD_H_
#define _GEOMETRY_FIELD_H_

#include <Math/Dense>

namespace QuantumMechanics {

namespace Geometry {
	
template<typename SCALAR, long DIMENSIONS>
class Field {

	typedef SCALAR Scalar;
	typedef DIMENSIONS Dimensions;
	typedef std::vector<Scalar> ScalarVector;
	typedef std::vector< ScalarVector > FieldBase;
	typedef Matrix<MatrixXi::Index, Dimensions, 1> VectorDi;

private:
	FieldBase base;
	VectorDi base_dimensions;

public:
	inline const VectorDi &dimensions() const
	{
		return base_dimensions;
	}

	inline const MatrixXi::Index &dimensions(const long &d) const
	{
		return base_dimensions[d];
	}

	inline const FieldBase &data() const
	{
		return base;
	}

	Field() { }

	Field(const Field &f) :
		base(f.data()),
		base_dimensions(f.dimensions())
	{ }

	inline void setEmpty()
	{
		std::fill(base.begin(), base.end(), ScalarVector());
	}
};

template<typename SCALAR>
class Field<SCALAR, 1>{

	typedef SCALAR Scalar;
	typedef 1 Dimensions;
	typedef std::vector<Scalar> ScalarVector;
	typedef std::vector< ScalarVector > FieldBase;
	typedef Matrix<MatrixXi::Index, Dimensions, 1> VectorDi;

private:
	FieldBase base;
	VectorDi base_dimensions;

public:
	inline const VectorDi &dimensions() const
	{
		return base_dimensions;
	}

	inline const MatrixXi::Index &dimensions(const long &d) const
	{
		return base_dimensions[d];
	}

	inline const FieldBase &data() const
	{
		return base;
	}

	Field() { }

	Field(const Field &f) :
		base(f.data()),
		base_dimensions(f.dimensions())
	{ }

	Field(const long &d1) :
		base(d1),
		base_dimensions(d1)
	{ }

	inline const ScalarVector &at(const long &i) const
	{
		return base[i];
	}

	inline ScalarVector &at(const long &i)
	{
		return base[i];
	}

	inline void resize(const long &d1)
	{
		base.resize(d1);
		base_dimensions = VectorDi(d1);
	}

	inline void setEmpty()
	{
		std::fill(base.begin(), base.end(), ScalarVector());
	}
};

template<typename SCALAR>
class Field<SCALAR, 2>{

	typedef SCALAR Scalar;
	typedef 2 Dimensions;
	typedef std::vector<Scalar> ScalarVector;
	typedef std::vector< ScalarVector > FieldBase;
	typedef Matrix<MatrixXi::Index, Dimensions, 1> VectorDi;

private:
	FieldBase base;
	VectorDi base_dimensions;

public:
	inline const VectorDi &dimensions() const
	{
		return base_dimensions;
	}

	inline const MatrixXi::Index &dimensions(const long &d) const
	{
		return base_dimensions[d];
	}

	inline const FieldBase &data() const
	{
		return base;
	}

	Field() { }

	Field(const Field &f) :
		base(f.data()),
		base_dimensions(f.dimensions())
	{ }

	Field(const long &d1, const long &d2) :
		base(d1 * d2),
		base_dimensions(d1, d2)
	{ }

	inline const ScalarVector &at(const long &i, const long &j) const
	{
		return base[i + j * base_dimensions[0]];
	}

	inline ScalarVector &at(const long &i, const long &j)
	{
		return base[i + j * base_dimensions[0]];
	}

	inline void resize(const long &d1, const long &d2)
	{
		base.resize(d1 * d2);
		base_dimensions = VectorDi(d1, d2);
	}

	inline void setEmpty()
	{
		std::fill(base.begin(), base.end(), ScalarVector());
	}
};

template<typename SCALAR>
class Field<SCALAR, 3>{

	typedef SCALAR Scalar;
	typedef 3 Dimensions;
	typedef std::vector<Scalar> ScalarVector;
	typedef std::vector< ScalarVector > FieldBase;
	typedef Matrix<MatrixXi::Index, Dimensions, 1> VectorDi;

private:
	FieldBase base;
	VectorDi base_dimensions;

public:
	inline const VectorDi &dimensions() const
	{
		return base_dimensions;
	}

	inline const MatrixXi::Index &dimensions(const long &d) const
	{
		return base_dimensions[d];
	}

	inline const FieldBase &data() const
	{
		return base;
	}

	Field() { }

	Field(const Field &f) :
		base(f.data()),
		base_dimensions(f.dimensions())
	{ }

	Field(const long &d1, const long &d2, const long &d3) :
		base(d1 * d2 * d3),
		base_dimensions(d1, d2, d3)
	{ }

	inline const ScalarVector &at(const long &i, const long &j, const long &h) const
	{
		return base[i + j * base_dimensions[0] + h * base_dimensions[1]];
	}

	inline ScalarVector &at(const long &i, const long &j, const long &h)
	{
		return base[i + j * base_dimensions[0] + h * base_dimensions[1]];
	}

	inline void resize(const long &d1, const long &d2, const long &d3)
	{
		base.resize(d1 * d2 * d3);
		base_dimensions = VectorDi(d1, d2, d3);
	}

	inline void setEmpty()
	{
		std::fill(base.begin(), base.end(), ScalarVector());
	}
};

template<typename SCALAR>
inline Field<SCALAR, 1> gatherNearestNeighbors(const Field<SCALAR, 1> &f, const long &nn = 1)
{
	Field<SCALAR, 1> result(f);

	for (long I = 0; I < f.dimensions(0); I++)
	{
		result.at(I).clear();
		long size = 0;

		for (long i = min(0, I - nn); i < max(I + nn + 1, f.dimensions(0)); i++)
		{
			size += f.at(i).size();
		}

		result.at(I).reserve(size);

		for (long i = min(0, I - nn); i < max(I + nn + 1, f.dimensions(0)); i++)
		{
			result.at(I).append(f.at(i));
		}
	}

	return result;
}

template<typename SCALAR>
inline Field<SCALAR, 2> gatherNearestNeighbors(const Field<SCALAR, 2> &f, const long &nn = 1)
{
	Field<SCALAR, 1> result(f);

	for (long I = 0; I < f.dimensions(0); I++)
	{
		for (long J = 0; J < f.dimensions(1); J++)
		{
			result.at(I, J).clear();
			long size = 0;

			for (long i = min(0, I - nn); i < max(I + nn + 1, f.dimensions(0)); i++)
			{
				for (long j = min(0, J - nn); j < max(J + nn + 1, f.dimensions(1)); j++)
				{
					size += f.at(i, j).size();
				}
			}

			result.at(I, J).reserve(size);

			for (long i = min(0, I - nn); i < max(I + nn + 1, f.dimensions(0)); i++)
			{
				for (long j = min(0, J - nn); j < max(J + nn + 1, f.dimensions(1)); j++)
				{
					result.at(I, J).append(f.at(i, j));
				}
			}
		}
	}

	return result;
}

template<typename SCALAR>
inline Field<SCALAR, 3> gatherNearestNeighbors(const Field<SCALAR, 3> &f, const long &nn = 1)
{
	Field<SCALAR, 1> result(f);

	for (long I = 0; I < f.dimensions(0); I++)
	{
		for (long J = 0; J < f.dimensions(1); J++)
		{
			for (long H = 0; H < f.dimensions(2); H++)
			{

				result.at(I, J, H).clear();
				long size = 0;

				for (long i = min(0, I - nn); i < max(I + nn + 1, f.dimensions(0)); i++)
				{
					for (long j = min(0, J - nn); j < max(J + nn + 1, f.dimensions(1)); j++)
					{
						for (long h = min(0, H - nn); h < max(H + nn + 1, f.dimensions(2)); h++)
						{
							size += f.at(i, j, h).size();
						}
					}
				}

				result.at(I, J, H).reserve(size);

				for (long i = min(0, I - nn); i < max(I + nn + 1, f.dimensions(0)); i++)
				{
					for (long j = min(0, J - nn); j < max(J + nn + 1, f.dimensions(1)); j++)
					{
						for (long h = min(0, H - nn); h < max(H + nn + 1, f.dimensions(2)); h++)
						{
							result.at(I, J, H).append(f.at(i, j, h));
						}
					}
				}
			}
		}
	}

	return result;
}

}

}

#endif                                                      // _GEOMETRY_FIELD_H_