#include <mathimf.h>
#include <numeric>
#include "MKL_support_custom.hpp"
#include "range.hpp"
#define EIGEN_DENSEBASE_PLUGIN "denseplugin.hpp"
#define EIGEN_MATRIXBASE_PLUGIN "matrixplugin.hpp"
#include <Eigen/Dense>
#include "blockmatrix.hpp"


namespace Eigen
{

#define MATH_MAKE_MATRIX1_TYPEDEFS(Type, TypeSuffix)   \
	/** \ingroup matrixtypedefs */                                    \
	typedef Matrix<Type, 1, 1> Matrix1##TypeSuffix;  \
	/** \ingroup matrixtypedefs */                                    \
	typedef Matrix<Type, 1, 1>    Vector1##TypeSuffix;  \
	/** \ingroup matrixtypedefs */                              \
	typedef Matrix<Type, 1, 1>    RowVector1##TypeSuffix;

#define MATH_MAKE_FIXED_MATRIX1_TYPEDEFS(Type, TypeSuffix)         \
	/** \ingroup matrixtypedefs */                                    \
	typedef Matrix<Type, 1, Dynamic> Matrix1##X##TypeSuffix;  \
	/** \ingroup matrixtypedefs */                                    \
	typedef Matrix<Type, Dynamic, 1> Matrix##X1##TypeSuffix;

#define MATH_MAKE_ARRAY1_TYPEDEFS(Type, TypeSuffix)   \
	/** \ingroup matrixtypedefs */                                    \
	typedef Array<Type, 1, 1> Array11##TypeSuffix;

#define MATH_MAKE_FIXED_ARRAY1_TYPEDEFS(Type, TypeSuffix)         \
	/** \ingroup matrixtypedefs */                                    \
	typedef Array<Type, 1, Dynamic> Array1##X##TypeSuffix;  \
	/** \ingroup matrixtypedefs */                                    \
	typedef Array<Type, Dynamic, 1> Array##X1##TypeSuffix;

#define MATH_MAKE_TYPEDEFS_ALL_MATRIX1(Type, TypeSuffix) \
	MATH_MAKE_MATRIX1_TYPEDEFS(Type, TypeSuffix) \
	MATH_MAKE_FIXED_MATRIX1_TYPEDEFS(Type, TypeSuffix) \
	MATH_MAKE_ARRAY1_TYPEDEFS(Type, TypeSuffix) \
	MATH_MAKE_FIXED_ARRAY1_TYPEDEFS(Type, TypeSuffix)

	MATH_MAKE_TYPEDEFS_ALL_MATRIX1(int, i)
		MATH_MAKE_TYPEDEFS_ALL_MATRIX1(float, f)
		MATH_MAKE_TYPEDEFS_ALL_MATRIX1(double, d)
		MATH_MAKE_TYPEDEFS_ALL_MATRIX1(std::complex<float>, cf)
		MATH_MAKE_TYPEDEFS_ALL_MATRIX1(std::complex<double>, cd)

#undef MATH_MAKE_MATRIX1_TYPEDEFS
#undef MATH_MAKE_FIXED_MATRIX1_TYPEDEFS
#undef MATH_MAKE_ARRAY1_TYPEDEFS
#undef MATH_MAKE_FIXED_ARRAY1_TYPEDEFS
#undef MATH_MAKE_TYPEDEFS_ALL_MATRIX1

}	// namespace Math

namespace Math = Eigen;

using Math::Array;
using Math::Matrix;

#define MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, SizeSuffix) \
	using Math::Array##SizeSuffix##SizeSuffix##TypeSuffix; \
	using Math::Array##SizeSuffix##TypeSuffix; \
	using Math::Array##SizeSuffix##X##TypeSuffix; \
	using Math::Array##X##SizeSuffix##TypeSuffix;

#define MATH_USING_MATRIX_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, SizeSuffix) \
	using Math::Matrix##SizeSuffix##TypeSuffix; \
	using Math::BlockMatrix##SizeSuffix##TypeSuffix; \
	using Math::Vector##SizeSuffix##TypeSuffix; \
	using Math::RowVector##SizeSuffix##TypeSuffix;

#define MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE(TypeSuffix) \
	MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 1) \
	MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 2) \
	MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 3) \
	MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 4) \
	using Math::ArrayX##TypeSuffix; \
	using Math::ArrayXX##TypeSuffix; \
	MATH_USING_MATRIX_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 1) \
	MATH_USING_MATRIX_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 2) \
	MATH_USING_MATRIX_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 3) \
	MATH_USING_MATRIX_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, 4) \
	MATH_USING_MATRIX_TYPEDEFS_FOR_TYPE_AND_SIZE(TypeSuffix, X) \

#define MATH_USING_ARRAY_TYPEDEFS \
	MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE(i) \
	MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE(f) \
	MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE(d) \
	MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE(cf) \
	MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE(cd)

MATH_USING_ARRAY_TYPEDEFS

#undef MATH_USING_ARRAY_TYPEDEFS_FOR_TYPE_AND_SIZE
#undef MATH_USING_MATRIX_TYPEDEFS_FOR_TYPE_AND_SIZE
#undef MATH_USING_ARRAY_TYPEDEFS