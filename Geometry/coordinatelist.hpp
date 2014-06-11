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

class CoordinateList;

template<typename DERIVED>
class CoordinateListBase
{
public:
	typedef DERIVED Derived;

	inline Derived &derived()
	{
		return (Derived)*this;
	}

	virtual inline CoordinateList<Derived::Scalar, Derived::Columns> coordinatelist() const = 0;

	virtual inline Derived::ListMatrix coordinatematrix() const = 0;
};

template<typename SCALAR, typename COLS>
class CoordinateList : public CoordinateListBase< CoordinateList<SCALAR, COLS> >
{
public:
	typedef SCALAR Scalar;
	typedef COLS Columns;
	typedef CoordinateListBase< CoordinateList<SCALAR, COLS> > Base;
	typedef Matrix<SCALAR, Dynamic, COLS> ListMatrix;
	typedef typename ListMatrix::RowXpr CoordinateRow;
	typedef typename ListMatrix::ConstRowXpr ConstCoordinateRow;

private:
	ListMatrix list;

public:
	virtual inline CoordinateList coordinatelist() const
	{
		return *this;
	}

	virtual inline ListMatrix coordinatematrix() const
	{
		return list;
	}

	CoordinateList() {}

	CoordinateList(const Base &base) :
		list(base.coordinatematrix())
	{ }

	CoordinateList(const CoordinateList &other) :
		list(other.list)
	{ }

	CoordinateList(const long &rows) :
		list(rows, COLS)
	{ }

	CoordinateList &operator = (const Base &base)
	{
		list = base.coordinatematrix();
	}

	CoordinateList &operator = (const CoordinateList &other)
	{
		list = other.list;
	}

	const long &size() const
	{
		return list.rows();
	}

	ConstCoordinateRow operator [] (const long &i) const
	{
		return list.row(i);
	}

	CoordinateRow operator [] (const long &i)
	{
		return list.row(i);
	}
};

typedef CoordinateList<long, 2> CoordinateList2i;
typedef CoordinateList<long, 3> CoordinateList3i;
typedef CoordinateList<long, 4> CoordinateList4i;
typedef CoordinateList<long, 5> CoordinateList5i;

typedef CoordinateList<double, 2> CoordinateList2d;
typedef CoordinateList<double, 3> CoordinateList3d;
typedef CoordinateList<double, 4> CoordinateList4d;
typedef CoordinateList<double, 5> CoordinateList5d;


template<typename COORDSLIST>
class CoordinateListMerge : public CoordinateListBase< CoordinateListMerge<COORDSLIST> >
{
public:
	typedef COORDSLIST::Scalar Scalar;
	typedef COORDSLIST::Columns Columns;
	typedef CoordinateListBase< CoordinateListMerge<COORDSLIST>> Base;
	typedef COORDSLIST CoordsList;
	typedef CoordsList::ListMatrix ListMatrix;

private:
	std::vector<const CoordsList *> mergelist;

public:
	const long &size() const
	{
		const long size = mergelist.size();
		long rows = 0;

		for (long i = 0; i < size; i++)
		{
			rows += sizes[i];
		}

		return rows;
	}

	virtual inline ListMatrix coordinatematrix() const
	{
		const long size = mergelist.size();
		long rows = 0;
		std::vector<long> offsets(size);
		std::vector<long> sizes(size);

		for (long i = 0; i < size; i++)
		{
			offsets[i] = rows;
			sizes[i] = mergelist[i]->size();
			rows += sizes[i];
		}

		ListMatrix result(rows, Columns);

		for (long i = 0; i < mergelist.size(); i++)
			result.block(offsets[i], 0, sizes[i], Columns) = mergelist[i]->coordinatematrix();

		return result;
	}

	virtual inline CoordinateList<Scalar, Columns> coordinatelist() const
	{
		return CoordinateList<Scalar, Columns>(coordinatematrix());
	}


protected:
	static inline std::vector<const CoordsList *> tolist(const CoordsList &l1, const CoordsList &l2)
	{
		std::vector<const CoordsList *> result(2);

		result[0] = &l1;
		result[1] = &l2;

		return result;
	}

	static inline std::vector<const CoordsList *> tolist(const CoordsList &l1, const CoordsList &l2, const CoordsList &l3)
	{
		std::vector<const CoordsList *> result(3);

		result[0] = &l1;
		result[1] = &l2;
		result[2] = &l3;

		return result;
	}

	static inline std::vector<const CoordsList *> tolist(const CoordsList &l1, const CoordsList &l2, const CoordsList &l3, const CoordsList &l4)
	{
		std::vector<const CoordsList *> result(4);

		result[0] = &l1;
		result[1] = &l2;
		result[2] = &l3;
		result[3] = &l4;

		return result;
	}

	static inline std::vector<const CoordsList *> tolist(const CoordsList &l1, const CoordsList &l2, const CoordsList &l3, const CoordsList &l4, const CoordsList &l5)
	{
		std::vector<const CoordsList *> result(5);

		result[0] = &l1;
		result[1] = &l2;
		result[2] = &l3;
		result[3] = &l4;
		result[4] = &l5;

		return result;
	}

public:
	CoordinateListMerge(const std::vector<const CoordsList *> &l) :
		mergelist(l)
	{ }

	CoordinateListMerge(const CoordsList &l1, const CoordsList &l2) :
		mergelist(tolist(l1, l2))
	{ }

	CoordinateListMerge(const CoordsList &l1, const CoordsList &l2, const CoordsList &l3) :
		mergelist(tolist(l1, l2, l3))
	{ }

	CoordinateListMerge(const CoordsList &l1, const CoordsList &l2, const CoordsList &l3, const CoordsList &l4) :
		mergelist(tolist(l1, l2, l3, l4))
	{ }

	CoordinateListMerge(const CoordsList &l1, const CoordsList &l2, const CoordsList &l3, const CoordsList &l4, const CoordsList &l5) :
		mergelist(tolist(l1, l2, l3, l4, l5))
	{ }

	CoordinateListMerge(const long &N) :
		mergelist(N)
	{ }

	void setPart(const long &i, const CoordsList &l)
	{
		mergelist[i] = l;
	}
};

template<typename COORDSLIST>
CoordinateListMerge<COORDSLIST> merge(const COORDSLIST &l1, const COORDSLIST &l2)
{
	return CoordinateListMerge<COORDSLIST>(l1, l2);
}

template<typename COORDSLIST>
CoordinateListMerge<COORDSLIST> merge(const COORDSLIST &l1, const COORDSLIST &l2, const COORDSLIST &l3)
{
	return CoordinateListMerge<COORDSLIST>(l1, l2, l3);
}

template<typename COORDSLIST>
CoordinateListMerge<COORDSLIST> merge(const COORDSLIST &l1, const COORDSLIST &l2, const COORDSLIST &l3, const COORDSLIST &l4)
{
	return CoordinateListMerge<COORDSLIST>(l1, l2, l3, l4);
}

template<typename COORDSLIST>
CoordinateListMerge<COORDSLIST> merge(const COORDSLIST &l1, const COORDSLIST &l2, const COORDSLIST &l3, const COORDSLIST &l4, const COORDSLIST &l5)
{
	return CoordinateListMerge<COORDSLIST>(l1, l2, l3, l4, l5);
}


template<typename COORDSLIST>
class CoordinateListRepeat : public CoordinateListBase< CoordinateListRepeat<COORDSLIST> >
{
public:
	typedef COORDSLIST::Scalar Scalar;
	typedef COORDSLIST::Columns Columns;
	typedef CoordinateListBase< CoordinateListRepeat<COORDSLIST>> Base;
	typedef COORDSLIST CoordsList;
	typedef CoordsList::ListMatrix ListMatrix;

private:
	const CoordsList &list;

	ListMatrix displacementlist;

public:
	const long &size() const
	{
		return list.size() * displacementlist.rows();
	}

	virtual inline ListMatrix coordinatematrix() const
	{
		const long rows = list.size();
		const long size = displacementlist.rows();
		const long allrows = rows * size;

		ListMatrix result(allrows, Columns);

		for (long i = 0; i < size; i++)
		{
			result.block(rows * i, 0, rows, Columns) = list->coordinatematrix();
			result.block(rows * i, 0, rows, Columns).rowwise() += displacementlist.row(i);
		}

		return result;
	}

	virtual inline CoordinateList<Scalar, Columns> coordinatelist() const
	{
		return CoordinateList<Scalar, Columns>(coordinatematrix());
	}

protected:
	static inline ListMatrix deleteDuplicateRows(const ListMatrix &fulllist)
	{
		const long rows = fulllist.rows();
		std::vector<long> goodrows;
		goodrows.reserve(rows);

		for (long i = 1; i < rows; i++)
		{
			if ( !(fulllist.topRows(i) == fulllist.row(i)).any() )
				goodrows.append(i);
		}

		const long new_rows = goodrows.size();

		ListMatrix result(new_rows, Columns);

		for (long i = 1; i < new_rows; i++)
		{
			result.row(i) = fulllist.row(goodrows[i]);
		}

		return result;
	}

	static inline ListMatrix makedisplacements(const Matrix<Scalar, 1, Columns> &vectors, const Matrix<long, 1, 1> &repeats)
	{
		ListMatrix result(repeats[0], Columns);

		for (int i = 0; i < repeats[0]; i++)
		{
			result.row() = vectors[0] * i;
		}

		return result;
	}

	static inline ListMatrix makedisplacements(const Matrix<Scalar, 2, Columns> &vectors, const Matrix<long, 2, 1> &repeats)
	{
		ListMatrix result(repeats[0] * repeats[1], Columns);

		for (int i = 0; i < repeats[0]; i++)
		{
			for (int j = 0; j < repeats[1]; j++)
			{
				result.row() = vectors[0] * i + vectors[1] * j;
			}
		}

		return deleteDuplicateRows(result);
	}

	static inline ListMatrix makedisplacements(const Matrix<Scalar, 3, Columns> &vectors, const Matrix<long, 3, 1> &repeats)
	{
		ListMatrix result(repeats[0] * repeats[1] * repeats[2], Columns);

		for (int i = 0; i < repeats[0]; i++)
		{
			for (int j = 0; j < repeats[1]; j++)
			{
				for (int h = 0; h < repeats[2]; h++)
				{
					result.row() = vectors[0] * i + vectors[1] * j + vectors[2] * h;
				}
			}
		}

		return deleteDuplicateRows(result);
	}

public:
	CoordinateListRepeat(const COORDSLIST &l, const Matrix<Scalar, 1, Columns> &vectors, const Matrix<long, 1, 1> &repeats) :
		list(l),
		displacementlist(makedisplacements(vectors, repeats))
	{ }

	CoordinateListRepeat(const COORDSLIST &l, const Matrix<Scalar, 2, Columns> &vectors, const Matrix<long, 2, 1> &repeats) :
		list(l),
		displacementlist(makedisplacements(vectors, repeats))
	{ }

	CoordinateListRepeat(const COORDSLIST &l, const Matrix<Scalar, 3, Columns> &vectors, const Matrix<long, 3, 1> &repeats) :
		list(l),
		displacementlist(makedisplacements(vectors, repeats))
	{ }
};

template<typename COORDSLIST>
CoordinateListRepeat<COORDSLIST> repeat(const COORDSLIST &list, const Matrix<COORDSLIST::Scalar, 1, COORDSLIST::Columns> &vectors, const Matrix<long, 1, 1> &repeats)
{
	return CoordinateListRepeat<COORDSLIST>(list, vectors, repeats);
}

template<typename COORDSLIST>
CoordinateListRepeat<COORDSLIST> repeat(const COORDSLIST &list, const Matrix<COORDSLIST::Scalar, 2, COORDSLIST::Columns> &vectors, const Matrix<long, 2, 1> &repeats)
{
	return CoordinateListRepeat<COORDSLIST>(list, vectors, repeats);
}

template<typename COORDSLIST>
CoordinateListRepeat<COORDSLIST> repeat(const COORDSLIST &list, const Matrix<COORDSLIST::Scalar, 3, COORDSLIST::Columns> &vectors, const Matrix<long, 3, 1> &repeats)
{
	return CoordinateListRepeat<COORDSLIST>(list, vectors, repeats);
}

#include <algorithm>    // std::sort

template<typename COORDSLIST, typename SORTFUNC>
class CoordinateListSort : public CoordinateListBase< CoordinateListSort<COORDSLIST, SORTFUNC> >
{
public:
	typedef COORDSLIST::Scalar Scalar;
	typedef COORDSLIST::Columns Columns;
	typedef CoordinateListBase< CoordinateListSort<COORDSLIST, SORTFUNC>> Base;
	typedef COORDSLIST CoordsList;
	typedef SORTFUNC SortFunction;
	typedef CoordsList::ListMatrix ListMatrix;

private:
	const CoordsList &list;
	ListMatrix coords;

	SortFunction sort_function;

protected:
	inline std::function<bool(VectorXi::Index, VectorXi::Index)> get_indices_sorting_function() const
	{
		coords = list.coordinatematrix();

		return[&](VectorXi::Index i, VectorXi::Index j) {
			return sort_function(coords.row(i), coords.row(j));
		};
	}

public:
	const long &size() const
	{
		return list.size();
	}

	virtual inline ListMatrix coordinatematrix() const
	{
		const long rows = list.size();
		VectorXi indices = VectorXi::LinSpaced(Sequential, 0, rows, 1);

		std::sort(indices.data(), indices.data() + rows, get_indices_sorting_function());

		ListMatrix result(rows, Columns);

		for (long i = 0; i < size; i++)
		{
			result.row(i) = coords.row(indices[i]);
		}

		return result;
	}

	virtual inline CoordinateList<Scalar, Columns> coordinatelist() const
	{
		return CoordinateList<Scalar, Columns>(coordinatematrix());
	}

public:
	CoordinateListSort(const COORDSLIST &list, const SortFunction &func) :
		list(list),
		sort_function(func)
	{ }
};

template<typename COORDSLIST, typename SORTFUNC>
CoordinateListSort<COORDSLIST, SORTFUNC> sort(const COORDSLIST &list, const SORTFUNC &func)
{
	return CoordinateListSort<COORDSLIST, SORTFUNC>(list, func);
}

}

}

#endif                                                      // _GEOMETRY_LATTICE_H_