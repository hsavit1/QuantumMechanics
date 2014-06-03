/*
---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/

#ifndef EIGEN_BLOCKMATRIX_H
#define EIGEN_BLOCKMATRIX_H

namespace Eigen {


	template<typename _Derived>
	class BlockMatrixBase
	{
	public:
		typedef _Derived Derived;

		virtual inline Derived &derived() = 0;

		virtual inline const Derived &derived() const = 0;
	};

	template<typename _Scalar, int _Rows, int _Cols>
	class BlockMatrix : public BlockMatrixBase< BlockMatrix<_Scalar, _Rows, _Cols> >
	{
	public:
		typedef Matrix<_Scalar, _Rows, _Cols> Base;
		typedef BlockMatrix<_Scalar, _Rows, _Cols> ReturnType;
		typedef typename Block<Base> ReturnBlockType;
		typedef typename Block<const Base> ConstBlockType;
		typedef typename const Block<const Base> ConstReturnBlockType;

		typedef Base::Index Index;

		inline ReturnType &derived()
		{
			return *this;
		}

		inline const ReturnType &derived() const
		{
			return *this;
		}

	private:

		Base *base;
		bool owner;
		ArrayXi row_sizes, column_sizes;
		ArrayXi row_offsets, column_offsets;

		long block_rows_offset;
		long block_cols_offset;
		long block_rows_count;
		long block_cols_count;


	protected:
		inline static ArrayXi from_zero_cum_sum(const ArrayXi &vector)
		{
			ArrayXi result(vector.size());
			result[0] = 0;
			
			if (vector.size() > 1) // A cumulative sum (first,last,destination_start).
				std::partial_sum(&vector[0], &vector[vector.size() - 1], &result[1]);

			return result;
		}

	public:
		inline const VectorBlock<const ArrayXi> blockRowSizes() const		{ return row_sizes.segment(block_rows_offset, block_rows_count); }
		inline const VectorBlock<const ArrayXi> blockColSizes() const		{ return column_sizes.segment(block_cols_offset, block_cols_count); }

		inline const VectorBlock<const ArrayXi> blockRowOffsets() const		{ return row_offsets.segment(block_rows_offset, block_rows_count); }
		inline const VectorBlock<const ArrayXi> blockColOffsets() const		{ return column_offsets.segment(block_cols_offset, block_cols_count); }

		inline const long blockRows() const									{ return blockRowSizes().size(); }
		inline const long blockCols() const									{ return blockColSizes().size(); }

		inline const long blockRowSize(long i) const						{ i = (i >= 0) ? i : blockRows() + i; return blockRowSizes()[i]; }
		inline const long blockColSize(long i) const						{ i = (i >= 0) ? i : blockCols() + i; return blockColSizes()[i]; }

		inline const long blockRowOffset(long i) const						{ i = (i >= 0) ? i : blockRows() + i; return blockRowOffsets()[i]; }
		inline const long blockColOffset(long i) const						{ i = (i >= 0) ? i : blockCols() + i; return blockColOffsets()[i]; }

		inline long rows() const
		{
			return blockRowSizes().sum();
		}

		inline long cols() const
		{
			return blockColSizes().sum();
		}

		void setBlocks(const ArrayXi &isotropic_sizes)
		{
			if (!owner)
				return;

			row_sizes = isotropic_sizes;
			column_sizes = isotropic_sizes;
			row_offsets = from_zero_cum_sum(isotropic_sizes);
			column_offsets = row_offsets;

			block_rows_count = isotropic_sizes.size();
			block_cols_count = isotropic_sizes.size();
		}

		void setBlocks(const ArrayXi &rows, const ArrayXi &columns)
		{
			if (!owner)
				return;

			row_sizes = rows;
			column_sizes = columns;
			row_offsets = from_zero_cum_sum(rows);
			column_offsets = from_zero_cum_sum(columns);

			block_rows_count = rows.size();
			block_cols_count = columns.size();
		}

		void resetBlocks()
		{
			if (!owner)
				return;

			row_sizes = ArrayXi::Constant(1, base->rows());
			column_sizes = ArrayXi::Constant(1, base->cols());
			row_offsets = ArrayXi::Zero(1);
			column_offsets = ArrayXi::Zero(1);

			block_rows_offset = (0);
			block_cols_offset = (0);
			block_rows_count = (1);
			block_cols_count =(1);
		}
		
		inline ReturnBlockType matrix()
		{
			return ReturnBlockType(
				*base,
				blockRowOffset(0),
				blockColOffset(0),
				rows(),
				cols()
				);
		}

		inline const ConstReturnBlockType matrix() const
		{
			return ConstBlockType(
				*base,
				blockRowOffset(0),
				blockColOffset(0),
				rows(),
				cols()
				);
		}

		static inline ReturnType reference(const BlockMatrix& owner, long blockRowsOffset, long blockColsOffset, long blockRowsCount, long blockColsCount)
		{
			ReturnType result;

			if (result.owner && result.base)
			{
				delete result.base;
				result.owner = false;
				result.base = nullptr;
			}

			result.base = owner.base;
			result.owner = false;
			result.row_sizes = (owner.row_sizes);
			result.column_sizes = (owner.column_sizes);
			result.row_offsets = (owner.row_offsets);
			result.column_offsets = (owner.column_offsets);

			result.block_rows_offset = blockRowsOffset + owner.block_rows_offset;
			result.block_cols_offset = blockColsOffset + owner.block_cols_offset;
			result.block_rows_count = blockRowsCount;
			result.block_cols_count = blockColsCount;

			return result;
		}

		static inline ReturnType copy(const BlockMatrix& owner)
		{
			ReturnType result(owner.matrix());
			result.setBlocks(owner.blockRowSizes(), owner.blockColSizes());
			return result;
		}

		inline ReturnType &withBlocks(const BlockMatrix &other)
		{
			if (other.rows() == rows() && other.cols() == cols())
			{
				if (block_rows_offset == 0 && block_rows_count == row_sizes.size())
					row_sizes = (other.blockRowSizes());
				else
				{
					long total_block_count = row_sizes.size() - block_rows_count + other.blockRows();

					ArrayXi old_start = row_sizes.segment(0, block_rows_offset);
					ArrayXi old_end = row_sizes.segment(block_rows_offset + block_rows_count, row_sizes.size());

					block_rows_count = other.blockRows();
					row_sizes.resize(total_block_count);

					row_sizes.segment(0, block_rows_offset) = old_start;
					row_sizes.segment(block_rows_offset, block_rows_count) = other.blockRowSizes();
					row_sizes.segment(block_rows_offset + block_rows_count, row_sizes.size()) = old_end;
				}

				if (block_cols_offset == 0 && block_cols_count == column_sizes.size())
					column_sizes = (other.blockColSizes());
				else
				{
					long total_block_count = column_sizes.size() - block_cols_count + other.blockCols();

					ArrayXi old_start = column_sizes.segment(0, block_cols_offset);
					ArrayXi old_end = column_sizes.segment(block_cols_offset + block_cols_count, column_sizes.size());

					block_cols_count = other.blockCols();
					column_sizes.resize(total_block_count);

					column_sizes.segment(0, block_cols_offset) = old_start;
					column_sizes.segment(block_cols_offset, block_cols_count) = other.blockColSizes();
					column_sizes.segment(block_cols_offset + block_cols_count, column_sizes.size()) = old_end;
				}

				row_offsets = from_zero_cum_sum(row_sizes);
				column_offsets = from_zero_cum_sum(column_sizes);
			}

			return *this;
		}

		BlockMatrix(void) : 
			base(new Base),
			owner(true),
			row_sizes(ArrayXi::Zero(1)),
			column_sizes(ArrayXi::Zero(1)),
			row_offsets(ArrayXi::Zero(1)),
			column_offsets(ArrayXi::Zero(1)),
			block_rows_offset(0), 
			block_cols_offset(0),
			block_rows_count(1),
			block_cols_count(1)
		{ }

		BlockMatrix(const BlockMatrix& other) :
			base(				(other.owner) ? new Base(other.matrix())					: other.base				),
			owner(				 other.owner																			),
			row_sizes(			(other.owner) ? other.blockRowSizes()						: other.row_sizes			),
			column_sizes(		(other.owner) ? other.blockColSizes()						: other.column_sizes		),
			row_offsets(		(other.owner) ? from_zero_cum_sum(other.blockRowSizes())	: other.row_offsets			),
			column_offsets(		(other.owner) ? from_zero_cum_sum(other.blockColSizes())	: other.column_offsets		),
			block_rows_offset(	(other.owner) ? 0											: other.block_rows_offset	),
			block_cols_offset(	(other.owner) ? 0											: other.block_cols_offset	),
			block_rows_count(	(other.owner) ? other.blockRows()							: other.block_rows_count	),
			block_cols_count(	(other.owner) ? other.blockCols()							: other.block_cols_count	)
		{ }

		// This constructor allows you to construct MyVectorType from Eigen expressions
		template<typename OtherDerived>
		BlockMatrix(const MatrixBase<OtherDerived>& other) :
			base(new Base(other.derived())),
			owner(true),
			row_sizes(ArrayXi::Constant(1, other.rows())),
			column_sizes(ArrayXi::Constant(1, other.cols())),
			row_offsets(ArrayXi::Zero(1)),
			column_offsets(ArrayXi::Zero(1)),
			block_rows_offset(0),
			block_cols_offset(0),
			block_rows_count(1),
			block_cols_count(1)
		{ }

		~BlockMatrix() 
		{
			if (owner && base)
			{
				delete base;
				base = nullptr;
			}
			//else if (!owner)
				//std::cout << "A reference BlockMatrix has been destroyed." << std::endl;
		}

		void assignToReference(const BlockMatrix& other)
		{
			base = other.base;
			owner = other.owner;
			row_sizes = other.row_sizes;
			column_sizes = other.column_sizes;
			row_offsets = other.row_offsets;
			column_offsets = other.column_offsets;
			block_rows_offset = other.block_rows_offset;
			block_cols_offset = other.block_cols_offset;
			block_rows_count = other.block_rows_count;
			block_cols_count = other.block_cols_count;
		}

		ReturnType & operator= (const BlockMatrix& other)
		{
			if (owner && base->size() == 0)
			{
				assignToReference(other);
			}
			if (owner && (other.owner || other.base != base))
			{
				*base = other.matrix();

				row_sizes = other.blockRowSizes();
				column_sizes = other.blockColSizes();
				row_offsets = from_zero_cum_sum(other.blockRowSizes());
				column_offsets = from_zero_cum_sum(other.blockColSizes());

				block_rows_offset = 0;
				block_cols_offset = 0;
				block_rows_count = other.blockRows();
				block_cols_count = other.blockCols();
			}
			if (!owner && (blockRowSizes() == other.blockRowSizes()).all() && (blockColSizes() == other.blockColSizes()).all())
			{
				matrix() = other.matrix();
			}
			else if (owner && other.base == base)
				std::cout << "A reference BlockMatrix has been assign to an alias, not allowed." << std::endl;
			else if (!owner)
				std::cout << "A reference BlockMatrix has been assign to a block matrix of different sizes, not allowed." << std::endl;

			return *this;
		}

		// This method allows you to assign Eigen expressions to MyVectorType
		template<typename OtherDerived>
		ReturnType & operator= (const MatrixBase <OtherDerived>& other)
		{
			if (!owner && other.rows() == rows() && other.cols() == cols())
			{
				matrix() = other.derived();
			}
			else if (owner && other.rows() == rows() && other.cols() == cols() && other.rows() == base->rows() && other.cols() == base->cols())
			{
				*base = other.derived();
			}
			else if (owner)
			{
				*base = other.derived();
				resetBlocks();
			}
			return *this;
		}

		operator Base() const
		{
			return matrix();
		}

		inline bool isSquare(bool also_square_block_view = true) const
		{
			return rows() == cols() && (!also_square_block_view || blockRows() == blockCols());
		}

		inline ReturnBlockType block(Index blockRow, Index blockColumn)
		{
			blockRow = (blockRow >= 0) ? blockRow : blockRows() + blockRow;
			blockColumn = (blockColumn >= 0) ? blockColumn : blockCols() + blockColumn;
			return ReturnBlockType(*base, blockRowOffset(blockRow), blockColOffset(blockColumn), blockRowSize(blockRow), blockColSize(blockColumn));
		}

		inline ConstReturnBlockType block(Index blockRow, Index blockColumn) const
		{
			blockRow = (blockRow >= 0) ? blockRow : blockRows() + blockRow;
			blockColumn = (blockColumn >= 0) ? blockColumn : blockCols() + blockColumn;
			return ConstBlockType(*base, blockRowOffset(blockRow), blockColOffset(blockColumn), blockRowSize(blockRow), blockColSize(blockColumn));
		}

		inline BlockMatrix blocks(Index startBlockRow, Index startBlockColumn, Index blockRowCount, Index blockColCount)
		{
			startBlockRow = (startBlockRow >= 0) ? startBlockRow : blockRows() + startBlockRow;
			startBlockColumn = (startBlockColumn >= 0) ? startBlockColumn : blockCols() + startBlockColumn;
			return reference(
				*this,
				((blockRowCount > 0) ? startBlockRow : startBlockRow + blockRowCount),
				((blockColCount > 0) ? startBlockColumn : startBlockColumn + blockColCount),
				((blockRowCount > 0) ? blockRowCount : -blockRowCount),
				((blockColCount > 0) ? blockColCount : -blockColCount)
				);
		}

		inline const BlockMatrix blocks(Index startBlockRow, Index startBlockColumn, Index blockRowCount, Index blockColCount) const
		{
			startBlockRow = (startBlockRow >= 0) ? startBlockRow : blockRows() + startBlockRow;
			startBlockColumn = (startBlockColumn >= 0) ? startBlockColumn : blockCols() + startBlockColumn;
			return reference(
				*this,
				((blockRowCount > 0) ? startBlockRow : startBlockRow + blockRowCount),
				((blockColCount > 0) ? startBlockColumn : startBlockColumn + blockColCount),
				((blockRowCount > 0) ? blockRowCount : -blockRowCount),
				((blockColCount > 0) ? blockColCount : -blockColCount)
				);
		}

		/*
		Now we reimplement matrix function from the Eigen::Matrix class.
		*/

		void setZero()
		{
			matrix().setZero();
		}

		void setIdentity()
		{
			matrix().setIdentity();
		}

		template<typename LocalReturnType = ReturnType>
		EIGEN_STRONG_INLINE LocalReturnType asZero() const
		{
			return LocalReturnType::Zero(rows(), cols());
		}

		template<>
		EIGEN_STRONG_INLINE ReturnType asZero<ReturnType>() const
		{
			ReturnType result = ReturnType::copy(*this);
			result.setZero();
			return result;
		}

		template<typename LocalReturnType = ReturnType>
		EIGEN_STRONG_INLINE LocalReturnType asIdentity() const
		{
			return LocalReturnType::Identity(rows(), cols());
		}

		template<>
		EIGEN_STRONG_INLINE ReturnType asIdentity<ReturnType>() const
		{
			ReturnType result = ReturnType::copy(*this);
			result.setIdentity();
			return result;
		}

		const ConstBlockType::AdjointReturnType adjoint() const
		{
			return matrix().adjoint();
		}

		void adjointInPlace()
		{
			matrix().adjointInPlace();
		}

		Base::Scalar trace() const
		{
			return matrix().trace();
		}

		const internal::inverse_impl<ConstBlockType> inverse() const
		{
			return matrix().inverse();
		}

		template<typename OtherDerived>
		ReturnType& operator+=(const MatrixBase<OtherDerived>& other)
		{
			matrix() += other.derived();
			return *this;
		}

		template<typename OtherDerived>
		ReturnType& operator-=(const MatrixBase<OtherDerived>& other)
		{
			matrix() -= other.derived();
			return *this;
		}

		template<typename OtherDerived>
		const typename ProductReturnType<ConstBlockType, OtherDerived>::Type
			operator*(const MatrixBase<OtherDerived> &other) const
		{
				return matrix()*other.derived();
		}

		template<typename OtherDerived>
		const typename ProductReturnType<ConstBlockType, OtherDerived::ConstBlockType>::Type
			operator*(const BlockMatrixBase<OtherDerived> &other) const
		{
				return matrix()*other.derived().matrix();
		}

		template<typename OtherDerived>
		EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_difference_op<Base::Scalar>, const ConstBlockType, const OtherDerived>
			operator-(const MatrixBase<OtherDerived> &other) const
		{
				return CwiseBinaryOp<internal::scalar_difference_op<Base::Scalar>, const ConstBlockType, const OtherDerived>(matrix(), other.derived());
		}

		template<typename OtherDerived>
		EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_difference_op<Base::Scalar>, const ConstBlockType, const OtherDerived::ConstBlockType>
			operator-(const BlockMatrixBase<OtherDerived> &other) const
		{
				return CwiseBinaryOp<internal::scalar_difference_op<Base::Scalar>, const ConstBlockType, const OtherDerived::ConstBlockType>(matrix(), other.derived().matrix());
		}

		template<typename OtherDerived>
		inline ReturnType&
			operator*=(const EigenBase<OtherDerived> &other)
		{
				other.derived().applyThisOnTheRight(matrix());
				return *this;
		}

	};


	template<typename Derived, typename OtherDerived>
	EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_difference_op<Derived::Scalar>, const Derived, const OtherDerived::ConstBlockType>
	operator-(const MatrixBase<Derived> &first, const BlockMatrixBase<OtherDerived> &other)
	{
		return CwiseBinaryOp<internal::scalar_difference_op<Derived::Scalar>, const Derived, const OtherDerived::ConstBlockType>(first.derived(), other.derived().matrix());
	}

	template<typename Derived, typename OtherDerived>
	const typename ProductReturnType<Derived, OtherDerived::ConstBlockType>::Type
	operator*(const MatrixBase<Derived> &first, const BlockMatrixBase<OtherDerived> &other)
	{
		return first.derived()*other.derived().matrix();
	}

#define EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, Size, SizeSuffix)   \
	/** \ingroup matrixtypedefs */                                    \
	typedef BlockMatrix<Type, Size, Size> BlockMatrix##SizeSuffix##TypeSuffix;

#define EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, Size)         \
	/** \ingroup matrixtypedefs */                                    \
	typedef BlockMatrix<Type, Size, Dynamic> BlockMatrix##Size##X##TypeSuffix;  \
	/** \ingroup matrixtypedefs */                                    \
	typedef BlockMatrix<Type, Dynamic, Size> BlockMatrix##X##Size##TypeSuffix;

#define EIGEN_MAKE_TYPEDEFS_ALL_SIZES(Type, TypeSuffix) \
	EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 2, 2) \
	EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 3, 3) \
	EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 4, 4) \
	EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, Dynamic, X) \
	EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 2) \
	EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 3) \
	EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 4)

EIGEN_MAKE_TYPEDEFS_ALL_SIZES(int, i)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(float, f)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(double, d)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(std::complex<float>, cf)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(std::complex<double>, cd)

#undef EIGEN_MAKE_TYPEDEFS_ALL_SIZES
#undef EIGEN_MAKE_TYPEDEFS
#undef EIGEN_MAKE_FIXED_TYPEDEFS

};

#endif // EIGEN_BLOCKED_MATRIX_H