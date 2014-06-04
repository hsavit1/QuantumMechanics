/*
---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/

typedef Matrix<RealScalar, Dynamic, 1> HermitianEigenvaluesReturnType;
typedef HermitianEigenvaluesReturnType HermitianEigenvaluesType;

inline HermitianEigenvaluesReturnType hermitianEigenvalues(range r) const
{
	typedef typename Derived::PlainObject PlainObject;
	PlainObject M(derived());

	r.fit_indices_to_size(M.rows());

	const char range_token = (r.type_value == range::full_range) ? 'A' : ((r.type_value == range::value_range) ? 'V' : 'I');
	const char job_token = 'N'; // Normal (only eigenvalues).
	const char upper_lower_token = 'U';

	
	// Problem dimensions:
	const lapack_int major_dim_order = (IsRowMajor) ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
	const lapack_int major_dim_length = (IsRowMajor) ? M.rows() : M.cols();
	const lapack_int size_vectors = (range_token == 'I') ? r.end_index - r.begin_index + 1 : major_dim_length;
	lapack_int value_count;

	HermitianEigenvaluesType values = Matrix<Scalar, Dynamic, 1>(size_vectors);

	// Additional eigenvector dimensions:
	MKL_Complex16* vectors = nullptr;
	lapack_int lead_dim_vectors = 1; // must be >= 1, or an error occurs. Only used for vector calculation.
	lapack_int* vectors_suppliements = new lapack_int[2 * major_dim_length];

	lapack_int info = LAPACKE_zheevr(
		major_dim_order,					// param. 0
		job_token,							// param. 1
		range_token,						// param. 2
		upper_lower_token,					// param. 3
		major_dim_length,					// param. 4
		M.data(),							// param. 5
		major_dim_length,					// param. 6
		r.lowest_value,						// param. 7
		r.highest_value,					// param. 8
		r.begin_index + 1,					// param. 9
		r.end_index + 1,					// param. 10
		0.,									// param. 11
		&value_count,						// param. 12
		values.data(),						// param. 13
		vectors,							// param. 14
		lead_dim_vectors,					// param. 15
		vectors_suppliements				// param. 16
		);

	delete[] vectors_suppliements;

	if (info != 0) {
		values.fill(nan(0));
	}

	if (value_count < lead_dim_vectors)
		values.conservativeResize(value_count);

	return values;
}

typedef Matrix<std::complex<RealScalar>, internal::traits<Derived>::ColsAtCompileTime, Dynamic> EigenvectorsType;
typedef std::pair<HermitianEigenvaluesType, EigenvectorsType> EigenvectorsReturnType;


inline EigenvectorsReturnType hermitianEigenvectors(range r) const
{
	typedef typename Derived::PlainObject PlainObject;
	PlainObject M(derived());

	r.fit_indices_to_size(M.rows());

	const char range_token = (r.type_value == range::full_range) ? 'A' : ((r.type_value == range::value_range) ? 'V' : 'I');
	const char job_token = 'V';
	const char upper_lower_token = 'U';

	const lapack_int lead_dim_vectors = M.rows();
	const lapack_int size_vectors = (range_token == 'I') ? r.end_index - r.begin_index + 1 : lead_dim_vectors;

	EigenvectorsReturnType result;
	HermitianEigenvaluesType &values = result.first;
	EigenvectorsType &eigenvectors = result.second;
	values = Matrix<RealScalar, Dynamic, 1>(size_vectors);
	eigenvectors = Matrix<Scalar, Dynamic, Dynamic>(lead_dim_vectors, size_vectors);

	// Problem dimensions:
	const lapack_int major_dim_order = (IsRowMajor) ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
	const lapack_int major_dim_length = (IsRowMajor) ? M.rows() : M.cols();
	lapack_int value_count;

	// Additional eigenvector dimensions:
	lapack_int* vectors_suppliements = new lapack_int[2 * lead_dim_vectors];
	
	lapack_int info = LAPACKE_zheevr(
		major_dim_order,				// param. 0
		job_token,						// param. 1
		range_token,					// param. 2
		upper_lower_token,				// param. 3
		lead_dim_vectors,				// param. 4
		M.data(),						// param. 5
		major_dim_length,				// param. 6
		r.lowest_value,					// param. 7
		r.highest_value,				// param. 8
		r.begin_index + 1,				// param. 9
		r.end_index + 1,				// param. 10
		0.,								// param. 11
		&value_count,					// param. 12
		values.data(),					// param. 13
		eigenvectors.data(),			// param. 14
		lead_dim_vectors,				// param. 15
		vectors_suppliements			// param. 16
		);

	delete[] vectors_suppliements;

	if (info != 0) {
		values.fill(nan(0));
	}

	if (value_count < lead_dim_vectors)
	{
		values.conservativeResize(value_count);
		eigenvectors.conservativeResize(lead_dim_vectors, value_count);
	}
	
	return result;
}