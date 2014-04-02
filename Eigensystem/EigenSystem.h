/*
 * EigenSystem.h
 *
 *  Created on: 05/02/2014
 *      Author: gregersen
 */

#ifndef EIGENSYSTEM_H_
#define EIGENSYSTEM_H_

#include <iostream>
#include <functional>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

namespace QuantumMechanics {

namespace EigenSystem {

struct range {
	enum range_type {
		full_range, index_range, mid_index_range, value_range
	} value;

	long begin_index, end_index;
	double lowest_value, highest_value;

	range() :
			value(full_range), begin_index(0), end_index(0), lowest_value(0.), highest_value(
					0.) {
	}
	;

	range(range_type type) :
			value(full_range), begin_index(0), end_index(0), lowest_value(0.), highest_value(
					0.) {
	}
	;

	range(range_type type, long begin, long end) :
			value(type), begin_index(
					(type == index_range || type == mid_index_range) ?
							begin : 0), end_index(
					(type == index_range || type == mid_index_range) ? end : 0), lowest_value(
					(type == value_range) ? begin : 0), highest_value(
					(type == value_range) ? end : 0) {
	}
	;

	range(range_type type, double lowest, double highest) :
			value(type), begin_index(
					(type == index_range || type == mid_index_range) ?
							lowest : 0), end_index(
					(type == index_range || type == mid_index_range) ?
							highest : 0), lowest_value(
					(type == value_range) ? lowest : 0), highest_value(
					(type == value_range) ? highest : 0) {
	}
	;

	static range full() {
		return range();
	}

	static range span(long begin, long end) {
		return range(index_range, begin, end);
	}

	static range lowest(long count) {
		return range(index_range, 0, count - 1);
	}

	static range highest(long count) {
		return range(index_range, -1, -count);
	}

	static range middle(long count) {
		return range(mid_index_range, -(count - 1) / 2, count / 2);
	}

	static range middle_span(long begin, long end) {
		return range(mid_index_range, begin, end);
	}

	static range values(double lowest, double highest) {
		return range(value_range, lowest, highest);
	}

	inline void fit_indices_to_size(const size_t &size) {
		if (value == full_range || value == value_range)
			return;

		if (value == mid_index_range) {
			begin_index += size / 2;
			end_index += size / 2;
			value == index_range;
		}

		while (value == index_range && begin_index < 0)
			begin_index += size;

		while (value == index_range && end_index < 0)
			end_index += size;
	}

	inline bool operator==(const range &other) {
		if(value != other.value)
			return false;

		if(value == full_range)
			return true;
		else if(value == index_range || value == mid_index_range)
		{
			if(begin_index == other.begin_index && end_index == other.end_index)
				return false;
			else
				return true;
		}
		else if(value == value_range)
		{
			if(lowest_value == other.lowest_value && highest_value == other.highest_value)
				return false;
			else
				return true;
		}

		return false;
	}


	inline bool operator!=(const range &other) {
		return !( (*this) == other );
	}
};

class HermitianSolver {
public:
	HermitianSolver();
	HermitianSolver(const MatrixXcd &M);
	HermitianSolver(size_t n, MatrixXcd * const M);
	HermitianSolver(size_t n, MatrixXcd * const M, const size_t &size);
	HermitianSolver(size_t n, const std::function<MatrixXcd(int)> &M, const size_t &size);
	virtual ~HermitianSolver();

	enum compute_action {
		EigenvaluesOnly,
		EigenvaluesAndVectors
	};

	void compute(compute_action action, range compute_range = range::full());

	MatrixXd eigenvalues();
	MatrixXd eigenvalues(range compute_range);

	static VectorXd eigenvalues(const MatrixXcd &M, range compute_range = range::full());

	std::vector<MatrixXcd> eigenvectors();
	std::vector<MatrixXcd> eigenvectors(range compute_range);

protected:
	void compute_eigenvalues();
	void compute_eigenvectors();

private:
	const size_t matrices_count;
	const size_t matrices_size;
	const MatrixXcd * const matrices_array;
	const std::function<MatrixXcd(int)> matrices_function;

	range computed_range;
	MatrixXd computed_eigenvalues;
	std::vector<MatrixXcd> computed_eigenvectors;
};

} /* namespace EigenSystem */

} /* namespace QuantumMechanics */

#endif /* EIGENSYSTEM_H_ */
