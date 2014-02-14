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

	void fit_indices_to_order(const size_t &order) {
		if (value == full_range || value == value_range)
			return;

		if (value == mid_index_range) {
			begin_index += order / 2;
			end_index += order / 2;
			value == index_range;
		}

		while (value == index_range && begin_index < 0)
			begin_index += order;

		while (value == index_range && end_index < 0)
			end_index += order;
	}
};

}
;
/* namespace EigenSystem */

class HermitianEigenSolver {
public:
	HermitianEigenSolver();
	HermitianEigenSolver(size_t n, MatrixXcd * const M);
	HermitianEigenSolver(size_t n, MatrixXcd * const M, const size_t &order);
	HermitianEigenSolver(size_t n, const std::function<MatrixXcd(int)> &M,
			const size_t &order);
	HermitianEigenSolver(const std::vector<MatrixXcd> &M);
	virtual ~HermitianEigenSolver();

	ArrayXXd eigenvalues(EigenSystem::range range =
			EigenSystem::range::full()) const;

	static ArrayXd eigenvalues(const MatrixXcd &M, EigenSystem::range range =
			EigenSystem::range::full());

private:
	const size_t points;
	const size_t order;
	MatrixXcd * const M_array;
	const std::function<MatrixXcd(int)> M_function;
};

} /* namespace QuantumMechanics */

#endif /* EIGENSYSTEM_H_ */
