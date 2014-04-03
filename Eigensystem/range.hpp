/*
Header file for Eigensystem:range

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _EIGENSYSTEM_RANGE_H_
# define _EIGENSYSTEM_RANGE_H_

namespace QuantumMechanics{
  
namespace Eigensystem{

struct range {
	enum range_type {
		full_range, index_range, mid_index_range, value_range
	} type_value;

	long begin_index, end_index;
	double lowest_value, highest_value;

	range() :
			type_value(full_range), begin_index(0), end_index(0), lowest_value(0.), highest_value(
					0.) {
	}
	;

	range(range_type type) :
			type_value(full_range), begin_index(0), end_index(0), lowest_value(0.), highest_value(
					0.) {
	}
	;

	range(range_type type, long begin, long end) :
			type_value(type), begin_index(
					(type == index_range || type == mid_index_range) ?
							begin : 0), end_index(
					(type == index_range || type == mid_index_range) ? end : 0), lowest_value(
					(type == value_range) ? begin : 0), highest_value(
					(type == value_range) ? end : 0) {
	}
	;

	range(range_type type, double lowest, double highest) :
			type_value(type), begin_index(
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
		return range(index_range, -count, -1);
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
		if (type_value == full_range || type_value == value_range)
			return;

		if (type_value == mid_index_range) {
			begin_index += size / 2;
			end_index += size / 2;
			type_value = index_range;
		}

		while (type_value == index_range && begin_index < 0)
			begin_index += size;

		while (type_value == index_range && end_index < 0)
			end_index += size;
	}

	inline bool operator==(const range &other) const {
		if(type_value != other.type_value)
			return false;

		if(type_value == full_range)
			return true;
		else if(type_value == index_range || type_value == mid_index_range)
		{
			if(begin_index == other.begin_index && end_index == other.end_index)
				return true;
			else
				return false;
		}
		else if(type_value == value_range)
		{
			if(lowest_value == other.lowest_value && highest_value == other.highest_value)
				return true;
			else
				return false;
		}

		return false;
	}


	inline bool operator!=(const range &other) const {
		return !( (*this) == other );
	}
};

}

}

#endif                                                      // _EIGENSYSTEM_RANGE_H_