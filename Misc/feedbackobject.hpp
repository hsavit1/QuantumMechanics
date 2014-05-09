/*
Header file for QuantumMechanics::FeedbackObject:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/
#ifndef _FEEDBACKOBJECT_H_
#define _FEEDBACKOBJECT_H_

#include <functional>
#include <tbb/tbb.h>

namespace QuantumMechanics {

class FeedbackObject {

	std::function<void(double)> feedback_function;

	typedef tbb::atomic<double> progress_counter;
	tbb::enumerable_thread_specific<progress_counter> local_counters;
	// zero_allocator is essential here.
	tbb::concurrent_vector<progress_counter*, tbb::zero_allocator<progress_counter*> > local_counter_pointers;

public:
	FeedbackObject() : feedback_function(nullptr) { } 
	virtual ~FeedbackObject() { }

	void enableFeedback(std::function<void(double)> function) {
		feedback_function = function;
	}

protected:
	void addToProgress(double delta) {
		bool exists;
		auto& i = local_counters.local(exists);
		i += delta;
		if (!exists)
			// First time we've seen this local counter.
			local_counter_pointers.push_back(&i);
	}

	double getProgress() {
		double sum = 0;
		size_t n = local_counter_pointers.size();
		for (size_t i = 0; i<n; ++i)
			// "if" deals with timing hold where slot in LocalCounterPointers was allocated but not initialized.
			if (auto* j = local_counter_pointers[i])
				sum += *j;
		return sum;
	}

	// Can be called asynchronously.
	void clearProgress() {
		size_t n = local_counter_pointers.size();
		for (size_t i = 0; i<n; ++i)
			// "if" deals with timing hold where slot in LocalCounterPointers was allocated but not initialized.
			if (auto* j = local_counter_pointers[i])
				*j = 0;
	}

public:
	void updateFeedback(double delta) {
		if (feedback_function)
		{
			addToProgress(delta);
			feedback_function(getProgress());
		}
	}

	void resetFeedback() {
		clearProgress();
	}
};

};

#endif //namespace _FEEDBACKOBJECT_H_