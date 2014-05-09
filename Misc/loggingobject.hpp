/*
Header file for QuantumMechanics::LoggingObject:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/
#ifndef _LOGGINGOBJECT_H_
#define _LOGGINGOBJECT_H_

#include <iostream>
#include <string>

namespace QuantumMechanics {

class LoggingObject {
	bool logging_enabled;

	static std::ostream null_stream;
	std::string objectIdenifier;

public:
	LoggingObject(const std::string &identifier, const bool &enabled = false) :
		logging_enabled(enabled),
		objectIdenifier(identifier)
	{ }

	void enable() {
		logging_enabled = true;
	}

	void disable() {
		logging_enabled = false;
	}

	virtual ~LoggingObject() { }

protected:
	std::ostream & log()
	{
		if (logging_enabled)
			return ( std::clog << ((objectIdenifier.empty()) ? "Message:" : objectIdenifier << " message: ") );
		else
			return null_stream;
	}

	std::ostream & logAppend()
	{
		if (logging_enabled)
			return std::clog;
		else
			return null_stream;
	}

public:
	std::ostream &operator()() {
		return log();
	}

	std::ostream &append() {
		return logAppend();
	}
};

std::ostream LoggingObject::null_stream(0);

};

#endif //namespace _LOGGINGOBJECT_H_