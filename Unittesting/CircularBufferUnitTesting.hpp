#ifndef CIRCULARBUFFER_UNITTESTING_H_
#define CIRCULARBUFFER_UNITTESTING_H_

#include <Misc/circularbuffer.hpp>
#include <Misc/parallelfor.hpp>
#include <functional> // for std::function
#include <iostream> // for std::cout
#include <sstream> // for std::sstream

namespace QuantumMechanics {

namespace MiscCircularBuffer {

namespace Unittesting {

	void test_sequential(std::function<void(std::string, bool)> assert_function) {

		CircularBuffer<100000> buffer;

		std::cout << "Buffer created!" << std::endl;

		size_t value = buffer.write("test", 5);

		assert_function("The MiscCircularBuffer could not retreive the correct written byte count.", value == 5);

		value = buffer.write("test2", 6);

		assert_function("The MiscCircularBuffer could not retreive the correct written byte count.", value == 6);

		char *data = new char[11];

		value = buffer.read(data, 11);

		assert_function("The MiscCircularBuffer could not retreive the correct read byte count.", value == 11);

		std::cout << "Data: " << data << data + 5 << std::endl;

		delete[] data;
	}

	void test_parallel(std::function<void(std::string, bool)> assert_function) {

		CircularBuffer<100000> buffer;

		std::cout << "Buffer created!" << std::endl;

		size_t value = buffer.write("test", 5);

		assert_function("The MiscCircularBuffer could not retreive the correct written byte count.", value == 5);

		parallelfor(long i = 0; i < 30; i++)
		{
			std::stringstream s;
			s << "test " << i / 2;

			if (!(i % 2))
			{
				value = buffer.write(s.str().c_str(), s.str().length() + 1);
				assert_function("The MiscCircularBuffer could not retreive the correct written byte count.", value == s.str().length() + 1);
			}
			else
			{
				buffer.wait_for_read();
				size_t valid_size = std::min(s.str().length() + 1, buffer.size());
				char *data = new char[valid_size];
				value = buffer.read(data, valid_size);
				assert_function("The MiscCircularBuffer could not retreive the correct read byte count.", value == s.str().length() + 1);
				std::cout << "Data: " << data << std::endl;
				delete[] data;
			}
		}
	}
		
	void test_all(std::function<void(std::string, bool)> assert_function) {

		std::cout << "CircularBuffer unittesting: test_sequential() ?" << std::endl;
		test_sequential(assert_function);
		std::cout << "Done! [CircularBuffer unittesting: test_sequential()]" << std::endl;

		std::cout << std::endl;

		std::cout << "CircularBuffer unittesting: test_parallel() ?" << std::endl;
		test_parallel(assert_function);
		std::cout << "Done! [CircularBuffer unittesting: test_parallel()]" << std::endl;

		std::cout << std::endl;

	}


} /* namespace UnitTesting */


} /* namespace MiscCircularBuffer */


} /* namespace QuantumMechanics */

#endif /* CIRCULARBUFFER_UNITTESTING_H_ */