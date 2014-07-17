#include "CircularBufferUnitTesting.hpp"
#include <string> // for std::string

using namespace QuantumMechanics::MiscCircularBuffer;

int main()
{
	std::function<void(std::string,bool)> assert_function = [&](std::string msg, bool assessment)
	{
		if(assessment == true)
			return;
		
		std::cout << "Assessment failed! Message:" << msg << std::endl;
	};
	
    std::cout << "Starting Eigensystem Unittesting:" << std::endl << std::endl;
	Unittesting::test_all(assert_function);
    std::cout  << std::endl << "Done with Eigensystem Unittesting!" << std::endl;
    return 0;
}