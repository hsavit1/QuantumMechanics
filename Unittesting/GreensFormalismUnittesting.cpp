#include "GreensFormalismUnittesting.hpp"

using namespace QuantumMechanics::GreensFormalism;

int main()
{
	std::function<void(std::string,bool)> assert_function = [&](std::string msg, bool assessment)
	{
		if(assessment == true)
			return;
		
		std::cout << "Assessment failed! Message:" << msg << std::endl;
	};
	
    std::cout << "Starting GreensFormalism Unittesting:" << std::endl << std::endl;
	Unittesting::test_all(assert_function);
    std::cout  << std::endl << "Done with GreensFormalism Unittesting!" << std::endl;
    return 0;
}