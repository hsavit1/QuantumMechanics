
#ifndef GREENSFORMALISM_UNITTESTING_H_
#define GREENSFORMALISM_UNITTESTING_H_

#include <Eigensystem/HermitianSolver>

using namespace Eigen;

namespace QuantumMechanics {

namespace GreensFormalism {

namespace Unittesting {

MatrixXcd random_hermitian(long n) {
	MatrixXcd result = MatrixXcd::Random(n,n);
	return (result += result.adjoint().eval());
}

void test_full_greens_inversion(std::function<void(std::string,bool)> assert_function) {

	MatrixXcd M = random_hermitian(10);

	GreensSolver solver(M);
	solver.compute(FullMatrix);
	
	assert_function("The GreensFormalism::GreensSolver could not solve a random hermitian 10x10 matrix.", solver.solution().isApprox(M.inverse()));
}

void test_all(std::function<void(std::string,bool)> assert_function) {

	std::cout << "GreensFormalism unittesting: test_full_greens_inversion() ?" << std::endl;
	test_full_greens_inversion(assert_function);
	std::cout << "Done! [GreensFormalism unittesting: test_full_greens_inversion()]" << std::endl;
	
	std::cout << std::endl;
}

} /* namespace UnitTesting */

} /* namespace GreensFormalism */

} /* namespace QuantumMechanics */

#endif /* GREENSFORMALISM_UNITTESTING_H_ */
