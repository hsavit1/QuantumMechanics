#ifndef EIGENSYSTEM_UNITTESTING_H_
#define EIGENSYSTEM_UNITTESTING_H_

#include <Math/Dense>

using namespace Eigen;

namespace QuantumMechanics {

namespace Eigensystem {

namespace Unittesting {

	MatrixXcd random_hermitian(long n) {
		MatrixXcd result = MatrixXcd::Random(n, n);
		return (result += result.adjoint().eval());

	}

	void test_full_range(std::function<void(std::string, bool)> assert_function) {

		MatrixXcd M = random_hermitian(10);

		auto result = M.hermitianEigenvectors(range::full());

		MatrixXd D = result.first;
		MatrixXcd V = result.second;

		std::cout << D << std::endl << std::endl << V;

		assert_function("The return format of Eigensystem::HermitianSolver is incorrect. The number of eigenvalues is not consistent with matrix size.", D.rows() == 10);
		assert_function("The return format of Eigensystem::HermitianSolver is incorrect. The number of (or size of all) eigenvectors is not consistent with matrix size.", V.rows() == 10 && V.cols() == 10);

		assert_function("The Eigensystem::HermitianSolver could not solve a random hermitian 10x10 matrix.", M.isApprox(V * D * V.inverse()));

	}

	void test_ranges(std::function<void(std::string, bool)> assert_function) {

		MatrixXcd M = random_hermitian(10);

		auto result = M.hermitianEigenvectors(range::full());

		MatrixXd D = result.first;
		MatrixXcd V = result.second;

		assert_function("The Eigensystem::HermitianSolver could not solve a random hermitian 10x10 matrix.", (V * D).isApprox(M * V));

		// # range test 1!

		result = M.hermitianEigenvectors(range::lowest(4));

		MatrixXd D1 = result.first;
		MatrixXcd V1 = result.second;

		assert_function("The Eigensystem::HermitianSolver could not find the lowest 4 values.", D1.isApprox(D.topLeftCorner(4, 4)));
		assert_function("The Eigensystem::HermitianSolver could not find the lowest 4 vectors.", (V1 * D1).isApprox(M * V1));

		// # range test 2!

		result = M.hermitianEigenvectors(range::highest(4));

		MatrixXd D2 = result.first;
		MatrixXcd V2 = result.second;

		assert_function("The Eigensystem::HermitianSolver could not find the highest 4 values.", D2.isApprox(D.bottomRightCorner(4, 4)));
		assert_function("The Eigensystem::HermitianSolver could not find the highest 4 vectors.", (V2 * D2).isApprox(M * V2));

		// # range test 3!

		result = M.hermitianEigenvectors(range::span(3, 8));

		MatrixXd D3 = result.first;
		MatrixXcd V3 = result.second;

		assert_function("The Eigensystem::HermitianSolver could not find the span(3,8) values.", D3.isApprox(D.block(3, 3, 6, 6)));
		assert_function("The Eigensystem::HermitianSolver could not find the span(3,8) values.", (V3 * D3).isApprox(M * V3));

		// # range test 4!

		result = M.hermitianEigenvectors(range::middle(4));

		MatrixXd D4 = result.first;
		MatrixXcd V4 = result.second;
		
		assert_function("The Eigensystem::HermitianSolver could not find the middle 4 values.", D4.isApprox(D.block(4, 4, 4, 4)));
		assert_function("The Eigensystem::HermitianSolver could not find the middle 4 vectors.", (V4 * D4).isApprox(M * V4));

		// # range test 5!

		double low = (D(3, 3) + D(4, 4)) / 2.;
		double high = (D(7, 7) + D(8, 8)) / 2.;

		result = M.hermitianEigenvectors(range::values(low, high));

		MatrixXd D5 = result.first;
		MatrixXcd V5 = result.second;

		assert_function("The Eigensystem::HermitianSolver could not find the 4 values between a low/high range.", D5.isApprox(D.block(4, 4, 4, 4)));
		assert_function("The Eigensystem::HermitianSolver could not find the 4 vectors between a low/high range.", (V5 * D5).isApprox(M * V5));

	}
	
	void test_all(std::function<void(std::string, bool)> assert_function) {

		std::cout << "Eigensystem unittesting: test_full_range() ?" << std::endl;
		test_full_range(assert_function);
		std::cout << "Done! [Eigensystem unittesting: test_full_range()]" << std::endl;

		std::cout << std::endl;

		std::cout << "Eigensystem unittesting: test_ranges() ?" << std::endl;
		test_ranges(assert_function);
		std::cout << "Done! [Eigensystem unittesting: test_ranges()]" << std::endl;

		std::cout << std::endl;

	}


} /* namespace UnitTesting */


} /* namespace EigenSystem */


} /* namespace QuantumMechanics */

#endif /* EIGENSYSTEM_UNITTESTING_H_ */