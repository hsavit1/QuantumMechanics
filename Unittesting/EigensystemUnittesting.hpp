
#ifndef EIGENSYSTEM_UNITTESTING_H_
#define EIGENSYSTEM_UNITTESTING_H_

#include <Eigensystem/HermitianSolver>

using namespace Eigen;

namespace QuantumMechanics {

namespace Eigensystem {

namespace Unittesting {

MatrixXcd random_hermitian(long n) {
	MatrixXcd result = MatrixXcd::Random(n,n);
	return (result += result.adjoint().eval());
}

void test_full_range(std::function<void(std::string,bool)> assert_function) {

	MatrixXcd M = random_hermitian(10);

	HermitianSolver solver(M);
	solver.compute(HermitianSolver::EigenvaluesAndVectors);

	MatrixXd D = solver.eigenvalues().asDiagonal();
	MatrixXcd V = solver.eigenvectors().front();
	
	assert_function("The return format of Eigensystem::HermitianSolver is incorrect. The number of eigenvalues is not consistent with matrix size.", D.rows() == 10 && D.cols() == 10);
	assert_function("The return format of Eigensystem::HermitianSolver is incorrect. The number of (or size of all) eigenvectors is not consistent with matrix size.", V.rows() == 10 && V.cols() == 10);

	assert_function("The Eigensystem::HermitianSolver could not solve a random hermitian 10x10 matrix.", M.isApprox(V * D * V.inverse() ));
}

void test_ranges(std::function<void(std::string,bool)> assert_function) {

	MatrixXcd M = random_hermitian(10);

	HermitianSolver solver(M);
	solver.compute(HermitianSolver::EigenvaluesAndVectors);

	MatrixXd D = solver.eigenvalues().asDiagonal();
	MatrixXcd V = solver.eigenvectors().front();

	assert_function("The Eigensystem::HermitianSolver could not solve a random hermitian 10x10 matrix.", (V * D).isApprox(M * V) );

	// # range test 1!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::lowest(4));

	MatrixXd D1 = solver.eigenvalues().asDiagonal();
	MatrixXcd V1 = solver.eigenvectors().front();

	assert_function("The Eigensystem::HermitianSolver could not find the lowest 4 values.", D1.isApprox(D.topLeftCorner(4,4)));
	assert_function("The Eigensystem::HermitianSolver could not find the lowest 4 vectors.", (V1 * D1).isApprox(M * V1));

	// # range test 2!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::highest(4));

	MatrixXd D2 = solver.eigenvalues().asDiagonal();
	MatrixXcd V2 = solver.eigenvectors().front();
	
	assert_function("The Eigensystem::HermitianSolver could not find the highest 4 values.", D2.isApprox(D.bottomRightCorner(4,4)));
	assert_function("The Eigensystem::HermitianSolver could not find the highest 4 vectors.", (V2 * D2).isApprox(M * V2));

	// # range test 3!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::span(3,8));

	MatrixXd D3 = solver.eigenvalues().asDiagonal();
	MatrixXcd V3 = solver.eigenvectors().front();

	assert_function("The Eigensystem::HermitianSolver could not find the span(3,8) values.", D3.isApprox(D.block(3,3,6,6)));
	assert_function("The Eigensystem::HermitianSolver could not find the span(3,8) values.", (V3 * D3).isApprox(M * V3));

	// # range test 4!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::middle(4));

	MatrixXd D4 = solver.eigenvalues().asDiagonal();
	MatrixXcd V4 = solver.eigenvectors().front();

	assert_function("The Eigensystem::HermitianSolver could not find the middle 4 values.", D4.isApprox(D.block(4,4,4,4)));
	assert_function("The Eigensystem::HermitianSolver could not find the middle 4 vectors.", (V4 * D4).isApprox(M * V4));

	// # range test 5!

	double low = (D(3,3) + D(4,4)) / 2.;
	double high = (D(7,7) + D(8,8)) / 2.;

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::values(low,high));

	MatrixXd D5 = solver.eigenvalues().asDiagonal();
	MatrixXcd V5 = solver.eigenvectors().front();
	
	assert_function("The Eigensystem::HermitianSolver could not find the 4 values between a low/high range.", D5.isApprox(D.block(4,4,4,4)));
	assert_function("The Eigensystem::HermitianSolver could not find the 4 vectors between a low/high range.", (V5 * D5).isApprox(M * V5));
}

void test_ranges_multiple_matrices(std::function<void(std::string,bool)> assert_function) {
	
	MatrixXcd *M = new MatrixXcd[4];
	
	M[0] = random_hermitian(10);
	M[1] = random_hermitian(10);
	M[2] = random_hermitian(10);
	M[3] = random_hermitian(10);
	
	HermitianSolver solver(4,M);
	
	std::function<bool(MatrixXcd,VectorXd,MatrixXcd)> test_single_result = [](MatrixXcd m, VectorXd d, MatrixXcd v)
	{
		return (v * d.head(v.cols()).asDiagonal()).isApprox(m * v);
	};
	
	std::function<bool()> test_all_results = [&]()
	{
		bool result = true;
		for(int i = 0; i < 4 && result == true; i++)
		{
			result &= test_single_result(M[i],solver.eigenvalues().col(i),solver.eigenvectors()[i]);
		}
		
		return result;
	};

	solver.compute(HermitianSolver::EigenvaluesAndVectors);

	assert_function("The Eigensystem::HermitianSolver could not solve a random hermitian 10x10 matrix.", test_all_results() );

	// # range test 1!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::lowest(4));

	assert_function("The Eigensystem::HermitianSolver could not find the lowest 4 values and vectors.", test_all_results());

	// # range test 2!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::highest(4));

	assert_function("The Eigensystem::HermitianSolver could not find the highest 4 values and vectors.", test_all_results());

	// # range test 3!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::span(3,8));

	assert_function("The Eigensystem::HermitianSolver could not find the span(3,8) values and vectors.", test_all_results());

	// # range test 4!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::middle(4));

	assert_function("The Eigensystem::HermitianSolver could not find the middle 4 values and vectors.", test_all_results());

	// # range test 5!

	double low = -1;
	double high = 1;

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::values(low,high));
	
	assert_function("The Eigensystem::HermitianSolver could not find the values and vectors between a low/high range.", test_all_results());
}

void test_ranges_multiple_matrices_varying_sizes(std::function<void(std::string,bool)> assert_function) {
	
	MatrixXcd *M = new MatrixXcd[4];
	
	M[0] = random_hermitian(10);
	M[1] = random_hermitian(11);
	M[2] = random_hermitian(12);
	M[3] = random_hermitian(10);
	
	HermitianSolver solver(4,M,9);
	
	std::function<bool(MatrixXcd,VectorXd,MatrixXcd)> test_single_result = [&](MatrixXcd m, VectorXd d, MatrixXcd v)
	{
		return (v * d.head(v.cols()).asDiagonal()).isApprox(m * v);
	};
	
	std::function<bool()> test_all_results = [&]()
	{
		bool result = true;
		for(int i = 0; i < 4 && result == true; i++)
		{
			result &= test_single_result(M[i].block(0,0,9,9),solver.eigenvalues().col(i),solver.eigenvectors()[i]);
		}
		
		return result;
	};

	solver.compute(HermitianSolver::EigenvaluesAndVectors);

	assert_function("The Eigensystem::HermitianSolver could not solve a random hermitian 10x10 matrix.", test_all_results() );

	// # range test 1!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::lowest(4));

	assert_function("The Eigensystem::HermitianSolver could not find the lowest 4 values and vectors.", test_all_results());

	// # range test 2!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::highest(4));

	assert_function("The Eigensystem::HermitianSolver could not find the highest 4 values and vectors.", test_all_results());

	// # range test 3!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::span(3,8));

	assert_function("The Eigensystem::HermitianSolver could not find the span(3,8) values and vectors.", test_all_results());

	// # range test 4!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::middle(4));

	assert_function("The Eigensystem::HermitianSolver could not find the middle 4 values and vectors.", test_all_results());

	// # range test 5!

	double low = -1;
	double high = 1;

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::values(low,high));
	
	assert_function("The Eigensystem::HermitianSolver could not find the values and vectors between a low/high range.", test_all_results());
}

void test_all(std::function<void(std::string,bool)> assert_function) {

	std::cout << "Eigensystem unittesting: test_full_range() ?" << std::endl;
	test_full_range(assert_function);
	std::cout << "Done! [Eigensystem unittesting: test_full_range()]" << std::endl;
	
	std::cout << std::endl;
	
	std::cout << "Eigensystem unittesting: test_ranges() ?" << std::endl;
	test_ranges(assert_function);
	std::cout << "Done! [Eigensystem unittesting: test_ranges()]" << std::endl;
	
	std::cout << std::endl;
	
	std::cout << "Eigensystem unittesting: test_ranges_multiple_matrices() ?" << std::endl;
	test_ranges_multiple_matrices(assert_function);
	std::cout << "Done! [Eigensystem unittesting: test_ranges_multiple_matrices()]" << std::endl;
	
	std::cout << std::endl;
	
	std::cout << "Eigensystem unittesting: test_ranges_multiple_matrices_varying_sizes() ?" << std::endl;
	test_ranges_multiple_matrices_varying_sizes(assert_function);
	std::cout << "Done! [Eigensystem unittesting: test_ranges_multiple_matrices_varying_sizes()]" << std::endl;
}

} /* namespace UnitTesting */

} /* namespace EigenSystem */

} /* namespace QuantumMechanics */

#endif /* EIGENSYSTEM_UNITTESTING_H_ */
