
#ifndef EIGENSYSTEM_UnitTesting_H_
#define EIGENSYSTEM_UnitTesting_H_

#include "../EigenSystem/EigenSystem.h"

using namespace Eigen;

namespace QuantumMechanics {

namespace EigenSystem {

namespace UnitTesting {

MatrixXcd random_hermitian(long n) {
	MatrixXcd result = MatrixXcd::Random(n,n);
	return (result += result.adjoint().eval());
}

void test_full_range(std::function<void(std::string,bool)> assert_function) {

	MatrixXcd M = random_hermitian(10);

	HermitianSolver solver(M);
	solver.compute(HermitianSolver::EigenvaluesAndVectors);

	MatrixXd D = solver.eigenvalues().asDiagonal();
	MatrixXcd V = solver.eigenvectors();

	assert_function("The EigenSystem::HermitianSolver could not solve a random hermitian 10x10 matrix.", M.isApprox(V * D * V.inverse() ));
}

void test_ranges(std::function<void(std::string,bool)> assert_function) {

	MatrixXcd M = random_hermitian(10);

	HermitianSolver solver(M);
	solver.compute(HermitianSolver::EigenvaluesAndVectors);

	MatrixXd D = solver.eigenvalues().asDiagonal();
	MatrixXcd V = solver.eigenvectors();

	assert_function("The EigenSystem::HermitianSolver could not solve a random hermitian 10x10 matrix.", M.isApprox(V * D * V.inverse() ));

	// # range test 1!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::lowest(4));

	MatrixXd D1 = solver.eigenvalues().asDiagonal();
	MatrixXcd V1 = solver.eigenvectors();

	assert_function("The EigenSystem::HermitianSolver could not find the lowest 4 values.", D1.isApprox(D.topLeftCorner(4,4)));
	assert_function("The EigenSystem::HermitianSolver could not find the lowest 4 vectors.", V1.isApprox(V.leftCols(4)));

	// # range test 2!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::highest(4));

	MatrixXd D2 = solver.eigenvalues().asDiagonal();
	MatrixXcd V2 = solver.eigenvectors();

	assert_function("The EigenSystem::HermitianSolver could not find the highest 4 values.", D2.isApprox(D.bottomRightCorner(4,4)));
	assert_function("The EigenSystem::HermitianSolver could not find the highest 4 vectors.", V2.isApprox(V.rightCols(4)));

	// # range test 3!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::span(3,8));

	MatrixXd D3 = solver.eigenvalues().asDiagonal();
	MatrixXcd V3 = solver.eigenvectors();

	assert_function("The EigenSystem::HermitianSolver could not find the highest 4 values.", D3.isApprox(D.block(3,3,6,6)));
	assert_function("The EigenSystem::HermitianSolver could not find the highest 4 vectors.", V3.isApprox(V.cols(3,8)));

	// # range test 4!

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::middle(4));

	MatrixXd D4 = solver.eigenvalues().asDiagonal();
	MatrixXcd V4 = solver.eigenvectors();

	assert_function("The EigenSystem::HermitianSolver could not find the highest 4 values.", D4.isApprox(D.block(4,4,4,4)));
	assert_function("The EigenSystem::HermitianSolver could not find the highest 4 vectors.", V4.isApprox(V.cols(4,7)));

	// # range test 5!

	double low = (solver.eigenvalues()(3,0) + solver.eigenvalues()(4,0)) / 2.;
	double high = (solver.eigenvalues()(7,0) + solver.eigenvalues()(8,0)) / 2.;

	solver.compute(HermitianSolver::EigenvaluesAndVectors, range::values(low,high));

	MatrixXd D5 = solver.eigenvalues().asDiagonal();
	MatrixXcd V5 = solver.eigenvectors();

	assert_function("The EigenSystem::HermitianSolver could not find the highest 4 values.", D4.isApprox(D.block(4,4,4,4)));
	assert_function("The EigenSystem::HermitianSolver could not find the highest 4 vectors.", V4.isApprox(V.cols(4,7)));
}

} /* namespace UnitTesting */

} /* namespace EigenSystem */

} /* namespace QuantumMechanics */

#endif /* EIGENSYSTEM_UnitTesting_H_ */
