
#ifndef GREENSFORMALISM_UNITTESTING_H_
#define GREENSFORMALISM_UNITTESTING_H_

#include <GreensFormalism/GreensSolver>
#include <GreensFormalism/ChainSolver>
#include <LanduarFormalism/TwoLeadTransportSolver>

namespace QuantumMechanics {

namespace GreensFormalism {

namespace Unittesting {

	const bool log = false;

	BlockMatrixXcd random_hermitian(ArrayXi sizes) {

	if (log) std::cout << "A random hermitian matrix is made!" << std::endl;

	const long full_size = sizes.sum();

	if (log) std::cout << "We initialize with zeros." << std::endl;

	BlockMatrixXcd result = MatrixXcd::Zero(full_size,full_size);

	if (log) std::cout << "We define blocks." << std::endl;

	result.setBlocks(sizes);

	if (log) std::cout << "We randomly set the diagonal as well as lower and upper diagonal blocks." << std::endl;

	for(int i = 0; i < sizes.size(); i++)
	{
		if (log) std::cout << "Diagonal block " << i << " is set." << std::endl;

		result.block(i, i) = MatrixXcd::Random(sizes[i],sizes[i]);

		if( i < sizes.size() - 1)
		{
			if (log) std::cout << "Upper diagonal block " << i << " is set." << std::endl;
			result.block(i + 1, i) = MatrixXcd::Random(sizes[i + 1], sizes[i]);
			if (log) std::cout << "Lower diagonal block " << i << " is set." << std::endl;
			result.block(i, i + 1) = MatrixXcd::Random(sizes[i], sizes[i + 1]);
		}
	}

	if (log) std::cout << "The matrix is added to its adjoint." << std::endl;
	
	result += result.adjoint().eval();

	if (log) std::cout << "The matrix is finished." << std::endl;

	return result;
}

void test_full_greens_inversion(std::function<void(std::string,bool)> assert_function) {

	ArrayXi sizes = Array4i(2,2,2,2);
	BlockMatrixXcd M = random_hermitian(sizes);

	GreensSolver solver(M);
	//solver.enableLog();

	solver.compute(FullMatrix);
	
	assert_function("The GreensFormalism::GreensSolver could not solve a random hermitian 10x10 matrix.", solver.greensMatrix().matrix().isApprox(M.matrix().inverse()));
}

void test_partial_greens_inversion(std::function<void(std::string,bool)> assert_function) {

	ArrayXi sizes = Array4i(2,3,2,3);
	BlockMatrixXcd M = random_hermitian(sizes);

	GreensSolver solver(M);
	//solver.enableLog();

	solver.compute(LastBlock);
	
	assert_function("The GreensFormalism::GreensSolver could not solve the last block of a random hermitian 10x10 matrix.", solver.greensMatrix().matrix().isApprox(M.matrix().inverse().block(7, 7, 3, 3)));

	solver.compute(FirstBlock);
	
	assert_function("The GreensFormalism::GreensSolver could not solve the first block of a random hermitian 10x10 matrix.", solver.greensMatrix().matrix().isApprox(M.matrix().inverse().block(0, 0, 2, 2)));

	solver.compute(LastBlockColumn);

	assert_function("The GreensFormalism::GreensSolver could not solve the last block column of a random hermitian 10x10 matrix.", solver.greensMatrix().matrix().isApprox(M.matrix().inverse().block(0, 7, 10, 3)));

	solver.compute(FirstBlockColumn);
	
	assert_function("The GreensFormalism::GreensSolver could not solve the first block column of a random hermitian 10x10 matrix.", solver.greensMatrix().matrix().isApprox(M.matrix().inverse().block(0, 0, 10, 2), 1e-11)); // default precision is 1e-12!
}

void test_chain_surface_greens(std::function<void(std::string, bool)> assert_function) {

	ArrayXi sizes = Array4i(2, 2, 2, 2);
	BlockMatrixXcd h = random_hermitian(sizes);
	BlockMatrixXcd v = random_hermitian(sizes);

	ChainSolver solver(h,v);
	//solver.enableLog();

	solver.compute(SurfaceGreensMatrix);

	assert_function("The GreensFormalism::ChainSolver could not solve a random hermitian 10x10 hamilton matrix and a 10x10 hopping matrix.", solver.greensMatrix().matrix().size() > 0);
}

void test_all(std::function<void(std::string,bool)> assert_function) {

	std::cout << "GreensFormalism unittesting: test_full_greens_inversion() ?" << std::endl;
	test_full_greens_inversion(assert_function);
	std::cout << "Done! [GreensFormalism unittesting: test_full_greens_inversion()]" << std::endl;

	std::cout << std::endl;

	std::cout << "GreensFormalism unittesting: test_partial_greens_inversion() ?" << std::endl;
	test_partial_greens_inversion(assert_function);
	std::cout << "Done! [GreensFormalism unittesting: test_partial_greens_inversion()]" << std::endl;

	std::cout << std::endl;

	std::cout << "GreensFormalism unittesting: test_chain_surface_greens() ?" << std::endl;
	test_chain_surface_greens(assert_function);
	std::cout << "Done! [GreensFormalism unittesting: test_chain_surface_greens()]" << std::endl;

	std::cout << std::endl;
}

} /* namespace UnitTesting */

} /* namespace GreensFormalism */

} /* namespace QuantumMechanics */

#endif /* GREENSFORMALISM_UNITTESTING_H_ */
