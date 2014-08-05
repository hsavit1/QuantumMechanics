/*
Header file for QuantumMechanics::GreensFormalism::GreensSolver: 

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the 

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
 */
#ifndef _GREENSFORMALISM_GREENSSOLVER_H_
#define _GREENSFORMALISM_GREENSSOLVER_H_

#include <Math/Dense>
#include "../misc/LoggingObject"

#include <vector>

namespace QuantumMechanics {

namespace GreensFormalism {

	enum GreenMatrixSubType {
		FullMatrix,
		FirstBlock,
		LastBlock,
		FirstBlockColumn,
		LastBlockColumn
	};
	
class GreensSolver {

	const BlockMatrixXcd &H;
	MatrixXcd sigma;
	BlockMatrixXcd G;

	static LoggingObject log;

public:
	GreensSolver(const BlockMatrixXcd &M) : H(M), sigma(), G() {}

	GreensSolver(const MatrixXcd &M) : H(M), sigma(), G() {}

	static inline void enableLog()
	{
		log.enable();
	}

protected:
	void compute_full_matrix() 
	{
		const long block_count = (H.isSquare() || H.blockRows() < H.blockCols() ? H.blockRows() : H.blockCols());

		log() << "Preparing to calculate the full solution of " << block_count << "-by-" << block_count << " blocks." << std::endl;

		sigma = H;

		log() << "The reduced sigma has been set to zeros." << std::endl;

		G = H.inverse();

		log() << "The solution is saved." << std::endl;
	}

	void compute_last_block()
	{
		const long block_count = (H.isSquare() || H.blockRows() < H.blockCols() ? H.blockRows() : H.blockCols());

		log() << "Preparing to calculate the last block out of " << block_count << "-by-" << block_count << " blocks." << std::endl;
		
		sigma = H.block(0, 0).asZero();

		log() << "The algorithm wil recursively find the self-energy of the left cells." << std::endl;

		for (long b = 0; b < block_count - 1; b++)
			sigma = H.block(b + 1, b) * (H.block(b, b) - sigma).inverse() * H.block(b, b + 1);

		log() << "The final self-energy became:" << std::endl << std::endl << sigma << std::endl << std::endl;

		if (block_count > 1)
			G = (H.block(-1, -1) - sigma).inverse();
		else
			G = sigma.inverse();

		log() << "The solution is saved." << std::endl;
	}
	
	void compute_first_block() 
	{
		const long block_count = (H.isSquare() || H.blockRows() < H.blockCols() ? H.blockRows() : H.blockCols());

		log() << "Preparing to calculate the last block out of " << block_count << "-by-" << block_count << " blocks." << std::endl;

		sigma = H.block(-1, -1).asZero();

		log() << "The algorithm wil recursively find the self-energy of the left cells." << std::endl;

		for (long b = -1; b >= -(block_count - 1); b--)
			sigma = H.block(b - 1, b) * (H.block(b, b) - sigma).inverse() * H.block(b, b - 1);

		log() << "The final self-energy became:" << std::endl << std::endl << sigma << std::endl << std::endl;

		if (block_count > 1)
			G = (H.block(0, 0) - sigma).inverse();
		else
			G = sigma.inverse();

		log() << "The solution is saved." << std::endl;
	}
	
	void compute_first_block_column() 
	{
		const long block_count = (H.isSquare() || H.blockRows() < H.blockCols() ? H.blockRows() : H.blockCols());

		log() << "Preparing to calculate the first block column out of " << block_count << "-by-" << block_count << " blocks." << std::endl;

		std::vector<MatrixXcd> g(block_count, MatrixXcd());

		sigma = H.block(-1, -1).asZero();

		log() << "The algorithm wil recursively find the self-energy of the right cells while saving intermediate isolated greens matrices." << std::endl;

		for (long b = -1; b > -block_count; b--)
		{
			g[-b - 1] = (H.block(b, b) - sigma).inverse();
			sigma = H.block(b - 1, b) * g[-b - 1] * H.block(b, b - 1);
		}

		log() << "The final self-energy became:" << std::endl << std::endl << sigma << std::endl << std::endl;

		log() << "The solution is is a column block " << H.blockRows() << "-by-1 matrix." << std::endl;

		G = H.blocks(0, 0, block_count, 1).asZero();

		log() << "The solution is calculated from the intermediate greens matrices." << std::endl;

		G.block(0, 0) = (H.block(0, 0) - sigma).inverse();

		log() << "Block 0 is calculated." << std::endl;

		for (long b = 1; b < block_count; b++)
		{
			log() << "Block "<< b << " is calculated." << std::endl;
			G.block(b, 0) = - g[block_count - 1 - b] * H.block(b, b - 1) * G.block(b - 1, 0);
		}

		log() << "The solution is finished." << std::endl;
	}

	void compute_last_block_column()
	{
		const long block_count = (H.isSquare() || H.blockRows() < H.blockCols() ? H.blockRows() : H.blockCols());

		log() << "Preparing to calculate the last block column out of " << block_count << "-by-" << block_count << " blocks." << std::endl;

		std::vector<MatrixXcd> g(block_count, MatrixXcd());

		sigma = H.block(0, 0).asZero();

		log() << "The algorithm wil recursively find the self-energy of the right cells while saving intermediate isolated greens matrices." << std::endl;

		for (long b = 0; b < block_count - 1; b++)
		{
			g[b] = (H.block(b, b) - sigma).inverse();
			sigma = H.block(b + 1, b) * g[b] * H.block(b, b + 1);
		}

		log() << "The final self-energy became:" << std::endl << std::endl << sigma << std::endl << std::endl;

		log() << "The solution is is a column block " << H.blockRows() << "-by-1 matrix." << std::endl;

		G = H.blocks(0, -1, block_count, 1).asZero();

		log() << "The solution is calculated from the intermediate greens matrices." << std::endl;

		G.block(-1, 0) = (H.block(-1, -1) - sigma).inverse();

		log() << "Block 0 is calculated." << std::endl;

		for (long b = -2; b >= -block_count; b--)
		{	
			log() << "Block "<< -b - 1 << " is calculated." << std::endl;
			G.block(b, 0) = - g[block_count + b] * H.block(b, b + 1) * G.block(b + 1, 0);
		}

		log() << "The solution is finished." << std::endl;
	}
	
public:
	void compute(const GreenMatrixSubType &action)
	{
		switch(action)
		{
		case FullMatrix:
			compute_full_matrix();
			break;
		case FirstBlock:
			compute_first_block();
			break;
		case LastBlock:
			compute_last_block();
			break;
		case FirstBlockColumn:
			compute_first_block_column();
			break;
		case LastBlockColumn:
			compute_last_block_column();
			break;
		}
	}

	const MatrixXcd &reducedSigma() {
		return sigma;
	}


	const BlockMatrixXcd &greensMatrix() const {
		return G;
	}
};

LoggingObject GreensSolver::log("GreensFormalism::GreensSolver", false);

}

}


#endif  
