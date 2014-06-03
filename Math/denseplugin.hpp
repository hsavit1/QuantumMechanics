/*
---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/


template<typename LocalReturnType = EvalReturnType>
EIGEN_STRONG_INLINE LocalReturnType asZero() const
{
	return LocalReturnType::Zero(rows(), cols());
}

template<typename LocalReturnType = EvalReturnType>
EIGEN_STRONG_INLINE LocalReturnType asIdentity() const
{
	return LocalReturnType::Identity(rows(), cols());
}