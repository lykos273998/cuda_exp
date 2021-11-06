#include <iostream>
#include <memory>
#include <cmath>
#include <vector>
#include "../../eigen/Eigen/Dense"
#include "../../eigen/Eigen/Sparse"

template<typename T>
void DiscreteLaplacian1D(size_t n,Eigen::SparseMatrix<T> &mat){
	std::vector<Eigen::Triplet<T>> triplets(4+(n-2)*3);
	int k=0;
	triplets[0] = Eigen::Triplet<T>(0,0,-1.);
	++k;
	triplets[k] = Eigen::Triplet<T>(0,1,2);
	++k;
	for(int i = 1; i < n-1; ++i){
		triplets[k] = Eigen::Triplet<T>(i,i-1,-1);
		++k;
		triplets[k] = Eigen::Triplet<T>(i,i,2);
		++k;
		triplets[k] = Eigen::Triplet<T>(i,i+1,-1);
		++k;
	}	
	
	
	triplets[k] = Eigen::Triplet<T>(n-1,n-2,-1);
	++k;	
	triplets[k] = Eigen::Triplet<T>(n-1,n-1,2);
	mat.setFromTriplets(triplets.begin(),triplets.end());
}

template<typename T>
void DiscreteLaplacian2D(size_t n,Eigen::SparseMatrix<T> &mat){
	std::vector<Eigen::Triplet<T>> triplets;
	for(int i = 0; i < n*n; ++i){
		int j = i;
		triplets.push_back(Eigen::Triplet<T>(i,j,4));
		
		
		j = i+1;
		if(j <= n*n - 1) triplets.push_back(Eigen::Triplet<T>(i,j,-1.));
		
			
		j = i+n;
		if(j <= n*n - 1) triplets.push_back(Eigen::Triplet<T>(i,j,-1.));
		
		
		j = i-1;
		if(j >= 0) triplets.push_back(Eigen::Triplet<T>(i,j,-1.));
		
			
		j = i-n;
		if(j >= 0) triplets.push_back(Eigen::Triplet<T>(i,j,-1.));
	}

	mat.setFromTriplets(triplets.begin(), triplets.end());
}


int main()
{
	int n = 5;
	Eigen::SparseMatrix<double> L(n*n,n*n);
	DiscreteLaplacian2D(n,L);
}
