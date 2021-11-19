#include <stdio.h>
#include <stdlib.h>
 #define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


typedef struct{

	size_t nx;
	size_t ny;
	size_t nz;

	double ***data;
	double **second_level;
	double *last_level;
}Matrix3D_double;

void get_matrix_dimension(size_t * my_size, size_t* starting_position, int my_rank, int num_procs, size_t tot_size){
	(*my_size)=tot_size/num_procs;
	(*starting_position) = (*my_size)*my_rank + min(my_rank,tot_size % num_procs);
	if(my_rank < tot_size % num_procs){
		(*my_size)++;
	}


}

void allocate_matrix3D(Matrix3D_double *mat, size_t nx, size_t ny, size_t nz){
	(*mat).nx = nx;
	(*mat).ny = ny;
	(*mat).nz = nz;
	(*mat).last_level = (double*)malloc(nx*ny*nz*sizeof(double));
	(*mat).second_level = (double**)malloc(nx*ny*sizeof(double*));
	(*mat).data = (double***)malloc(nx*sizeof(double**));
	for(size_t i=0; i<nx; ++i){
		for(size_t j=0; j<ny; ++j){
			(*mat).second_level[i*ny+j] = &(*mat).last_level[(i*ny+j)*nz];
			
		}
		(*mat).data[i] = &(*mat).second_level[i*ny];
	}
}

void sum_matrix_3D(Matrix3D_double* A, Matrix3D_double* B, Matrix3D_double* C){
	if(A.nx != B.nx || A.ny != B.ny || A.nz != B.nz){

		printf("Mismatched sizes of matrices A and B in sum\n");
		exit(1)
		return }


}
void free_matrix3D(Matrix3D_double *mat){
	free( (*mat).last_level);
	free((*mat).second_level);
	free((*mat).data);
}
