/**
 * This code has no warranty, it is free
 * to use and distribute.
 * It is a study exercise and has to be taken like that.
 *
 * 11-11-2021
 * Author: Francesco Tomba, DSSC student @ UNITS
 * */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <limits.h>
#include "tensor_utils.h"

int main(int argc, char** argv){
	size_t nx,ny,nz;

	if(argc != 4){
		printf("usage ./exec: [N_X] [N_Y] [N_Z]\n");
		return 1;	
	}
	else{
		nx = atoll(argv[1]);
		ny = atoll(argv[2]);
		nz = atoll(argv[3]);
		printf("Got matrix dimensions %lu X %lu X %lu\n", nx, ny, nz);
	}
	
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;		
	int mpi_rank;
	int nnodes;
	MPI_Comm_rank(comm,&mpi_rank);


	//creating 3d communicator with periodic boundary conditions
	MPI_Comm grid_3d;
	int dims[3];
	int periods[3] = {1,1,1};
	MPI_Comm_size(comm,&nnodes);
	MPI_Dims_create(nnodes,3,dims);
	MPI_Cart_create(comm,3,dims,periods,1,&grid_3d);

	int my_coords[3];
	MPI_Cart_coords(grid_3d,mpi_rank, 3,  my_coords);
	nx = dims[0]*nx;
	ny = dims[1]*ny;
	nz = dims[2]*nz;
	printf("Hello, world from processor %d in coordinates: [%d %d %d]\n",mpi_rank,my_coords[0], my_coords[1], my_coords[2]);
	Matrix3D_double mat;
	if(mpi_rank == 0){
		int v = 0;
		allocate_matrix3D(&mat,nx,ny,nz);
		for(size_t i = 0; i < nx; ++i){
			
			for(size_t j = 0; j < ny; ++j){
				
				for(size_t k = 0; k < nz; ++k){
					mat.data[i][j][k] = v;// rand()/(double)INT_MAX;
					++v;
				}
			}
		}
	//	printf("[1,1,1] element in matrix 3D for processor %d is %lf\n",mpi_rank,mat.data[1][1][1]);
	}	
	
	size_t my_nx, my_ny, my_nz;
	size_t sp_x, sp_y, sp_z;
	get_matrix_dimension(&my_nx, &sp_x, my_coords[0], dims[0], nx);
	get_matrix_dimension(&my_ny, &sp_y, my_coords[1], dims[1], ny);
	get_matrix_dimension(&my_nz, &sp_z, my_coords[2], dims[2], nz);



	printf("From processor %d dimensions are: [%ld %ld %ld]\n",mpi_rank,my_nx,my_ny,my_nz);

	MPI_Datatype mat_block;
	//Not so efficient, rather better using size_t type but mpi for reasons raises 
	//warings
	int* block_lenghts = (int*)malloc(my_nx*my_ny*sizeof(int));
	int* displacements = (int*)malloc(my_nx*my_ny*sizeof(int));

	for(size_t i = 0; i < my_nx; ++i){
		
		for(size_t j = 0; j < my_ny; ++j){
			displacements[i*my_ny+j] = ((i+sp_x)*ny + j + sp_y)*nz + sp_z;
			block_lenghts[i*my_ny+j] = my_nz;
		}
	}

	MPI_Type_indexed(my_nx*my_ny, block_lenghts, displacements, MPI_DOUBLE, &mat_block);

	MPI_Type_commit(&mat_block);
	
	double* matrix_chunk = (double*)malloc(my_nx*my_ny*my_nz*sizeof(double));
	MPI_Scatter(mat.last_level, 1, mat_block, matrix_chunk, my_nx*my_ny*my_nz, MPI_DOUBLE, 0, grid_3d); 
	printf("proc: %d recieved in [1,1,1] %lf\n",mpi_rank, matrix_chunk[1]);
	free(matrix_chunk);
	//free_matrix3D(&mat);	
	MPI_Finalize();	
	return 0;
}
