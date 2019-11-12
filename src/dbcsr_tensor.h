#ifndef DBCSR_TENSOR_H
#define DBCSR_TENSOR_H

#include <mpi.h>
#include <stdbool.h> /* we need bool from C99 */

#:include "tensors/dbcsr_tensor.fypp"

#ifdef __cplusplus
extern "C" {
#endif
    void c_dbcsr_t_pgrid_create(MPI_Fint* fcomm, int* c_dims, int dims_size, void** c_pgrid, 
		int* c_map1_2d, int map1_2d_size, int* c_map2_2d, int map2_2d_size, int* nsplit, int* dimsplit);
																			   //^            ^
																			   // needed by reference, to check if present 
	void c_dbcsr_t_distribution_new(void** c_dist, void* c_pgrid, int* c_map1_2d, int map1_2d_size,
		int* c_map2_2d, int map2_2d_size, ${extern_varlist_and_size("c_nd_dist")}$, bool* own_comm);
		
	void c_dbcsr_t_create_new(void** c_tensor, const char* c_name, void* c_dist, int* c_map1_2d, int c_map1_2d_size,
       int* c_map2_2d, int c_map2_2d_size, int* data_type, ${extern_varlist_and_size("c_blk_size")}$);

		
	
#ifdef __cplusplus
}
#endif

#endif // DBCSR_H
