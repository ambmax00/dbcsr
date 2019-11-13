#ifndef DBCSR_TENSOR_H
#define DBCSR_TENSOR_H

#include <mpi.h>
#include <stdbool.h> /* we need bool from C99 */

#:include "tensors/dbcsr_tensor.fypp"
#:set ndims = range(2,maxrank+1)

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

	void c_dbcsr_t_contract(int* c_alpha, void* c_tensor_1, void* c_tensor_2, int* c_beta, void* c_tensor_3, 
                               int* c_contract_1, int c_contract_1_size,
                               int* c_notcontract_1, int c_notcontract_1_size,
                               int* c_contract_2, int c_contract_2_size,
                               int* c_notcontract_2, int c_notcontract_2_size,
                               int* c_map_1, int c_map_1_size,
                               int* c_map_2, int c_map_2_size,
                               int* c_bounds_1, int c_bounds_1_size,
                               int* c_bounds_2, int c_bounds_2_size,
                               int* c_bounds_3, int c_bounds_3_size,
                               bool* c_optimize_dist, 
                               void** c_pgrid_opt_1, 
                               void** c_pgrid_opt_2, 
                               void** c_pgrid_opt_3,
                               double* filter_eps, long long int* flop, 
                               int* move_data, int* unit_nr, bool* log_verbose);
                               
    int c_ndims_tensor(void* tensor);
                               
#:for dsuffix, ctype in c_dtype_float_list
#:for ndim in ndims    
                    
      void c_dbcsr_t_get_${ndim}$d_block_${dsuffix}$ (void* c_tensor, int tensor_dim, int* c_ind, int* c_sizes, 
		${ctype}$** c_block, bool* c_found);
		
	  void c_dbcsr_t_put_${ndim}$d_block_${dsuffix}$ (void* c_tensor, int tensor_dim, int* c_ind, int*c_sizes, 
        ${ctype}$* c_block, bool* c_summation, ${ctype}$* c_scale);
	
#:endfor
#:endfor      	
	
#ifdef __cplusplus
}
#endif

#:for dsuffix, ctype in c_dtype_float_list
static void c_dbcsr_t_get_block(void* c_tensor, int* c_ind, int* c_sizes, 
		${ctype}$** c_block, bool* c_found) {
	
	int tensor_dim = c_ndims_tensor(c_tensor);
	
	switch(tensor_dim) {
		#:for ndim in ndims
		case ${ndim}$: c_dbcsr_t_get_${ndim}$d_block_${dsuffix}$ (c_tensor, tensor_dim, 
		c_ind, c_sizes, c_block, c_found);
		#:endfor
	}
	
}
#:endfor

#:for dsuffix, ctype in c_dtype_float_list
static void c_dbcsr_t_put_block(void* c_tensor, int* c_ind, int* c_sizes, 
		${ctype}$* c_block, bool* c_summation, ${ctype}$* c_scale) {
	
	int tensor_dim = c_ndims_tensor(c_tensor);
	
	switch(tensor_dim) {
		#:for ndim in ndims
		case ${ndim}$: c_dbcsr_t_put_${ndim}$d_block_${dsuffix}$ (c_tensor, tensor_dim, 
		c_ind, c_sizes, c_block, c_summation, c_scale);
		#:endfor
	}
	
}
#:endfor

#endif // DBCSR_H
