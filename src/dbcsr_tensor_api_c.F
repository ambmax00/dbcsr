#:include "tensors/dbcsr_tensor.fypp"
#:set ndims = range(2,maxrank+1)
#:set ddims = range(1,maxrank+1)

MODULE dbcsr_tensor_api_c

   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_loc, c_ptr, c_double, C_NULL_CHAR, &
                              c_float, c_f_pointer, c_int, c_long_long, &
                              c_char, c_null_ptr, c_bool, &
                              c_float_complex, c_double_complex
   USE dbcsr_tensor_api
   USE dbcsr_kinds, ONLY: &
      ${uselist(dtype_float_prec)}$, default_string_length, int_8
   USE dbcsr_data_types, ONLY: dbcsr_scalar_type
   USE dbcsr_machine, ONLY: default_output_unit

   IMPLICIT NONE
   PRIVATE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param c_str ...
!> \param str ...
! **************************************************************************************************
   SUBROUTINE c_f_string(c_str, str)
      USE, INTRINSIC :: iso_c_binding, ONLY: c_ptr, c_f_pointer, c_char
      TYPE(c_ptr), INTENT(in) :: c_str
      CHARACTER(kind=c_char), POINTER :: arr(:)
      CHARACTER(:, kind=c_char), ALLOCATABLE, INTENT(out) :: str
      INTEGER(8) :: n, i
      INTERFACE
         ! steal std c library function rather than writing our own.
         FUNCTION strlen(s) bind(c, name='strlen')
            USE, INTRINSIC :: iso_c_binding, ONLY: c_ptr, c_size_t
            IMPLICIT NONE
            !----
            TYPE(c_ptr), INTENT(in), value :: s
            INTEGER(c_size_t) :: strlen
         END FUNCTION strlen
      END INTERFACE
      n = strlen(c_str)
      !****
      CALL c_f_pointer(c_str, arr, [n])
      ALLOCATE (CHARACTER(len=n) :: str)
      DO i = 1, n
         str(i:i) = arr(i)
      ENDDO
   END SUBROUTINE c_f_string

   SUBROUTINE c_dbcsr_t_pgrid_create(fcomm, c_dims, dims_size, c_pgrid, &
                                     c_map1_2d, map1_2d_size, c_map2_2d, map2_2d_size, nsplit, dimsplit) &
      bind(C, name="c_dbcsr_t_pgrid_create")

      INTEGER(kind=c_int), INTENT(in)                    :: fcomm
      INTEGER(kind=c_int), INTENT(in), value             :: dims_size
      TYPE(c_ptr), INTENT(out)                           :: c_pgrid
      INTEGER(kind=c_int), INTENT(inout), TARGET         :: c_dims(dims_size)
      INTEGER(kind=c_int), INTENT(in), VALUE             :: map1_2d_size
      INTEGER(kind=c_int), INTENT(in), TARGET, OPTIONAL  :: c_map1_2d(map1_2d_size)
      INTEGER(kind=c_int), INTENT(in), VALUE             :: map2_2d_size
      INTEGER(kind=c_int), INTENT(in), TARGET, OPTIONAL  :: c_map2_2d(map2_2d_size)
      INTEGER(kind=c_int), INTENT(in), OPTIONAL          :: nsplit, dimsplit

      INTEGER, POINTER                                   :: map1_2d(:), map2_2d(:), dims(:)
      TYPE(dbcsr_t_pgrid_type), POINTER                  :: pgrid

      dims => c_dims

      PRINT *, "SIZE: ", map1_2d_size

      IF (PRESENT(c_map1_2d) .AND. PRESENT(c_map2_2d)) THEN
         map1_2d => c_map1_2d
         map2_2d => c_map2_2d
         PRINT *, "PRESENT"
      ELSE
         map1_2d => null()
         map2_2d => null()
         PRINT *, "NOT PRESENT"
      ENDIF

      ALLOCATE (pgrid)

      CALL dbcsr_t_pgrid_create(fcomm, dims, pgrid, map1_2d, map2_2d, &
                                nsplit, dimsplit)

      c_pgrid = c_loc(pgrid)

      PRINT *, "PGRID", pgrid%mp_comm_2d

   END SUBROUTINE

! **************************************************************************************************
!> \brief destroy process grid
!> \param pgrid ...
!> \param keep_comm if .TRUE. communicator is not freed
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_pgrid_destroy(c_pgrid, c_keep_comm) &
   BIND(C, name="c_dbcsr_t_pgrid_destroy")
   
      TYPE(c_ptr), INTENT(IN), VALUE :: c_pgrid
      LOGICAL(kind=c_bool), INTENT(IN), OPTIONAL :: c_keep_comm
      TYPE(dbcsr_t_pgrid_type), POINTER :: pgrid
      LOGICAL, ALLOCATABLE           :: keep_comm
      
      
      CALL c_f_pointer(c_pgrid, pgrid)
      IF (PRESENT(c_keep_comm)) keep_comm = c_keep_comm
      
      CALL dbcsr_t_pgrid_destroy(pgrid, keep_comm)
      
   END SUBROUTINE
   

   SUBROUTINE c_dbcsr_t_distribution_new(c_dist, c_pgrid, c_map1_2d, map1_2d_size, &
                                         c_map2_2d, map2_2d_size, ${c_varlist_and_size("c_nd_dist")}$, c_own_comm) &
      bind(C, name="c_dbcsr_t_distribution_new")

      TYPE(c_ptr), INTENT(OUT)                        :: c_dist
      TYPE(c_ptr), INTENT(IN), value                  :: c_pgrid

      INTEGER(kind=c_int), INTENT(in), value             :: map1_2d_size
      INTEGER(kind=c_int), INTENT(in), TARGET            :: c_map1_2d(map1_2d_size)
      INTEGER(kind=c_int), INTENT(in), value             :: map2_2d_size
      INTEGER(kind=c_int), INTENT(in), TARGET            :: c_map2_2d(map2_2d_size)

      INTEGER(kind=c_int), INTENT(in), value             :: ${varlist("c_nd_dist", suffix="_size")}$
      INTEGER(kind=c_int), INTENT(IN), TARGET, OPTIONAL  :: ${c_arrlist_and_size("c_nd_dist")}$
      LOGICAL(kind=c_bool), INTENT(IN), OPTIONAL         :: c_own_comm

      INTEGER, POINTER                                   :: ${pointer_list("nd_dist")}$
      INTEGER, POINTER                                   :: map1_2d(:), map2_2d(:)
      TYPE(dbcsr_t_pgrid_type), POINTER                  :: pgrid
      TYPE(dbcsr_t_distribution_type), POINTER           :: dist
      LOGICAL, allocatable                               :: own_comm
      INTEGER :: i

      map1_2d => c_map1_2d
      map2_2d => c_map2_2d

#:for ddim in ddims
      IF (PRESENT(c_nd_dist_${ddim}$)) THEN
         nd_dist_${ddim}$ => c_nd_dist_${ddim}$
      ELSE
         nd_dist_${ddim}$ => null()
      ENDIF

#:endfor

      IF (PRESENT(c_own_comm)) THEN
         own_comm = c_own_comm
         PRINT *, "Bool is there."
      ENDIF

      DO i = 1, SIZE(c_nd_dist_1)
        PRINT*, c_nd_dist_1(i)
      END DO

      PRINT *, "HERE"

      ALLOCATE (dist)

      CALL c_f_pointer(c_pgrid, pgrid)
      
      

      PRINT *, "ENTERING"

      CALL dbcsr_t_distribution_new(dist, pgrid, map1_2d, map2_2d, &
                                    ${varlist("nd_dist")}$, own_comm)

      PRINT *, "EXITING"

      c_dist = c_loc(dist)

   END SUBROUTINE
   
! **************************************************************************************************
!> \brief Destroy tensor distribution
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_distribution_destroy(c_dist) BIND(C, name="c_dbcsr_t_distribution_destroy")
      TYPE(c_ptr), INTENT(IN), VALUE :: c_dist
      TYPE(dbcsr_t_distribution_type), POINTER :: dist
      
      CALL c_f_pointer(c_dist, dist);
      
      CALL dbcsr_t_distribution_destroy(dist)
      
   END SUBROUTINE
  

   SUBROUTINE c_dbcsr_t_create_new(c_tensor, c_name, c_dist, c_map1_2d, map1_2d_size, &
                                   c_map2_2d, map2_2d_size, c_data_type, ${c_varlist_and_size("c_blk_size")}$) &
      bind(C, name="c_dbcsr_t_create_new")

      TYPE(c_ptr), INTENT(OUT)                          :: c_tensor
      TYPE(c_ptr), INTENT(IN), value                    :: c_dist, c_name

      INTEGER(kind=c_int), INTENT(in), value             :: map1_2d_size
      INTEGER(kind=c_int), INTENT(in), TARGET            :: c_map1_2d(map1_2d_size)
      INTEGER(kind=c_int), INTENT(in), value             :: map2_2d_size
      INTEGER(kind=c_int), INTENT(in), TARGET            :: c_map2_2d(map2_2d_size)

      INTEGER(kind=c_int), INTENT(IN), OPTIONAL          :: c_data_type

      INTEGER(kind=c_int), INTENT(in), value             :: ${varlist("c_blk_size", suffix="_size")}$
      INTEGER(kind=c_int), INTENT(IN), TARGET, OPTIONAL  :: ${c_arrlist_and_size("c_blk_size")}$

      CHARACTER(:, kind=c_char), ALLOCATABLE             :: fname

      INTEGER, POINTER                                  :: map1_2d(:), map2_2d(:)
      INTEGER, POINTER                                  :: ${pointer_list("blk_size")}$
      TYPE(dbcsr_t_type), POINTER                        :: tensor
      TYPE(dbcsr_t_distribution_type), POINTER           :: dist
      INTEGER, ALLOCATABLE :: data_type

      map1_2d => c_map1_2d
      map2_2d => c_map2_2d

#:for ddim in ddims
      IF (PRESENT(c_blk_size_${ddim}$)) THEN
         blk_size_${ddim}$ => c_blk_size_${ddim}$
      ELSE
         blk_size_${ddim}$ => null()
      ENDIF

#:endfor

      IF (PRESENT(c_data_type)) THEN
        data_type = c_data_type
      END IF

      ALLOCATE (tensor)

      CALL c_f_pointer(c_dist, dist)
      CALL c_f_string(c_name, fname)

      CALL dbcsr_t_create(tensor, fname, dist, map1_2d, map2_2d, data_type, &
                          ${varlist("blk_size")}$)

      c_tensor = c_loc(tensor)

   END SUBROUTINE
   
! **************************************************************************************************
!> \brief Destroy a tensor
!> \param tensor ...
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_destroy(c_tensor) BIND(C, name="c_dbcsr_t_destroy")
      TYPE(c_ptr), INTENT(IN), VALUE :: c_tensor
      TYPE(dbcsr_t_type), POINTER :: tensor
      
      CALL c_f_pointer(c_tensor, tensor)
      
      CALL dbcsr_t_destroy(tensor)
      
   END SUBROUTINE

   ! **************************************************************************************************
!> \brief Contract tensors by multiplying matrix representations.
!>        tensor_3(map_1, map_2) := alpha * tensor_1(notcontract_1, contract_1)
!>                                        * tensor_2(contract_2, notcontract_2)
!>                                 + beta * tensor_3(map_1, map_2)
!> \param alpha ...
!> \param tensor_1 first tensor (in)
!> \param tensor_2 second tensor (in)
!> \param beta ...
!> \param tensor_3 contracted tensor (out)
!> \param contract_1 indices of tensor_1 to contract
!> \param notcontract_1 indices of tensor_1 not to contract
!> \param contract_2 indices of tensor_2 to contract (1:1 with contract_1)
!> \param notcontract_2 indices of tensor_2 not to contract
!> \param map_1 which indices of tensor_3 map to non-contracted indices of tensor_1 (1:1 with notcontract_1)
!> \param map_2 which indices of tensor_3 map to non-contracted indices of tensor_2 (1:1 with notcontract_2)
!> \param bounds_1 bounds corresponding to contract_1 AKA contract_2
!> \param bounds_2 bounds corresponding to notcontract_1
!> \param bounds_3 bounds corresponding to notcontract_2
!> \param optimize_dist Whether distribution should be optimized internally. In the current implementation
!>                      this guarantees optimal parameters only for dense matrices.
!> \param pgrid_opt_1 Optionally return optimal process grid for tensor_1. This can be used to choose
!>                    optimal process grids for subsequent tensor contractions with tensors of similar
!>                    shape and sparsity. Under some conditions, pgrid_opt_1 can not be returned, in this
!>                    case the pointer is not associated.
!> \param pgrid_opt_2 Optionally return optimal process grid for tensor_2.
!> \param pgrid_opt_3 Optionally return optimal process grid for tensor_3.
!> \param filter_eps As in DBCSR mm
!> \param flop As in DBCSR mm
!> \param move_data memory optimization: transfer data such that tensor_1 and tensor_2 may be empty
!>        on return
!> \param unit_nr output unit for logging
!> \param log_verbose verbose logging (for testing only)
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_contract(c_alpha, c_tensor_1, c_tensor_2, c_beta, c_tensor_3, &
                                 c_contract_1, contract_1_size, &
                                 c_notcontract_1, notcontract_1_size, &
                                 c_contract_2, contract_2_size, &
                                 c_notcontract_2, notcontract_2_size, &
                                 c_map_1, map_1_size, c_map_2, map_2_size, &
                                 c_bounds_1, bounds_1_size, c_bounds_2, bounds_2_size, &
                                 c_bounds_3, bounds_3_size, &
                                 c_optimize_dist, c_pgrid_opt_1, c_pgrid_opt_2, c_pgrid_opt_3, &
                                 c_filter_eps, c_flop, c_move_data, c_unit_nr, c_log_verbose) &
      BIND(C, name="c_dbcsr_t_contract")

      TYPE(c_ptr), INTENT(IN)  ::  c_alpha
      TYPE(c_ptr), INTENT(INOUT) :: c_tensor_1, c_tensor_2
      TYPE(c_ptr), INTENT(IN) :: c_beta
      TYPE(c_ptr), INTENT(INOUT)  :: c_tensor_3

      INTEGER(kind=c_int), INTENT(IN), VALUE :: contract_1_size, contract_2_size
      INTEGER(kind=c_int), INTENT(IN), VALUE :: notcontract_1_size, notcontract_2_size
      INTEGER(kind=c_int), INTENT(IN), VALUE :: map_1_size, map_2_size, bounds_1_size, &
                                                                         bounds_2_size, bounds_3_size

      INTEGER(kind=c_int), INTENT(IN), TARGET :: c_contract_1(contract_1_size), c_contract_2(contract_2_size)
      INTEGER(kind=c_int), INTENT(IN), TARGET :: c_map_1(map_1_size), c_map_2(map_2_size)
      INTEGER(kind=c_int), INTENT(IN), TARGET :: c_notcontract_1(notcontract_1_size), c_notcontract_2(notcontract_2_size)

      INTEGER(kind=c_int), INTENT(IN), DIMENSION(2, contract_1_size), OPTIONAL, TARGET :: c_bounds_1
      INTEGER(kind=c_int), INTENT(IN), DIMENSION(2, notcontract_1_size), OPTIONAL, TARGET :: c_bounds_2
      INTEGER(kind=c_int), INTENT(IN), DIMENSION(2, notcontract_2_size), OPTIONAL, TARGET :: c_bounds_3

      LOGICAL(kind=c_bool), INTENT(IN), OPTIONAL :: c_optimize_dist
      TYPE(c_ptr), INTENT(OUT), OPTIONAL :: c_pgrid_opt_1, c_pgrid_opt_2, c_pgrid_opt_3

      REAL(kind=c_double), INTENT(IN), OPTIONAL  :: c_filter_eps
      INTEGER(kind=c_long_long), INTENT(IN), OPTIONAL :: c_flop
      LOGICAL(kind=c_bool), INTENT(IN), OPTIONAL :: c_move_data
      INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: c_unit_nr
      LOGICAL(kind=c_bool), INTENT(IN), OPTIONAL :: c_log_verbose

      TYPE(dbcsr_scalar_type), POINTER           :: alpha
      TYPE(dbcsr_t_type), POINTER                 :: tensor_1, tensor_2
      TYPE(dbcsr_scalar_type), POINTER           :: beta
      INTEGER, DIMENSION(:), POINTER              :: contract_1, contract_2, map_1, map_2
      INTEGER, DIMENSION(:), POINTER              :: notcontract_1, notcontract_2
      TYPE(dbcsr_t_type), POINTER             :: tensor_3
      INTEGER, DIMENSION(:, :), POINTER          :: bounds_1 => null()
      INTEGER, DIMENSION(:, :), POINTER          :: bounds_2 => null()
      INTEGER, DIMENSION(:, :), POINTER          :: bounds_3 => null()
      LOGICAL, POINTER                  :: optimize_dist => null()
      TYPE(dbcsr_t_pgrid_type), POINTER     :: pgrid_opt_1 => null(), &
                                               pgrid_opt_2 => null(), pgrid_opt_3 => null()

      REAL(KIND=real_8), ALLOCATABLE        :: filter_eps
      INTEGER(KIND=int_8), ALLOCATABLE    :: flop
      LOGICAL, ALLOCATABLE                 :: move_data
      INTEGER, ALLOCATABLE                 :: unit_nr
      LOGICAL, ALLOCATABLE                  :: log_verbose

      CALL c_f_pointer(c_alpha, alpha)
      CALL c_f_pointer(c_tensor_1, tensor_1)
      CALL c_f_pointer(c_tensor_2, tensor_2)
      CALL c_f_pointer(c_beta, beta)
      CALL c_f_pointer(c_tensor_3, tensor_3)

      contract_1 => c_contract_1
      contract_2 => c_contract_2
      map_1 => c_map_1
      map_2 => c_map_2
      notcontract_1 => c_notcontract_1
      notcontract_2 => c_notcontract_2

#:for dim in range(1,4)

      IF (PRESENT(c_bounds_${dim}$)) THEN
         bounds_${dim}$ => c_bounds_${dim}$
      ENDIF

      ALLOCATE (pgrid_opt_${dim}$)

#:endfor

#:set list = ['filter_eps', 'move_data', 'unit_nr', 'log_verbose']
#:for var in list

      IF (PRESENT(c_${var}$)) THEN
         ${var}$ = c_${var}$
      ENDIF

#:endfor

      CALL dbcsr_t_contract(alpha, tensor_1, tensor_2, beta, tensor_3, &
                            contract_1, notcontract_1, &
                            contract_2, notcontract_2, &
                            map_1, map_2, &
                            bounds_1, bounds_2, bounds_3, &
                            optimize_dist, pgrid_opt_1, pgrid_opt_2, pgrid_opt_3, &
                            filter_eps, flop, move_data, unit_nr, log_verbose)

#:for dim in range(1,4)

      c_pgrid_opt_${dim}$ = c_loc(pgrid_opt_${dim}$)

#:endfor

   END SUBROUTINE
   
   FUNCTION c_ndims_tensor(c_tensor) BIND(C, name="c_ndims_tensor")
      TYPE(c_ptr), INTENT(IN), VALUE :: c_tensor
      TYPE(dbcsr_t_type), POINTER :: tensor
      INTEGER :: c_ndims_tensor

      CALL c_f_pointer(c_tensor, tensor)
       
      c_ndims_tensor = tensor%nd_index%ndim_nd
   END FUNCTION
   
#:for dparam, dtype, dsuffix, basetype, cftype in cf_dtype_float_list
#:for ndim in ndims
! **************************************************************************************************
!> \brief allocate and get block
!> \param tensor ...
!> \param ind block index
!> \param block block to get
!> \param found whether block was found
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_allocate_and_get_${ndim}$d_block_${dsuffix}$ (c_tensor, tensor_dim, c_ind, c_block, c_found) &
   bind(C, name="c_dbcsr_t_allocate_and_get_${ndim}$d_block_${dsuffix}$")
      
       TYPE(c_ptr), INTENT(IN), VALUE :: c_tensor
       INTEGER(kind=c_int), INTENT(IN), VALUE :: tensor_dim
       INTEGER(kind=c_int), DIMENSION(tensor_dim), INTENT(IN), TARGET :: c_ind
       TYPE(c_ptr), INTENT(OUT) :: c_block
       TYPE(c_ptr), INTENT(OUT) :: c_found
        
       TYPE(dbcsr_t_type), POINTER :: tensor
       INTEGER, DIMENSION(:), POINTER :: ind
       ${dtype}$, DIMENSION(${shape_colon(ndim)}$), ALLOCATABLE, TARGET :: block
       LOGICAL, POINTER :: found
       
       ind => c_ind
       !block => c_block
       
       CALL c_f_pointer(c_tensor, tensor)
       ALLOCATE(found)
   
       CALL dbcsr_t_get_block(tensor, ind+1, block, found) 
   
       c_block = c_loc(block)
       c_found = c_loc(found)

   END SUBROUTINE
#:endfor
#:endfor
   
#:for dparam, dtype, dsuffix, basetype, cftype in cf_dtype_float_list
#:for ndim in ndims
! **************************************************************************************************
!> \brief Template for dbcsr_t_get_block.
!> \param tensor ...
!> \param ind block index
!> \param sizes block size
!> \param block block to get
!> \param found whether block was found
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_get_${ndim}$d_block_${dsuffix}$ (c_tensor, tensor_dim, c_ind, c_sizes, c_block, c_found) &
   BIND(C, name="c_dbcsr_t_get_${ndim}$d_block_${dsuffix}$")
   
       TYPE(c_ptr), INTENT(IN), VALUE :: c_tensor
       INTEGER(kind=c_int), INTENT(IN), VALUE :: tensor_dim
       INTEGER(kind=c_int), DIMENSION(tensor_dim), INTENT(IN), TARGET :: c_ind
       INTEGER(kind=c_int), DIMENSION(tensor_dim), INTENT(IN), TARGET :: c_sizes
       ${basetype}$(kind=${cftype}$), DIMENSION(${arrlist("c_sizes", nmax=ndim)}$), &
            TARGET, INTENT(OUT) :: c_block
       TYPE(c_ptr), INTENT(OUT) :: c_found
        
       TYPE(dbcsr_t_type), POINTER :: tensor
       INTEGER, DIMENSION(:), POINTER :: ind
       INTEGER, DIMENSION(:), POINTER :: sizes
       ${dtype}$, DIMENSION(${shape_colon(ndim)}$), POINTER  :: block
       LOGICAL, POINTER :: found
       
       ind => c_ind
       sizes => c_sizes
       block => c_block
       
       CALL c_f_pointer(c_tensor, tensor)
       ALLOCATE(found)
   
       CALL dbcsr_t_get_block(tensor, ind+1, sizes, block, found) 
   
       c_found = c_loc(found)
   
   END SUBROUTINE
#:endfor
#:endfor

#:for dparam, dtype, dsuffix, basetype, cftype in cf_dtype_float_list
#:for ndim in ndims
! **************************************************************************************************
!> \brief Template for dbcsr_t_put_block.
!> \param tensor ...
!> \param ind block index
!> \param sizes block size
!> \param block block to put
!> \param summation whether block should be summed to existing block
!> \param scale scaling factor
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_put_${ndim}$d_block_${dsuffix}$ (c_tensor, tensor_dim, c_ind, c_sizes, &
        c_block, c_summation, c_scale) BIND(C, name="c_dbcsr_t_put_${ndim}$d_block_${dsuffix}$")
        
        
      TYPE(c_ptr), INTENT(IN), VALUE  :: c_tensor
      INTEGER(kind=c_int), INTENT(IN), VALUE :: tensor_dim
      INTEGER(kind=c_int), DIMENSION(tensor_dim), INTENT(IN), TARGET :: c_ind
      INTEGER(kind=c_int), DIMENSION(tensor_dim), INTENT(IN), TARGET :: c_sizes
      ${basetype}$(kind=${cftype}$), INTENT(IN), &
        DIMENSION(${arrlist("c_sizes", nmax=ndim)}$), TARGET :: c_block
      LOGICAL(kind=c_bool), INTENT(IN), OPTIONAL :: c_summation 
      ${basetype}$(kind=${cftype}$), INTENT(IN), OPTIONAL :: c_scale
       
      TYPE(dbcsr_t_type), POINTER :: tensor
      INTEGER, DIMENSION(:), POINTER :: ind
      INTEGER, DIMENSION(:), POINTER :: sizes
      ${dtype}$, DIMENSION(${shape_colon(ndim)}$), POINTER  :: block
      LOGICAL, ALLOCATABLE :: summation
      ${dtype}$, ALLOCATABLE :: scale
      
      INTEGER :: i
      
      PRINT*, "EVERYTHING IS DEFINED:"
      
      ind => c_ind
      
      IF (PRESENT(c_summation)) THEN
        summation = c_summation
      ENDIF
      
      IF (PRESENT(c_scale)) THEN
        scale = c_scale
      ENDIF
      
      CALL c_f_pointer(c_tensor, tensor) 
      sizes => c_sizes
      block => c_block
      
      PRINT*, "IND"
      DO i = 1, SIZE(ind)
        PRINT*, ind(i)
      END DO
      
      PRINT*, "SIZES"
      DO i = 1, SIZE(sizes)
        PRINT*, sizes(i)
      END DO
    
      PRINT*, "CALLING..."
      
      CALL dbcsr_t_put_block(tensor, ind+1, sizes, &
        block, summation, scale)
   
   END SUBROUTINE
#:endfor
#:endfor

! **************************************************************************************************
!> \brief Generalization of dbcsr_get_stored_coordinates for tensors.
!> \param tensor ...
!> \param ind_nd ...
!> \param processor ...
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_get_stored_coordinates(c_tensor, tensor_dim, c_ind_nd, c_processor) &
   BIND(C, name="c_dbcsr_t_get_stored_coordinates")
   
      TYPE(c_ptr), INTENT(IN), VALUE :: c_tensor
      INTEGER(kind=c_int), INTENT(IN), VALUE :: tensor_dim
      INTEGER(kind=c_int), DIMENSION(tensor_dim), INTENT(IN), TARGET :: c_ind_nd
      INTEGER(kind=c_int), INTENT(OUT) :: c_processor

      TYPE(dbcsr_t_type), POINTER :: tensor
      INTEGER, DIMENSION(:), POINTER  :: ind_nd
      
       INTEGER, DIMENSION(tensor_dim) :: nblks_total, nfull_total, &
                                                                         nblks_local, nfull_local, &
                                                                         my_ploc, pdims
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("blks_local")}$
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("proc_dist")}$
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("blk_size")}$
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("blk_offset")}$
      TYPE(dbcsr_t_distribution_type)         :: distribution

      INTEGER, DIMENSION(tensor_dim)                  :: pdims_tmp, my_ploc_tmp
      INTEGER :: i
        
      CALL c_f_pointer(c_tensor, tensor)
      ind_nd => c_ind_nd
      
      PRINT*, "HEEEEEY"
      
      CALL dbcsr_t_get_info(tensor, nblks_total, &
                               nfull_total, &
                               nblks_local, &
                               nfull_local, &
                               pdims, &
                               my_ploc, &
                               ${varlist("blks_local")}$, &
                               ${varlist("proc_dist")}$, &
                               ${varlist("blk_size")}$, &
                               ${varlist("blk_offset")}$, &
                               distribution)
                               
       PRINT*, "nblks_total"
       DO i = 1, SIZE(nblks_total)
        PRINT*, nblks_total(i)
       END DO
       
       PRINT*, "nfull_total"
       DO i = 1, SIZE(nfull_total)
        PRINT*, nfull_total(i)
       END DO
       
       PRINT*, "nblks_local"
       DO i = 1, SIZE(nblks_local)
        PRINT*, nblks_local(i)
       END DO
       
       PRINT*, "nfull_local"
       DO i = 1, SIZE(nfull_local)
        PRINT*, nfull_local(i)
       END DO
       
       PRINT*, "pdims"
       DO i = 1, SIZE(pdims)
        PRINT*, pdims(i)
       END DO
       
       PRINT*, "my_ploc"
       DO i = 1, SIZE(my_ploc)
        PRINT*, my_ploc(i)
       END DO
       
#:for ddim in ddims

       PRINT*, "blks_local_${ddim}$"
       DO i = 1, SIZE(blks_local_${ddim}$)
        PRINT*, blks_local_${ddim}$(i)
       END DO
       
       PRINT*, "proc_dist_${ddim}$"
       DO i = 1, SIZE(proc_dist_${ddim}$)
        PRINT*, proc_dist_${ddim}$(i)
       END DO
       
       PRINT*, "blk_size_${ddim}$"
       DO i = 1, SIZE(blk_size_${ddim}$)
        PRINT*, blk_size_${ddim}$(i)
       END DO
       
       PRINT*, "blk_offset_${ddim}$"
       DO i = 1, SIZE(blk_offset_${ddim}$)
        PRINT*, blk_offset_${ddim}$(i)
       END DO
  
#:endfor    
        
       CALL dbcsr_t_get_stored_coordinates(tensor, ind_nd+1, c_processor)
      
   END SUBROUTINE
   
! **************************************************************************************************
!> \brief reserve tensor blocks using block indices
!> \param tensor ...
!> \param blk_ind_1 index of blocks to reserve in each dimension
!> \param blk_ind_2 ...
!> \param blk_ind_3 ...
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_reserve_blocks_index(c_tensor, ${c_varlist_and_size("c_blk_ind")}$) &
   BIND(C, name="c_dbcsr_t_reserve_blocks_index")
      
      TYPE(c_ptr), INTENT(IN), VALUE :: c_tensor
      INTEGER(kind=c_int), INTENT(in), value             :: ${varlist("c_blk_ind", suffix="_size")}$
      INTEGER(kind=c_int), INTENT(IN), TARGET, OPTIONAL  :: ${c_arrlist_and_size("c_blk_ind")}$
      
      TYPE(dbcsr_t_type), POINTER           :: tensor
      INTEGER, DIMENSION(:), POINTER :: ${varlist("blk_ind")}$
      
    
#:for ddim in ddims
      IF (PRESENT(c_blk_ind_${ddim}$)) THEN
         blk_ind_${ddim}$ => c_blk_ind_${ddim}$
      ELSE
         blk_ind_${ddim}$ => null()
      ENDIF
#:endfor

      CALL c_f_pointer(c_tensor, tensor)

      CALL dbcsr_t_reserve_blocks(tensor, ${varlist("blk_ind", suffix=" +1")}$)
    
   END SUBROUTINE  
   
! **************************************************************************************************
!> \brief Number of dimensions.
!> \param iterator ...
!> \return ...
!> \par Note: specification function below must be defined before it is used in
!>      the source due to a bug in the IBM XL Fortran compiler (compilation fails)
! **************************************************************************************************
   FUNCTION c_ndims_iterator(c_iterator) BIND(C, name="c_ndims_iterator")
   
      TYPE(c_ptr), INTENT(IN), VALUE :: c_iterator
   
      TYPE(dbcsr_t_iterator_type), POINTER :: iterator
      INTEGER(kind=c_int) :: c_ndims_iterator

      CALL c_f_pointer(c_iterator, iterator)

      c_ndims_iterator = iterator%nd_index%ndim_nd
      
   END FUNCTION   
   
! **************************************************************************************************
!> \brief Generalization of dbcsr_iterator_start for tensors.
!> \param iterator ...
!> \param tensor ...
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_iterator_start(c_iterator, c_tensor) BIND(C, name="c_dbcsr_t_iterator_start")
   
      TYPE(c_ptr), INTENT(OUT) :: c_iterator
      TYPE(c_ptr), INTENT(IN), VALUE :: c_tensor
   
      TYPE(dbcsr_t_iterator_type), POINTER           :: iterator
      TYPE(dbcsr_t_type), POINTER                     :: tensor
      
      ALLOCATE(iterator)
      CALL c_f_pointer(c_tensor, tensor)
      
      CALL dbcsr_t_iterator_start(iterator, tensor)
      
      c_iterator = c_loc(iterator)
      
    END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of dbcsr_iterator_stop for tensors.
!> \param iterator ...
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_iterator_stop(c_iterator) BIND(C, name="c_dbcsr_t_iterator_stop")
   
      TYPE(c_ptr), INTENT(IN), VALUE :: c_iterator
      TYPE(dbcsr_t_iterator_type), POINTER :: iterator
      
      CALL c_f_pointer(c_iterator, iterator)
      
      CALL dbcsr_t_iterator_stop(iterator)
      
   END SUBROUTINE

! **************************************************************************************************
!> \brief iterate over nd blocks of an nd rank tensor, index only (blocks must be retrieved by calling
!>        dbcsr_t_get_block on tensor).
!> \param iterator ...
!> \param ind_nd nd index of block
!> \param blk is this needed?
!> \param blk_p is this needed?
!> \param blk_size blk size in each dimension
!> \param blk_offset blk offset in each dimension
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_iterator_next_block(c_iterator, iterator_size, c_ind_nd, c_blk, c_blk_p, c_blk_size, c_blk_offset) &
   BIND(C, name="c_dbcsr_t_iterator_next_block")
   
      TYPE(c_ptr), INTENT(IN), VALUE :: c_iterator
      INTEGER(kind=c_int), INTENT(IN), VALUE :: iterator_size
      INTEGER(kind=c_int), DIMENSION(iterator_size), TARGET, &
        INTENT(OUT) :: c_ind_nd
      INTEGER(kind=c_int), INTENT(OUT) :: c_blk
      INTEGER(kind=c_int), INTENT(OUT), OPTIONAL :: c_blk_p
      INTEGER(kind=c_int), DIMENSION(iterator_size), INTENT(OUT), &
        TARGET, OPTIONAL :: c_blk_size, c_blk_offset
      
      TYPE(dbcsr_t_iterator_type), POINTER     :: iterator
      INTEGER, DIMENSION(:), POINTER           :: ind_nd
      INTEGER, ALLOCATABLE                 :: blk_p
      INTEGER, DIMENSION(:), POINTER       :: blk_size => null()
      INTEGER, DIMENSION(:), POINTER       :: blk_offset => null()
         
         
      CALL c_f_pointer(c_iterator, iterator)
      ind_nd => c_ind_nd   
         
      IF (PRESENT(c_blk_p)) blk_p = c_blk_p
      IF (PRESENT(c_blk_size)) blk_size => c_blk_size    
      IF (PRESENT(c_blk_offset)) blk_offset => c_blk_offset

      CALL dbcsr_t_iterator_next_block(iterator, ind_nd, c_blk, blk_p, blk_size, blk_offset)
      
      c_ind_nd = c_ind_nd - 1
     
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of dbcsr_iterator_blocks_left for tensors.
!> \param iterator ...
!> \return ...
! **************************************************************************************************
   FUNCTION c_dbcsr_t_iterator_blocks_left(c_iterator) &
   BIND(C, name="c_dbcsr_t_iterator_blocks_left")
    
      TYPE(c_ptr), INTENT(IN), VALUE :: c_iterator
      TYPE(dbcsr_t_iterator_type), POINTER :: iterator
      LOGICAL(kind=c_bool) :: c_dbcsr_t_iterator_blocks_left
      
      CALL c_f_pointer(c_iterator, iterator);
      
      c_dbcsr_t_iterator_blocks_left = dbcsr_t_iterator_blocks_left(iterator)
     
   END FUNCTION
   
! **************************************************************************************************
!> \brief Split tensor blocks into smaller blocks of maximum size PRODUCT(block_sizes).
!> \param tensor_in Input tensor
!> \param tensor_out Output tensor (split blocks)
!> \param block_sizes block sizes for each of the tensor dimensions
!> \param nodata don't copy data from tensor_in to tensor_out
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_split_blocks(c_tensor_in, tensor_dim, c_tensor_out, c_block_sizes, c_nodata) &
   BIND(C, name="c_dbcsr_t_split_blocks")
   
      TYPE(c_ptr), INTENT(IN), VALUE :: c_tensor_in
      TYPE(c_ptr), INTENT(OUT) :: c_tensor_out
      INTEGER(kind=c_int), INTENT(IN), VALUE :: tensor_dim
      INTEGER(kind=c_int), INTENT(IN), DIMENSION(tensor_dim), & 
      TARGET :: c_block_sizes
      LOGICAL(kind=c_bool), OPTIONAL :: c_nodata
   
      TYPE(dbcsr_t_type), POINTER               :: tensor_in
      TYPE(dbcsr_t_type), POINTER                 :: tensor_out
      INTEGER, DIMENSION(:), POINTER :: block_sizes
      LOGICAL, ALLOCATABLE                  :: nodata
      
      CALL c_f_pointer(c_tensor_in, tensor_in)
      ALLOCATE(tensor_out)
      block_sizes => c_block_sizes
      IF (PRESENT(c_nodata)) nodata = c_nodata
      
      CALL dbcsr_t_split_blocks(tensor_in, tensor_out, block_sizes, nodata)
      
      c_tensor_out = c_loc(tensor_out)

   END SUBROUTINE

#:for dparam, dtype, dsuffix, dbase, dctype in cf_dtype_float_list
! **************************************************************************************************
!> \brief As dbcsr_filter
!> \param tensor ...
!> \param eps ...
!> \param method ...
!> \param use_absolute
! **************************************************************************************************
   SUBROUTINE c_dbcsr_t_filter_${dsuffix}$ (c_tensor, c_eps, c_method, c_use_absolute) &
   BIND(C, name="c_dbcsr_t_filter_${dsuffix}$")
       
      TYPE(c_ptr), INTENT(IN), VALUE :: c_tensor
      ${dbase}$(kind=${dctype}$), INTENT(IN), VALUE :: c_eps
      INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: c_method
      LOGICAL(kind=c_bool), INTENT(IN), OPTIONAL ::  c_use_absolute
      
      TYPE(dbcsr_t_type), POINTER    :: tensor
      INTEGER, ALLOCATABLE        :: method
      LOGICAL, ALLOCATABLE        :: use_absolute

      CALL c_f_pointer(c_tensor, tensor)
      IF (PRESENT(c_method)) method = c_method
      IF (PRESENT(c_use_absolute)) use_absolute = c_use_absolute

      CALL c_dbcsr_t_filter(tensor, c_eps, method, use_absolute)

   END SUBROUTINE
#:endfor
   
END MODULE
