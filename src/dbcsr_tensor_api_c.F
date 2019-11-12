#:include "tensors/dbcsr_tensor.fypp"
#:set ndims = range(1,maxrank+1)

MODULE dbcsr_tensor_api_c

   USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_loc, c_ptr, c_double, C_NULL_CHAR, &
      c_f_pointer, c_int, c_char, c_null_ptr, c_bool
   USE dbcsr_tensor_api
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
      
      PRINT*, "SIZE: ", map1_2d_size

      IF (PRESENT(c_map1_2d) .AND. PRESENT(c_map2_2d)) THEN
        map1_2d => c_map1_2d
        map2_2d => c_map2_2d
        PRINT*, "PRESENT"
      ELSE 
        map1_2d => null()
        map2_2d => null()
        PRINT*, "NOT PRESENT"
      ENDIF
    
      ALLOCATE (pgrid)
      
      CALL dbcsr_t_pgrid_create(fcomm, dims, pgrid, map1_2d, map2_2d, &
      nsplit, dimsplit)
      
      c_pgrid = c_loc(pgrid)
      
      PRINT*, "PGRID", pgrid%mp_comm_2d
      
   END SUBROUTINE

!#:for ndim in ndims

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
      
      map1_2d => c_map1_2d
      map2_2d => c_map2_2d
      
#:for ndim in ndims
      IF (PRESENT(c_nd_dist_${ndim}$)) THEN
        nd_dist_${ndim}$ => c_nd_dist_${ndim}$
      ELSE 
        nd_dist_${ndim}$ => null()
      ENDIF
      
#:endfor
    
      IF (PRESENT(c_own_comm)) THEN
        own_comm = c_own_comm
        PRINT*, "Bool is there."
      ENDIF

      PRINT*, "HERE"

      ALLOCATE(dist)
      
      CALL c_f_pointer(c_pgrid, pgrid)
      
      PRINT*, "ENTERING"

      CALL dbcsr_t_distribution_new(dist, pgrid, map1_2d, map2_2d, &
      ${varlist("nd_dist")}$, own_comm)
      
      PRINT*, "EXITING"

      c_dist = c_loc(dist)

END SUBROUTINE

!#:endfor

SUBROUTINE c_dbcsr_t_create_new(c_tensor, c_name, c_dist, c_map1_2d, map1_2d_size, &
        c_map2_2d, map2_2d_size, data_type, ${c_varlist_and_size("c_blk_size")}$) &
        bind(C, name="c_dbcsr_t_create_new")
                                 
      TYPE(c_ptr), INTENT(OUT)                          :: c_tensor                           
      TYPE(c_ptr), INTENT(IN), value                    :: c_dist, c_name
        
      INTEGER(kind=c_int), INTENT(in), value             :: map1_2d_size
      INTEGER(kind=c_int), INTENT(in), TARGET            :: c_map1_2d(map1_2d_size)
      INTEGER(kind=c_int), INTENT(in), value             :: map2_2d_size
      INTEGER(kind=c_int), INTENT(in), TARGET            :: c_map2_2d(map2_2d_size)
      
      INTEGER(kind=c_int), INTENT(IN), OPTIONAL          :: data_type 
      
      INTEGER(kind=c_int), INTENT(in), value             :: ${varlist("c_blk_size", suffix="_size")}$
      INTEGER(kind=c_int), INTENT(IN), TARGET, OPTIONAL  :: ${c_arrlist_and_size("c_blk_size")}$
    
      CHARACTER(:, kind=c_char), ALLOCATABLE             :: fname
      
      INTEGER, POINTER                                  :: map1_2d(:), map2_2d(:)
      INTEGER, POINTER                                  :: ${pointer_list("blk_size")}$
      TYPE(dbcsr_t_type), POINTER                        :: tensor
      TYPE(dbcsr_t_distribution_type), POINTER           :: dist
    
      map1_2d => c_map1_2d
      map2_2d => c_map2_2d
      
#:for ndim in ndims
      IF (PRESENT(c_blk_size_${ndim}$)) THEN
        blk_size_${ndim}$ => c_blk_size_${ndim}$
      ELSE 
        blk_size_${ndim}$ => null()
      ENDIF
      
#:endfor

      ALLOCATE(tensor)
      
      CALL c_f_pointer(c_dist, dist)
      CALL c_f_string(c_name, fname)
      
      CALL dbcsr_t_create(tensor, fname, dist, map1_2d, map2_2d, data_type, &
                                 ${varlist("blk_size")}$)
      
      c_tensor = c_loc(tensor)
      
END SUBROUTINE
     
     


END MODULE
