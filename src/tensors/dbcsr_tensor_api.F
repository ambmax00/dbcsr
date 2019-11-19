!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   This is the start of a dbcsr_tensor_api, all publically needed functions
!>          are exported here. The others remain private to the library.
!>          Currently, this is the CP2K used set.
!>          Ultimately, a reduced subset and well defined api will remain,
!>          possibly grouped in to standard and expert api.
!>          Currently, this is work in progress.
!> \author  Alfio Lazzaro
! **************************************************************************************************
MODULE dbcsr_tensor_api
   USE dbcsr_tensor, ONLY: &
      dbcsr_t_contract, dbcsr_t_get_block, dbcsr_t_get_stored_coordinates, dbcsr_t_put_block, &
      dbcsr_t_reserve_blocks, dbcsr_t_copy_matrix_to_tensor, dbcsr_t_copy, &
      dbcsr_t_copy_tensor_to_matrix, dbcsr_t_need_contract
   USE dbcsr_tensor_block, ONLY: &
      dbcsr_t_iterator_blocks_left, dbcsr_t_iterator_next_block, dbcsr_t_iterator_start, &
      dbcsr_t_iterator_stop, dbcsr_t_iterator_type
   USE dbcsr_tensor_types, ONLY: &
      dbcsr_t_create, dbcsr_t_destroy, dbcsr_t_distribution_destroy, dbcsr_t_distribution_new, &
      dbcsr_t_distribution_type, dbcsr_t_nd_mp_comm, dbcsr_t_nd_mp_free, dbcsr_t_type, &
      dbcsr_t_pgrid_type, dbcsr_t_pgrid_create, dbcsr_t_pgrid_destroy, dbcsr_t_set, dbcsr_t_filter, &
      dbcsr_t_mp_environ_pgrid => mp_environ_pgrid, dbcsr_t_blk_sizes, dbcsr_t_get_info, &
      dbcsr_t_finalize, dbcsr_t_scale, dbcsr_t_get_nze, dbcsr_t_get_nze_total, &
      dbcsr_t_get_num_blocks, dbcsr_t_get_num_blocks_total, dbcsr_t_clear
   USE dbcsr_tensor_test, ONLY: &
      dbcsr_t_contract_test, dbcsr_t_checksum
   USE dbcsr_tensor_split, ONLY: &
      dbcsr_t_split_blocks
   USE dbcsr_tensor_index, ONLY: &
      dbcsr_t_get_mapping_info
   USE dbcsr_tensor_io, ONLY: &
      dbcsr_t_write_split_info, dbcsr_t_write_blocks

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor_api'

   PUBLIC :: dbcsr_t_contract, dbcsr_t_contract_test
   PUBLIC :: dbcsr_t_get_block ! done
   PUBLIC :: dbcsr_t_get_stored_coordinates ! done
   PUBLIC :: dbcsr_t_put_block ! done
   PUBLIC :: dbcsr_t_reserve_blocks ! done
   PUBLIC :: dbcsr_t_create !done
   PUBLIC :: dbcsr_t_destroy !done
   PUBLIC :: dbcsr_t_distribution_destroy ! done
   PUBLIC :: dbcsr_t_distribution_new ! done
   PUBLIC :: dbcsr_t_distribution_type ! skip
   PUBLIC :: dbcsr_t_nd_mp_comm !skip
   PUBLIC :: dbcsr_t_nd_mp_free !skip
   PUBLIC :: dbcsr_t_type !skip
   PUBLIC :: dbcsr_t_iterator_next_block !done
   PUBLIC :: dbcsr_t_iterator_blocks_left !done
   PUBLIC :: dbcsr_t_iterator_stop !done
   PUBLIC :: dbcsr_t_iterator_start !done
   PUBLIC :: dbcsr_t_iterator_type !skip
   PUBLIC :: dbcsr_t_split_blocks !done
   PUBLIC :: dbcsr_t_pgrid_type !skip
   PUBLIC :: dbcsr_t_pgrid_create ! done
   PUBLIC :: dbcsr_t_pgrid_destroy ! done
   PUBLIC :: dbcsr_t_set !done
   PUBLIC :: dbcsr_t_filter !done
   PUBLIC :: dbcsr_t_mp_environ_pgrid !skip
   PUBLIC :: dbcsr_t_copy_matrix_to_tensor !done
   PUBLIC :: dbcsr_t_blk_sizes ! done
   PUBLIC :: dbcsr_t_copy ! done
   PUBLIC :: dbcsr_t_copy_tensor_to_matrix !done
   PUBLIC :: dbcsr_t_get_info !done
   PUBLIC :: dbcsr_t_checksum ! skip
   PUBLIC :: dbcsr_t_finalize ! skip
   PUBLIC :: dbcsr_t_need_contract ! skip
   PUBLIC :: dbcsr_t_scale ! done
   PUBLIC :: dbcsr_t_get_num_blocks, dbcsr_t_get_num_blocks_total
   PUBLIC :: dbcsr_t_get_nze, dbcsr_t_get_nze_total
   PUBLIC :: dbcsr_t_clear
   PUBLIC :: dbcsr_t_get_mapping_info
   PUBLIC :: dbcsr_t_write_split_info
   PUBLIC :: dbcsr_t_write_blocks

END MODULE dbcsr_tensor_api
