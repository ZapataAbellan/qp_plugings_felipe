! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/oem/Documentos/UPMC/LCT/LCT_FELIPE/PROJECTS/quantum_package/qp2/src/bsplines_utils/EZFIO.cfg


BEGIN_PROVIDER [ integer, bsp_order  ]
  implicit none
  BEGIN_DOC
! B-spline order
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_bsplines_utils_bsp_order(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: bsp_order ] <<<<< ..'
      call ezfio_get_bsplines_utils_bsp_order(bsp_order)
    else
      print *, 'bsplines_utils/bsp_order not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( bsp_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read bsp_order with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ integer, bsp_lmax  ]
  implicit none
  BEGIN_DOC
! maximum angular momentum
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_bsplines_utils_bsp_lmax(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: bsp_lmax ] <<<<< ..'
      call ezfio_get_bsplines_utils_bsp_lmax(bsp_lmax)
    else
      print *, 'bsplines_utils/bsp_lmax not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( bsp_lmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read bsp_lmax with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ integer, bsp_glp  ]
  implicit none
  BEGIN_DOC
! Gauss-Legendre points
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_bsplines_utils_bsp_glp(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: bsp_glp ] <<<<< ..'
      call ezfio_get_bsplines_utils_bsp_glp(bsp_glp)
    else
      print *, 'bsplines_utils/bsp_glp not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( bsp_glp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read bsp_glp with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ double precision, bsp_box_size  ]
  implicit none
  BEGIN_DOC
! size of the simulation box
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_bsplines_utils_bsp_box_size(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: bsp_box_size ] <<<<< ..'
      call ezfio_get_bsplines_utils_bsp_box_size(bsp_box_size)
    else
      print *, 'bsplines_utils/bsp_box_size not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( bsp_box_size, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read bsp_box_size with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ character*(32), bsp_grid_type  ]
  implicit none
  BEGIN_DOC
! type of B-spline grid, can be [ linear | exponential | fischer | parabolic ]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_bsplines_utils_bsp_grid_type(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: bsp_grid_type ] <<<<< ..'
      call ezfio_get_bsplines_utils_bsp_grid_type(bsp_grid_type)
    else
      print *, 'bsplines_utils/bsp_grid_type not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( bsp_grid_type, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read bsp_grid_type with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ integer, bsp_number  ]
  implicit none
  BEGIN_DOC
! Total number of B-splines
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_bsplines_utils_bsp_number(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: bsp_number ] <<<<< ..'
      call ezfio_get_bsplines_utils_bsp_number(bsp_number)
    else
      print *, 'bsplines_utils/bsp_number not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( bsp_number, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read bsp_number with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
