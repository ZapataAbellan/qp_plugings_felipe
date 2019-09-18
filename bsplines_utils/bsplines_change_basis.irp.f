program change_basis
 implicit none
 BEGIN_DOC
! program that trick the AO basis in order to have:
!
! a new set of AO basis with only one primitive per AO,
!
! each AO is the same S function of exponent pi/2
!
! such that it is normalized to unity
 END_DOC

 integer :: ao_num_new
 PROVIDE ezfio_filename
 ao_num_new = bsp_dim
 ! First you change the size of the AO basis
 call routine_ao_num(ao_num_new)
 ! Then the number of primitives per AO
 call routine_ao_prim_num(ao_num_new)
 ! Then that only S functions
 call routine_ao_power(ao_num_new)
 ! Then what atom is attached to what AO, by default all AOs are attached to atom 1
 call routine_ao_nucl(ao_num_new)
 ! Then the coeficients of the AO are set to 1.d0 as they are composed of only 1 primitive
 call routine_ao_coef(ao_num_new)
 ! Then the exponent of the gaussians, set to pi/2 such that they are normalized to unity
 call routine_ao_expo(ao_num_new)

 !!!!!!!!!!!!!!!!!! INTEGRALS PART !!!!!!!!!!!!!!!
 ! here you write the overlap, kinetic and e-n potential integrals on disk
 ! and you specify to read all these integrals
 call write_bsplines_one_e_integrals

end

 subroutine routine_ao_num(ao_num_new)
  implicit none
  integer, intent(in) :: ao_num_new
  call ezfio_set_ao_basis_ao_num(ao_num_new)
 end

 subroutine routine_ao_power(ao_num_new)
  ! Then you set that there are only S functions
  implicit none
  integer, intent(in) :: ao_num_new
  integer, allocatable :: ao_power_new(:,:)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
  allocate(ao_power_new(ao_num_new,3))
  ao_power_new = 0
  call ezfio_set_ao_basis_ao_power(ao_power_new)

 end

 subroutine routine_ao_prim_num(ao_num_new)
  ! Then you set that there are only one primitive per AO
  implicit none
  integer, intent(in) :: ao_num_new
  integer, allocatable :: ao_prim_num_new(:)
  allocate(ao_prim_num_new(ao_num_new))
  ao_prim_num_new = 1
  call ezfio_set_ao_basis_ao_prim_num(ao_prim_num_new)
 end


 subroutine routine_ao_nucl(ao_num_new)
  ! Then you change the ao_nucl array
  implicit none
  integer, intent(in) :: ao_num_new
  integer, allocatable :: ao_nucl_new(:)
  allocate(ao_nucl_new(ao_num_new))
  ao_nucl_new = 1
  call ezfio_set_ao_basis_ao_nucl(ao_nucl_new)
 end

 subroutine routine_ao_coef(ao_num_new)
  ! Then you set the ao_coef to be 1.d0
  implicit none
  integer, intent(in) :: ao_num_new
  integer :: ao_prim_num_max_new
  ao_prim_num_max_new = 1
  double precision, allocatable :: ao_coef_new(:,:)
  allocate(ao_coef_new(ao_num_new,ao_prim_num_max_new))
  ao_coef_new = 1.d0
  call ezfio_set_ao_basis_ao_coef(ao_coef_new)
 end

 subroutine routine_ao_expo(ao_num_new)
 ! Then you set the ao_expo such that each AO is normalized to unity, which is a gaussian of exponent pi/2
  implicit none
  integer, intent(in) :: ao_num_new
  integer :: ao_prim_num_max_new
  double precision, allocatable :: ao_expo_new(:,:)
  ao_prim_num_max_new = 1
  allocate(ao_expo_new(ao_num_new,ao_prim_num_max_new))
  ao_expo_new(:,:) = dacos(-1.d0) * 0.5d0
  call ezfio_set_ao_basis_ao_expo(ao_expo_new)
 end

subroutine write_bsplines_one_e_integrals
 implicit none
 call ezfio_set_ao_one_e_ints_ao_integrals_overlap(bsp_s_ij)
 call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(bsp_t_ij)
 call ezfio_set_ao_one_e_ints_ao_integrals_e_n(bsp_v_ij)

 integer :: j,i
 print*,''
 print*,''
 print*,'AO overlap matrix'
 print*,''
 do i = 1, bsp_dim
  write(*,'(100(F18.14,X))')bsp_s_ij(i,:)
  do j = 1, bsp_dim
   if(dabs(bsp_s_ij(i,j)).gt.1.d+2)then
    print*,''
    print*,''
    print*,i,j,bsp_s_ij(i,j)
    print*,''
    print*,''
   endif
  enddo
 enddo
 print*,'AO overlap matrix'
 print*,'AO overlap matrix'
 
 call ezfio_set_ao_one_e_ints_io_ao_integrals_e_n('Read')
 call ezfio_set_ao_one_e_ints_io_ao_integrals_kinetic('Read')
 call ezfio_set_ao_one_e_ints_io_ao_integrals_overlap('Read')

end
