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
 !!!!!!!!!!!!!!!!!! INTEGRALS PART !!!!!!!!!!!!!!!
 ! here you write the overlap, kinetic and e-n potential integrals on disk
 ! and you specify to read all these integrals
 print*,''
 print*,''
 print*,''
 print*,'ao_num = ',ao_num
 print*,''
 print*,''
 print*,''
 call write_bsplines_one_e_integrals

end

subroutine write_bsplines_one_e_integrals
 implicit none
!provide block_t1_ij  block_t2_ij

 provide bsp_t_ij
 provide bsp_v_ij
 provide bsp_s_ij


 call ezfio_set_ao_one_e_ints_ao_integrals_overlap(bsp_s_ij)
 call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(bsp_t_ij)
 call ezfio_set_ao_one_e_ints_ao_integrals_e_n(bsp_v_ij)
 
 call ezfio_set_ao_one_e_ints_io_ao_integrals_e_n('Read')
 call ezfio_set_ao_one_e_ints_io_ao_integrals_kinetic('Read')
 call ezfio_set_ao_one_e_ints_io_ao_integrals_overlap('Read')

end
