subroutine write_bsplines_one_e_integrals
 implicit none
!provide block_t1_ij  block_t2_ij

 provide bsp_t_full
 provide bsp_vne_full
 provide bsp_s_full


 call ezfio_set_ao_one_e_ints_ao_integrals_overlap(bsp_s_full)
 call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(bsp_t_full)
 call ezfio_set_ao_one_e_ints_ao_integrals_e_n(bsp_vne_full)
 
 call ezfio_set_ao_one_e_ints_io_ao_integrals_e_n('Read')
 call ezfio_set_ao_one_e_ints_io_ao_integrals_kinetic('Read')
 call ezfio_set_ao_one_e_ints_io_ao_integrals_overlap('Read')

end
