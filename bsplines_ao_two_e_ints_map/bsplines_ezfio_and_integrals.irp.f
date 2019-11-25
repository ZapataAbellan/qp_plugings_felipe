program bsplines_ezfio

 implicit none
 integer :: ao_p,ao_q,ao_t,ao_u
 double precision, allocatable :: pq_integrals(:,:)
 double precision :: tmp,get_ao_two_e_integral_bsplines,tmp2
 
 allocate(pq_integrals(ao_dim,ao_dim))
 
 write(*,*)'two-electron integrals...'
 call change_basis 
 provide ao_two_e_integrals_bsplines_in_map 
 !ao_p = 1
 !ao_q = 1
 
 write(*,*)'ao_dim',ao_dim
 
 do ao_p=1, ao_dim
  do ao_q=1, ao_dim
   
   call give_all_bsplines_jl(ao_p,ao_q,pq_integrals)
   do ao_t=1, ao_dim
    do ao_u=1, ao_dim
     tmp  = bsp_vee_full(ao_t,ao_u,ao_q,ao_p)
     tmp2 = bsp_vee_full(ao_u,ao_t,ao_q,ao_p)
 
     if ( dabs(tmp2) .gt. 1.d-10 ) then
      if ( dabs(tmp2-tmp)/dabs(tmp) .gt. 1.d-5) then
 
      write(55,'(4(1x,I6),2(1x,ES20.12E02))')ao_p,ao_q,ao_t,ao_u,tmp2,tmp
 
      end if
     end if
 
    end do 
   end do
  end do
 end do

end


!     tmp = get_ao_two_e_integral_bsplines(ao_p,ao_t,ao_q,ao_u,ao_integrals_bsplines_map) 
!     tmp = pq_integrals(ao_t,ao_u)
