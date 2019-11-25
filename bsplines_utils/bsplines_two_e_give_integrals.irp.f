program test
implicit none
integer :: ao_p,ao_q,ao_t,ao_u
double precision, allocatable :: pq_integrals(:,:)
double precision :: tmp

allocate(pq_integrals(ao_dim,ao_dim))

write(*,*)'two-electron integrals...'

!ao_p = 1
!ao_q = 1

write(*,*)'ao_dim',ao_dim
!write(*,*)'ao_p, ao_q',ao_p,ao_q

!call give_two_e_integrals_pq(ao_p,ao_q,pq_integrals)

do ao_p=1, ao_dim
 do ao_q=1, ao_dim
  
  call give_two_e_integrals_pq(ao_p,ao_q,pq_integrals)
  
  do ao_t=1, ao_dim
   do ao_u=1, ao_dim
    
    tmp = pq_integrals(ao_t,ao_u)

    if ( dabs(tmp) .gt. 1.d-10 ) then
     if ( dabs(bsp_vee_full(ao_p,ao_q,ao_t,ao_u)-tmp)/dabs(tmp) .gt. 1.d-5) then

     write(55,'(4(1x,I6),2(1x,ES20.12E02))')ao_p,ao_q,ao_t,ao_u,pq_integrals(ao_t,ao_u),bsp_vee_full(ao_p,ao_q,ao_t,ao_u)

     end if
    end if

   end do 
  end do
 end do
end do

end


