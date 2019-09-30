program pouet
 implicit none
 integer :: i,j,k,l 
 double precision, allocatable :: s_mat_tmp(:,:)
 allocate(s_mat_tmp(ao_num,ao_num))
 do i = 1, ao_num
  do j = 1, ao_num
   ! compute the <i|j> = \sum_{kl} c_ik c_jl <k|l>
   double precision :: accu
   accu = 0.d0
   do k = 1, ao_num
    do l = 1, ao_num
     accu += ao_ortho_canonical_coef(l,i) * ao_ortho_canonical_coef(k,j) * ao_overlap(l,k)
    enddo
   enddo
   s_mat_tmp(j,i) = accu
  enddo 
 enddo
  print*,''
  print*,''
  print*,''
  do i = 1, ao_num
   write(*,'(100(F10.5,X))')s_mat_tmp(i,:)
  enddo

 double precision, allocatable :: H(:,:),eigvectors(:,:),eigvalues(:)
 allocate( H(ao_num,ao_num), eigvectors(ao_num,ao_num), eigvalues(ao_num) )
 H = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     H(j,i) += ( ao_kinetic_integrals(k,l) + ao_integrals_n_e(k,l) )  * ao_ortho_canonical_coef(l,i) * ao_ortho_canonical_coef(k,j) 
    enddo
   enddo
  enddo
  print*,'<i|H|i> = ',H(i,i)
 enddo
 call lapack_diagd(eigvalues,eigvectors,H,ao_num,ao_num)
 do i = 1, ao_num
 print*,'eigvalues = ',i,eigvalues(i) 
 enddo
 
 print*,''
 print*, 'two-electron integrals...'
 print*,''

! print*,'V_abcd',1,1,1,1,bsp_v_abcd(1,1,1,1)

end
