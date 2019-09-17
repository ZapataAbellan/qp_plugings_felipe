 BEGIN_PROVIDER [double precision, kin_int , (dim1,dim2)]
&BEGIN_PROVIDER [double precision, pot_int , (dim1,dim2)]
 implicit none
 integer :: i,j
 double precision :: tmp
 do i = 1, dim1
  do j = 1, dim2
    call pouet_2(i+j,tmp)
    kin_int(j,i) = tmp 
    pot_int(j,i) = tmp * 12.d0 
  enddo
 enddo

END_PROVIDER 
