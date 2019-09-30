 BEGIN_PROVIDER [double precision,  bsp_coef_norm  , (bsp_number)]

 implicit none

 !integers for loops...
 integer::i

 do i=1,bsp_number
    bsp_coef_norm(i)=1.d0/dsqrt(o_int(i,i))
 enddo

 END_PROVIDER 
