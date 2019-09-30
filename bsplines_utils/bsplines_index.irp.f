 BEGIN_PROVIDER [integer,  ao_bspline_index_r , (ao_dim)]
&BEGIN_PROVIDER [integer,  ao_bspline_index_l , (ao_dim)]
&BEGIN_PROVIDER [integer,  ao_bspline_index_m , (ao_dim)]
&BEGIN_PROVIDER [integer,  ao_to_bspline , (-bsp_lmax:bsp_lmax,0:bsp_lmax,bsp_dim)]

 implicit none
 
 integer:: ao_I
 integer:: i,l,m
 ao_I=0

 do i=1,bsp_dim
  do l=0,bsp_lmax
   do m=-l,l
    ao_I +=1 
    ao_bspline_index_r(ao_I) = i
    ao_bspline_index_l(ao_I) = l
    ao_bspline_index_m(ao_I) = m
    ao_to_bspline(m,l,i) = ao_I
   enddo
  enddo
 enddo 
 END_PROVIDER
