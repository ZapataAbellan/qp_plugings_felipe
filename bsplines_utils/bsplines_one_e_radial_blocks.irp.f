 BEGIN_PROVIDER [double precision,  radial_block_bsp_t1  , (bsp_dim,bsp_dim)]
&BEGIN_PROVIDER [double precision,  radial_block_bsp_t2  , (bsp_dim,bsp_dim)]
&BEGIN_PROVIDER [double precision,  radial_block_bsp_s   , (bsp_dim,bsp_dim)]
&BEGIN_PROVIDER [double precision,  radial_block_bsp_v   , (bsp_dim,bsp_dim)]

 implicit none

 !integers for loops...
 integer:: i,j
 !
 !computes the radial blocks of the atomic orbital matrices...
 !
 do j=2,bsp_number-1
    do i=2,bsp_number-1
       radial_block_bsp_t1(i-1,j-1)=-0.5d0*d2r_int(i,j)*bsp_coef_norm(i)*bsp_coef_norm(j)
       radial_block_bsp_t2(i-1,j-1)=r2_int(i,j)*bsp_coef_norm(i)*bsp_coef_norm(j)
       radial_block_bsp_v(i-1,j-1)=-dfloat(bsp_atomic_charge)*r1_int(i,j)*bsp_coef_norm(i)*bsp_coef_norm(j)
       radial_block_bsp_s(i-1,j-1)=o_int(i,j)*bsp_coef_norm(i)*bsp_coef_norm(j)
    enddo
 enddo
  
 END_PROVIDER
