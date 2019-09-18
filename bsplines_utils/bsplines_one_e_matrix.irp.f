 BEGIN_PROVIDER [double precision,  block_t1_ij  , (bsp_dim,bsp_dim)]
&BEGIN_PROVIDER [double precision,  block_t2_ij  , (bsp_dim,bsp_dim)]
&BEGIN_PROVIDER [double precision,  block_bsp_s_ij   , (bsp_dim,bsp_dim)]
&BEGIN_PROVIDER [double precision,  block_bsp_v_ij   , (bsp_dim,bsp_dim)]

 implicit none

 !integers for loops...
 integer:: i,j

 do j=2,bsp_number-1
    do i=2,bsp_number-1
       block_t1_ij(i-1,j-1)=-0.5d0*d2r_int(i,j)*s_coef_norm(i)*s_coef_norm(j)
       block_t2_ij(i-1,j-1)=r2_int(i,j)*s_coef_norm(i)*s_coef_norm(j)
       block_bsp_v_ij(i-1,j-1)=-dfloat(bsp_atomic_charge)*r1_int(i,j)*s_coef_norm(i)*s_coef_norm(j)
       block_bsp_s_ij(i-1,j-1)=o_int(i,j)*s_coef_norm(i)*s_coef_norm(j)
    enddo
 enddo
  
 END_PROVIDER 

BEGIN_PROVIDER [integer, bsp_total_dim]
 implicit none
 BEGIN_DOC
 ! final dimention of bspline basis
 END_DOC
 bsp_total_dim = bsp_dim*(1+bsp_lmax)
END_PROVIDER

 BEGIN_PROVIDER [double precision,  bsp_t_ij  , (bsp_total_dim,bsp_total_dim)]
 implicit none
 !angular momentum...
 integer:: l
 !angular momentum block...
 integer:: b
 !integers for loops...
 integer:: i,j
 b = 0
  do l=0,bsp_lmax
     do j=1,bsp_dim
        do i=1,bsp_dim
           bsp_t_ij(i+b,j+b)=block_t1_ij(i,j)+(dfloat(l*(l+1))/2.d0)*block_t2_ij(i,j)
        enddo
     enddo
     b = b + bsp_dim
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [double precision,  bsp_s_ij  , (bsp_total_dim,bsp_total_dim)]
 implicit none
 !angular momentum...
 integer:: l
 !angular momentum block...
 integer:: b
 !integers for loops...
 integer:: i,j
 b = 0
  do l=0,bsp_lmax
     do j=1,bsp_dim
        do i=1,bsp_dim
           bsp_s_ij(i+b,j+b)=block_bsp_s_ij(i,j)
        enddo
     enddo
     b = b + bsp_dim
  enddo
 END_PROVIDER
 

 BEGIN_PROVIDER [double precision,  bsp_v_ij  , (bsp_total_dim,bsp_total_dim)]
 implicit none
 !angular momentum...
 integer:: l
 !angular momentum block...
 integer:: b
 !integers for loops...
 integer:: i,j
 b = 0
  do l=0,bsp_lmax
     do j=1,bsp_dim
        do i=1,bsp_dim
           bsp_v_ij(i+b,j+b)=block_bsp_v_ij(i,j)
        enddo
     enddo
     b = b + bsp_dim
  enddo
 END_PROVIDER
