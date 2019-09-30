 BEGIN_PROVIDER [double precision,  bsp_s_full  , (ao_dim,ao_dim)]
 implicit none
 !
 double precision:: tmp_s
 !
 !b-spline index
 integer:: i,j
 !ao index
 integer::ao_I,ao_J
 !
 !angular index
 integer:: l, m 
 !
 !inizialization
 bsp_s_full=0.d0
 !
 do i=1,bsp_dim
  !
  do j=1,bsp_dim
   !
   !radial contribution...
   tmp_s = radial_block_bsp_s(i,j) 
   !
   !angular contribution...
   do l=0,bsp_lmax
    do m=-l,l
     !
     !atomic orbital index
     ao_I = ao_to_bspline(m,l,i) 
     ao_J = ao_to_bspline(m,l,j)
   
     bsp_s_full(ao_I,ao_J) = tmp_s

    enddo
   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision,  bsp_vne_full  , (ao_dim,ao_dim)]
 implicit none
 !
 !temporal value
 double precision:: tmp_vne
 !
 !b-spline index
 integer:: i,j
 !ao index
 integer::ao_I,ao_J

 !
 !angular index
 integer:: l, m 
 !
 !inizialization
 bsp_vne_full=0.d0
 !
 do i=1,bsp_dim
  !
  do j=1,bsp_dim
   !
   !radial contribution...
   tmp_vne = radial_block_bsp_v(i,j) 
   !
   !angular contribution...
   do l=0,bsp_lmax
    do m=-l,l
     !
     !atomic orbital index
     ao_I = ao_to_bspline(m,l,i) 
     ao_J = ao_to_bspline(m,l,j)
   
     bsp_vne_full(ao_I,ao_J) = tmp_vne

    enddo
   enddo
  enddo
 enddo

 END_PROVIDER 


 BEGIN_PROVIDER [double precision,  bsp_t_full  , (ao_dim,ao_dim)]
 implicit none
 !
 !temporal value
 double precision:: tmp_t1, tmp_t2
 !
 !b-spline index
 integer:: i,j
 !
 !angular index
 integer:: l, m 
!
 integer::ao_I,ao_J
 !
 !angular coefficient
 double precision:: aa
 !
 !inizialization
 bsp_t_full=0.d0
 !
 do i=1,bsp_dim
  !
  do j=1,bsp_dim
   !
   !radial contribution...
   tmp_t1 = radial_block_bsp_t1(i,j)
   tmp_t2 = radial_block_bsp_t2(i,j)
   !
   !angular contribution...
   do l=0,bsp_lmax
    do m=-l,l
     !
     !atomic orbital index
     ao_I = ao_to_bspline(m,l,i) 
     ao_J = ao_to_bspline(m,l,j)
     
     !anuglar coefficient
     aa = 0.5d0*dfloat(l*(l+1))
     bsp_t_full(ao_I,ao_J) = tmp_t1 + aa*tmp_t2

    enddo
   enddo
  enddo
 enddo

 END_PROVIDER 

