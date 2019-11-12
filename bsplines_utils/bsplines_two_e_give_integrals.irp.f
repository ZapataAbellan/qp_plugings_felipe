program test
implicit none
integer :: ao_p,ao_q,ao_t,ao_u
double precision, allocatable :: full_integral(:,:)

allocate(full_integral(ao_dim,ao_dim))

write(*,*)'two-electron integrals...'


end


subroutine give_two_e_integrals_qp(ao_p,ao_q,full_integral)
BEGIN_DOC
!We enter ao_q and ao_q which are the index of two atomic orbitals and we get out with
!full_integral(t,u) = (ao_p(1) a_q(1) | ao_t(2) ao_u(2) ) within the chemist notation
END_DOC
!
implicit none
!
!atomic orbitals (ao)
integer,intent(in) :: ao_p, ao_q
integer :: ao_t, ao_u
!full integral
double precision,intent(out) :: full_integral(ao_dim,ao_dim)
!
!ao radial 
integer :: rp, rq, rt, ru
!ao angular momentum
integer :: lp, lq, lt, lu
!ao magnetic momentum
integer :: mp, mq, mt, mu
!
!multipolar expansion 
integer :: k, kmax, kmin, mk
!
!Gaunt coefficient  <lm|C_mk^k|lpmp>
double precision :: gaunt
double precision :: angular
!
!Slater radial integral coefficient R^k(p,q,t,u)
double precision :: c
!
!integres for loops in the assembling algorithm...
integer ::i, j, ip, jp, ii, jj, iv, jv
!
!initialization
full_integral(:,:) = 0.d0
!
!ao_p
rp = ao_bspline_index_r(ao_p) 
lp = ao_bspline_index_l(ao_p)
mp = ao_bspline_index_m(ao_p)
!ao_q
rq = ao_bspline_index_r(ao_q)
lq = ao_bspline_index_l(ao_q)
mq = ao_bspline_index_m(ao_q)
!
!Loop over the radial grid...
do jv=1,bsp_nv
 jj=0
 do j=1,bsp_order
  do jp=1,bsp_order
   jj = jj + 1
   do iv=1,bsp_nv
    ii=0
    do i=1,bsp_order
     do ip=1,bsp_order
      ii = ii + 1

      if( iv < jv ) then
       c = roff1_kmax(ii,iv,k)*roff2_kmax(jj,jv,k) 
      else if( iv > jv ) then
       c = roff2_kmax(ii,iv,k)*roff1_kmax(jj,jv,k)
      else
       c = rdiag_kmax(ii,jj,iv,k)
      end if

      print*, 'rp,rq',((i +iv-1)-1),((ip+iv-1)-1)

     end do
    end do
   end do !loop : radial intervals r2
  end do 
 end do
end do !loop : radial intervals r1
       !


end
