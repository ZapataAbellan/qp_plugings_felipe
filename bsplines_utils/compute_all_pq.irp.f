
subroutine give_all_bsplines_jl(ao_p,ao_q,pq_integrals)
           
BEGIN_DOC
!We enter ao_p and ao_q which are the index of two atomic orbitals and we get out 
!full_integral(t,u) = (ao_p(1) a_q(1) | ao_t(2) ao_u(2) ) within chemist notation
END_DOC
!
implicit none
!
!atomic orbitals (ao)
integer,intent(in) :: ao_p, ao_q
integer :: ao_t, ao_u
!full integral
double precision,intent(out) :: pq_integrals(ao_dim,ao_dim)
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
pq_integrals(:,:) = 0.d0
!
!============================
!ao_p (electron 1)
rp = ao_bspline_index_r(ao_p) 
lp = ao_bspline_index_l(ao_p)
mp = ao_bspline_index_m(ao_p)
!ao_q (electron 1)
rq = ao_bspline_index_r(ao_q)
lq = ao_bspline_index_l(ao_q)
mq = ao_bspline_index_m(ao_q)
!============================
!
!Loop over the ao_t and ao_u orbitals
!
!ao_t (electron 2)
do lt=0,bsp_lmax
 do mt=-lt,lt
  !ao_u (electron 2)
  do lu=0,bsp_lmax
   do mu=-lu,lu
   !
   !Multipolar expansion limits...
   kmin = max(abs(lp-lq),abs(lt-lu))
   kmax = min(abs(lp+lq),abs(lt+lu))
   !Loop multipolar expansion...
   do k=kmin,kmax
    !angular coefficient for a given k... 
    angular = 0.d0
    do mk=-k,k
     angular +=  gaunt(lp,mp,k,mk,lq,mq) * gaunt(lt,mt,k,mk,lu,mu) 
    end do
    !Loop over the radial grid...
    !Attention: 
    !R^k is stored in the chemist notation...
    !
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
         !
         !=======================================================
         if (((i +iv-1).gt.1).and.((i +iv-1).lt.bsp_number)) then 
         if (((j +jv-1).gt.1).and.((j +jv-1).lt.bsp_number)) then
         if (((ip+iv-1).gt.1).and.((ip+iv-1).lt.bsp_number)) then
         if (((jp+jv-1).gt.1).and.((jp+jv-1).lt.bsp_number)) then
         !=======================================================
         if ((rp.eq.((i+iv-1)-1)).and.(rq.eq.((ip+iv-1)-1))) then
         !======================================================= 
         !
         ao_t = ao_to_bspline(mt,lt,(j +jv-1)-1) !(electron 2)
         ao_u = ao_to_bspline(mu,lu,(jp+jv-1)-1) !(electron 2)

         pq_integrals(ao_t,ao_u) += c * angular 

         end if 
          
         end if
         end if 
         end if 
         end if

         end do
        end do
       end do !loop : radial intervals r2
      end do 
     end do
    end do !loop : radial intervals r1
    !
   end do !loop: k
   !
   end do !loop:mu
  end do !loop:lu
  !
 end do !loop:mt
end do !loop:lt


end
