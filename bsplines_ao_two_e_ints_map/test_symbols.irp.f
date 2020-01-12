program test
BEGIN_DOC
!
END_DOC
implicit none

integer :: l, m, lp, mp, k
double precision :: c_k_1, c_k_2
double precision :: c_k_full_1, c_k_full_2

write(*,*)bsp_lmax,bsp_dlmax

!l=4
!m=1
!
!lp=2
!mp=2
!
!k=2
!
!
!c_k_1      = (-1.d0)**(l +m )*bsp_reduced_ck(l,lp,k)*bsp_racah_v(l,lp,k,-m ,mp,m-mp)
!c_k_full_1 = bsp_full_ck(l,m,lp,mp,k)
!
!c_k_2      = (-1.d0)**(lp+mp)*bsp_reduced_ck(lp,l,k)*bsp_racah_v(lp,l,k,-mp,m ,mp-m)
!c_k_full_2 = bsp_full_ck(lp,mp,l,m,k)
!
!write(*,*) c_k_1     , c_k_2     , c_k_1      - c_k_2      *(-1.d0)**(m-mp)
!write(*,*) c_k_full_1, c_k_full_2, c_k_full_1 - c_k_full_2 *(-1.d0)**(m-mp)  
write(99,'("  l  lp m  mp k  ck ")')

integer :: kmin, kmax

do l=0,bsp_lmax
 do lp=0,bsp_lmax
  do m=-l,l
   do mp=-lp,lp
    kmin=max(abs(m-mp),abs(l-lp))
    kmax=l+lp 
    do k=kmin,kmax
     if (mod((l+lp+k),2).eq.0) then
      write(99,'(5(1x,I2),1x,ES20.12E02)')l,lp,m,mp,k,bsp_full_ck(l,m,lp,mp,k)
     endif
    enddo
   enddo
  enddo
 enddo
enddo

stop
end program
