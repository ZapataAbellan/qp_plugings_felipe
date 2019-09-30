!       program angular
!       implicit none 
!       integer::n
!       integer::la,lc,lb,ld
!       integer::bsp_lmax
!       integer::k,kmin,kmax
!       integer::s2_ac,s2_bd
!       real(8)::special_gaunt
!
!       write(*,*)'lmax'
!       read(*,*)bsp_lmax
!      
!       n = 0
!       do la=0,bsp_lmax
!        do lc=0,bsp_lmax
!         do lb=0,bsp_lmax
!          do ld=0,bsp_lmax
!          !choosing the "k" expansion
!          kmin = max(abs(la-lc),abs(lb-ld))
!          kmax = min(abs(la+lc),abs(lb+ld))
!          do k=kmin,kmax
!
!             s2_ac = la + lc + k                      
!             s2_bd = lb + ld + k
!
!             if ((mod(s2_ac,2).eq.0).and.(mod(s2_bd,2).eq.0)) then
!             
!             n = n + 1 
!             print*, 'la,lb,lc,ld',la,lb,lc,ld
!             print*, 'k',k
!             print*, 'gaunt la,k,lc',special_gaunt(la,k,lc)
!             print*, 'gaunt lb,k,ld',special_gaunt(lb,k,ld)
!             print*,''
!             end if             
!          enddo
!          enddo
!         enddo
!        enddo
!       enddo
!       print*, ''
!       print*,'TOTAL INTEGRALES :',n
!       stop
!       end
!       
       
       function special_gaunt(l,k,lp)
       implicit none
       integer::l,lp,k
       real(8)::three0
       real(8)::special_gaunt



       call ThreeJBGN(dfloat(l),0.d0,dfloat(k),0.d0,&
        &dfloat(lp),0.d0,three0)

       special_gaunt = dsqrt(dfloat(2*l+1))*dsqrt(dfloat(2*lp+1))&
                    &*three0**2

       return
       end function special_gaunt
    


       subroutine ThreeJBGN(AJ,am,BJ,bm,J,m,threeJ)
!c
!c The routine is modified such that all inputs should
!c be real numbers. (And all values are given as the true value, 
!c i.e. NOT its double.)
!c
!c Calculates Clesbh-Gordan coefficient as described 
!c in the section about Testprogram and subroutines.
!c
!c Clebsch-Gordan coefficients, see
!c Brink & Satchler, Ref. [8] P.140-->special case
!c
!c Brink & Satchler, Ref. [8] P.34, eq.2.34-->general case
!c
      implicit none
      double precision cleb,threeJ
      integer ii,k,kmin,kmax
      double precision f1,f2,f3,t1,t2,t3,t4,t5,t6
      double precision d1,d2,d3,d4,g1,g2,g3,g4,g5,g6
      double precision AJ,am,BJ,bm,CJ,cm,J,m

      !.. Added by T.Kjellsson 2015
      double precision sqrtf3

      CJ=J
      cm=-m

 
      !check triangle inequality
      if (aj.gt.abs(bj+cj+.01).or.aj.lt.(abs(bj-cj)-.01)) then
         threeJ=0.0d0
         return
      end if

      f1=0.0d0
      kmin=int(max( 0.0d0,BJ-CJ-am,AJ-CJ+bm))
      kmax=int(min(AJ+BJ-CJ,AJ-am,BJ+bm))
      if (kmin .le.  kmax) then
         do ii=(kmin+1),(kmax+1)
            k=ii-1
            call factorial(dble(k),t1)
            call factorial(AJ+BJ-CJ-dble(k),t2)
            call factorial(AJ-am-dble(k),t3)
            call factorial(BJ+bm-dble(k),t4) 
            call factorial(CJ-BJ+am+dble(k),t5)
            call factorial(CJ-AJ-bm+dble(k),t6)
            f1=f1+(-1)**dble(k)/(t1*t2*t3*t4*t5*t6)
         end do
      end if
      
      call factorial(-CJ+AJ+BJ,d1)
      call factorial(CJ-AJ+BJ,d2)
      call factorial(CJ+AJ-BJ,d3)
      call factorial(1.0d0+CJ+AJ+BJ,d4)
      
      f2=d1*d2*d3/d4
      
      call factorial(CJ-cm,g1)
      call factorial(CJ+cm,g2)
      call factorial(AJ-am,g3)
      call factorial(AJ+am,g4)
      call factorial(BJ-bm,g5)
      call factorial(BJ+bm,g6)
      
      !.. Conditional statement inserted by T.Kjellsson
      !   f3 overflows maximum real number if AJ>49
      !..
      if(max(AJ,BJ,CJ) .gt. 49) then
         sqrtf3 = sqrt((1.d0+2.d0*CJ))*sqrt(g1)*sqrt(g2)*sqrt(g3)
         sqrtf3 = sqrtf3*sqrt(g4)*sqrt(g5)*sqrt(g6)
      else
         f3 = (1.d0+2.d0*CJ)*g1*g2*g3*g4*g5*g6
         sqrtf3 = sqrt(f3)
      end if

      if ((cm.ne.(am+bm)) .or. (CJ.lt.abs(AJ-BJ)) .or. (CJ.gt.(AJ+BJ))&
         &.or. (abs(am).gt.AJ) .or.(abs(bm).gt.BJ) .or.(abs(cm).gt.CJ)) then
         cleb=0.0d0
      else
         cleb=f1*sqrt(f2)*sqrtf3         
      end if
      
      if ((CJ.eq.0) .and. (cm.eq.0)) then 
         if (bm.eq.-am) then
            cleb=(-1.0d0)**(AJ-am)/sqrt(2.0d0*AJ+1.0d0)
         else
            cleb=0.0d0
         end if
      end if
     
      threeJ=cleb/((-1.0d0)**(-m+AJ-BJ)*dsqrt(2*J+1))

      return     
      end subroutine ThreeJBGN
      
     
       subroutine factorial(x,fact)

       implicit none
       double precision x,fact
       integer n,ii

       n=nint(x)
       fact=1.0d0
       if (n.lt.1) then
         fact=1.0d0
       else
         do ii=1,n
           fact=fact*ii
         end do
       end if

       return
       end subroutine factorial


   
