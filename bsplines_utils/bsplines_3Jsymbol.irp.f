
       double precision function gaunt(l,m,k,mk,lp,mp)
       
       !gaunt = < l m | C_mk^k | l' m' > 

       implicit none
       integer::k,mk
       integer::l,m
       integer::lp,mp
       double precision::three0,threeJ
       double precision::coef
       double precision::coef_1
       double precision::coef_2

       !3j-symbols...
       three0 = bsp_3j_symbol(l,0,k,0,lp,0)
       threeJ = bsp_3j_symbol(l,-m,k,mk,lp,mp)

       !Coefficients... 
       coef   = (-1.d0)**m
       coef_1 = dsqrt(dfloat(2*l+1))
       coef_2 = dsqrt(dfloat(2*lp+1))

       !Gaunt coefficient...                                  
       !< l m | C_mk^k | l' m' > = (-1)^m [(2l+1)(2lp+1)]^1/2 ||(l;k;lp)(-m;mk;mp)|| ||(l;k;lp)(0;0;0)||

       gaunt = coef*coef_1*coef_2*three0*threeJ

       return
  end function gaunt
    
   

 BEGIN_PROVIDER [double precision,  bsp_3j_symbol , (0:bsp_lmax,-bsp_lmax:bsp_lmax,0:2*bsp_lmax,-2*bsp_lmax:2*bsp_lmax,0:bsp_lmax,-bsp_lmax:bsp_lmax)]

 implicit none
 
 double precision :: three3J
 integer:: lp,mp,lq,mq,k,mk

 three3J = 0.d0

 do lp=0,bsp_lmax
  do mp=-lp,lp
   !
   do lq=0,bsp_lmax
    do mq=-lq,lp
     !
     do k=0,2*bsp_lmax
      do mk=-k,k
      
       call ThreeJBGN(dfloat(lp),dfloat(mp),dfloat(k),dfloat(mk),dfloat(lq),dfloat(mq),three3J)

       bsp_3j_symbol(lp,mp,k,mk,lq,mq) = three3J

      end do 
     end do
    end do
   end do
  end do
 end do 

 END_PROVIDER


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



