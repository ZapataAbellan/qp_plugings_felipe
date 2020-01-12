      BEGIN_PROVIDER [ double precision, bsp_full_ck, (0:bsp_dlmax,-bsp_dlmax:bsp_dlmax,0:bsp_dlmax,-bsp_dlmax:bsp_dlmax,0:bsp_dlmax)]
      BEGIN_DOC
      !Provides ALL the spherical-harmonic full matrix elements < l m | C^k_(m-m') | l' m' > 
      !Those matrix elements are named as the Gaunt coefficients
      END_DOC
      
      implicit none 

      integer :: k
      integer :: l, m
      integer :: lp, mp
      
      do l=0,bsp_dlmax
       do lp=0,bsp_dlmax

        do k=0,bsp_dlmax

         do m=-l,l
          do mp=-lp,lp

           bsp_full_ck(l,m,lp,mp,k) = (-1.d0)**(l+m)*bsp_reduced_ck(l,lp,k)*bsp_racah_v(l,lp,k,-m,mp,m-mp)

          enddo
         enddo

        enddo

       enddo
      enddo

      END_PROVIDER 
       
      BEGIN_PROVIDER [double precision, bsp_reduced_ck, (0:bsp_qlmax,0:bsp_qlmax,0:bsp_qlmax)]
      BEGIN_DOC
      !Provides ALL the spherical-harmonic reduced matrix elements < l || C^k || l' > 
      END_DOC

      implicit none

      integer :: l, lp, k
 
      bsp_reduced_ck(:,:,:) = 0.d0

      do l =0,bsp_qlmax
       do lp=0,bsp_qlmax
        do k =0,bsp_qlmax

         bsp_reduced_ck(l,lp,k) = (-1.d0)**(-l)*dsqrt(dfloat(2*l+1)*dfloat(2*lp+1))*bsp_racah_v(l,lp,k,0,0,0)
    
        enddo
       enddo
      enddo

      END_PROVIDER
      
     
      
      BEGIN_PROVIDER [double precision, bsp_racah_v, (0:bsp_qlmax,0:bsp_qlmax,0:bsp_qlmax,-bsp_qlmax:bsp_qlmax,-bsp_qlmax:bsp_qlmax,-bsp_qlmax:bsp_qlmax)]
     
      BEGIN_DOC
      !Provides ALL the Racah V(a,b,c;ma,mb,mc) functions
      END_DOC

      implicit none
       
      integer :: a, ma
      integer :: b, mb
      integer :: c, mc
       
      bsp_racah_v(:,:,:,:,:,:) = 0.d0

      do a=0,bsp_qlmax
       do ma=-a,a
        !
        do b=0,bsp_qlmax
         do mb=-b,b
          !
          do c=0,bsp_qlmax
           do mc=-c,c


            bsp_racah_v(a,b,c,ma,mb,mc) = (-1.d0)**(c+b-a)*bsp_3j(a,ma,b,mb,c,mc)

           end do
          end do
         end do
        end do 
       end do
      end do

      END_PROVIDER


      BEGIN_PROVIDER [double precision, bsp_3j , (0:bsp_qlmax,-bsp_qlmax:bsp_qlmax,0:bsp_qlmax,-bsp_qlmax:bsp_qlmax,0:bsp_qlmax,-bsp_qlmax:bsp_qlmax)]

       BEGIN_DOC
       !Provides ALL the 3j-Wigner symbols
       END_DOC

       implicit none
       
       integer :: l, m
       integer :: k, mk
       integer :: lp, mp
       
       double precision :: tmp

       bsp_3j(:,:,:,:,:,:) = 0.d0
  
       do l=0,bsp_qlmax
        do m=-l,l
         !
         do k=0,bsp_qlmax
          do mk=-k,k
           !
           do lp=0,bsp_qlmax
            do mp=-lp,lp
              
             call ThreeJBGN(dfloat(l),dfloat(m),dfloat(k),dfloat(mk),dfloat(lp),dfloat(mp),tmp)
             
             bsp_3j(l,m,k,mk,lp,mp) = tmp
           
            end do
           end do
          end do
         end do 
        end do
       end do

       END_PROVIDER

       subroutine ThreeJBGN(AJ,am,BJ,bm,J,m,threeJ)
!c
!c Clebsch-Gordan coefficients, see
!c Brink & Satchler, Ref. [8] P.140-->special case
!c Brink & Satchler, Ref. [8] P.34, eq.2.34-->general case
!c
      implicit none
      double precision cleb,threeJ
      integer ii,k,kmin,kmax
      double precision f1,f2,f3,t1,t2,t3,t4,t5,t6
      double precision d1,d2,d3,d4,g1,g2,g3,g4,g5,g6
      double precision AJ,am,BJ,bm,CJ,cm,J,m

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



