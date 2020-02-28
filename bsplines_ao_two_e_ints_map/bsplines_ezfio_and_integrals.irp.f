program bsplines_ezfio

 implicit none
 integer :: ao_p,ao_q,ao_t,ao_u
 integer :: rp,lp,mp
 integer :: rq,lq,mq
 integer :: rt,lt,mt
 integer :: ru,lu,mu

 integer :: k, kmax, kmin, mk
 double precision :: coef_1, coef_2, coef_3
 double precision :: fact
! double precision, allocatable :: pq_integrals(:,:)
! double precision :: tmp,get_ao_two_e_integral_bsplines,tmp2
 double precision :: tmp, tmp2  
! allocate(pq_integrals(ao_dim,ao_dim))
 
 write(*,*)'two-electron integrals...'
 call change_basis 
! provide ao_two_e_integrals_bsplines_in_map 
 !ao_p = 1
 !ao_q = 1
 
 
 write(*,*)'ao_dim',ao_dim
 write(*,*)'bsp_lmax',bsp_lmax 

 do ao_p=1, ao_dim
  do ao_q=1, ao_dim   
!   call give_all_bsplines_jl(ao_p,ao_q,pq_integrals)
   do ao_t=1, ao_dim
    do ao_u=1, ao_dim
! 
!         
         rp = ao_bspline_index_r(ao_p) 
         lp = ao_bspline_index_l(ao_p) 
         mp = ao_bspline_index_m(ao_p) 

         rq = ao_bspline_index_r(ao_q) 
         lq = ao_bspline_index_l(ao_q) 
         mq = ao_bspline_index_m(ao_q) 

         rt = ao_bspline_index_r(ao_t) 
         lt = ao_bspline_index_l(ao_t) 
         mt = ao_bspline_index_m(ao_t) 

         ru = ao_bspline_index_r(ao_u) 
         lu = ao_bspline_index_l(ao_u) 
         mu = ao_bspline_index_m(ao_u) 


         tmp  = bsp_vee_full(ao_p,ao_q,ao_t,ao_u)
         tmp2 = bsp_vee_full(ao_u,ao_t,ao_p,ao_q)

     if ( (dabs(tmp)-dabs(tmp2) .gt. 1.d-10)) then
! 
      write(55,'("ao_p, ao_q, ao_t, ao_u",4(1x,I6),3(1x,ES20.12E02))') ao_p, ao_q, ao_t, ao_u, tmp, tmp2
      write(55,'("rp  , rq  , rt  , ru  ",4(1x,I6))') rp, rq, rt, ru
      write(55,'("lp  , lq  , lt  , lu  ",4(1x,I6))') lp, lq, lt, lu
      write(55,'("mp  , mq  , mt  , mu  ",4(1x,I6))') mp, mq, mt, mu
      write(55,*)
!      write(55,'("Int:",1(1x,ES20.12E02))')tmp
      write(55,*)" "
!     
      end if
!     end if
!
!
    end do 
   end do
  end do
 end do

end


!     tmp = get_ao_two_e_integral_bsplines(ao_p,ao_t,ao_q,ao_u,ao_integrals_bsplines_map) 
!     tmp = pq_integrals(ao_t,ao_u)


!    do lp=0,bsp_lmax
!     do lq=0,bsp_lmax
!      do lt=0,bsp_lmax
!        do lu=0,bsp_lmax
!         !triangylar relation rules...
!         kmin = max(abs(lp-lq),abs(lt-lu))
!         kmax = min(abs(lp+lq),abs(lt+lu))
!         do k=kmin,kmax   
!          !triangular relations...
!          if ((mod(lp+lq+k,2).eq.0).and.(mod(lt+lu+k,2).eq.0)) then
!           !angular reduced elements...
!           coef_1 = dsqrt(dfloat(2*lp+1)*dfloat(2*lq+1))*bsp_3j(lp,0,lq,0,k,0)
!           coef_2 = dsqrt(dfloat(2*lt+1)*dfloat(2*lu+1))*bsp_3j(k,0,lt,0,lu,0)
!           !magnetic numbers...
!           do mp=-lp,lp
!            do mq=-lq,lq
!             do mt=-lt,lt
!              do mu=-lu,lu
!               !magnetic selection rules...
!               fact = 0.d0
!               if (((mp-mq)).eq.((mu-mt))) then
!                do mk=-k,k
!                 fact += (-1.d0)**(-mk-mp-mt)*bsp_3j(lp,-mp,lq,mq,k,mk)*bsp_3j(k,-mk,lt,-mt,lu,mu)
!                enddo
!                !total angular element...
!                coef_3 = fact*coef_1*coef_2 
!                write(55,'("k,lp,mp,lq,mq: (lp||lq)",5(1x,I2),1(1x,ES20.12E02))')k,lp,mp,lq,mq,coef_1
!                write(55,'("k,lt,mt,lu,mu: (lt||lu)",5(1x,I2),1(1x,ES20.12E02))')k,lt,mt,lu,mu,coef_2
!                write(55,'("k,lp,mp,lq,mq,lt,mt,lu,mu:",9(1x,I2),2(1x,ES20.12E02))')k,lp,mp,lq,mq,lt,mt,lu,mu,fact,coef_3
!                write(55,*)' '
!               end if !magnetic selection rules
!              end do !mu
!             end do !mt
!            end do !mq
!           end do !mp
!          end if !triangle relations...
!         end do !loop: k
!       
!       end do !loop:lu
!      end do !loop:lt
!     end do !loop:lq
!    end do !loop:lp
   



