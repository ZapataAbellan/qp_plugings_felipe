BEGIN_PROVIDER [double precision, bsp_vee_full ,  (ao_dim,ao_dim,ao_dim,ao_dim)]
BEGIN_DOC
!Provides ALL the two-electron matrix elements :                           
END_DOC
 
implicit none      

!atomic orbitals (ao)
integer:: ao_p,ao_q,ao_t,ao_u 
!ao angular momentum
integer:: lp,lq,lt,lu
!ao magnetic momentum
integer:: mp,mq,mt,mu

!multipolar expansion 
integer:: k,kmax,kmin,mk

!angular coefficient 
double precision:: angular 

!Slater radial integral coefficient R^k(p,q,t,u)
double precision:: c

!integres for loops in the assembling algorithm...
integer::i,j,ip,jp,ii,jj,iv,jv,contador

!elements in "chemist notation"...
!electron   (1,1,2,2)
!index      (p,q,t,u)
bsp_vee_full(:,:,:,:) = 0.d0
!=======================
!two-electron integrals:
!=======================
do lp=0,bsp_lmax
 do mp=-lp,lp
  !
  do lq=0,bsp_lmax
   do mq=-lq,lq
    !
    do lt=0,bsp_lmax
     do mt=-lt,lt
      !
      do lu=0,bsp_lmax
       do mu=-lu,lu
        !
        kmin = max(abs(mp-mq),abs(mu-mt),abs(lp-lq),abs(lt-lu))
        kmax = min(abs(lp+lq),abs(lt+lu))
        !
        do k=kmin,kmax   

         if ((mod(lp+lq+k,2).eq.0).and.(mod(lt+lu+k,2).eq.0)) then
          
          angular = bsp_full_ck(lp,mp,lq,mq,k)*bsp_full_ck(lu,mu,lt,mt,k)

          if (angular.ne.0.d0) then

           !=================================
           !radial integration...Rk(p,q|t,u):
           !=================================
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
                 !
                 if( iv < jv ) then
                  c = roff1_kmax(ii,iv,k)*roff2_kmax(jj,jv,k) 
                 else if( iv > jv ) then
                  c = roff2_kmax(ii,iv,k)*roff1_kmax(jj,jv,k)
                 else
                  c = rdiag_kmax(ii,jj,iv,k)
                 end if
                 !                                                                   
                 !we remove from the basis the first and the last b-spline
                 if (((i +iv-1).gt.1).and.((i +iv-1).lt.bsp_number)) then 
                 if (((j +jv-1).gt.1).and.((j +jv-1).lt.bsp_number)) then
                 if (((ip+iv-1).gt.1).and.((ip+iv-1).lt.bsp_number)) then
                 if (((jp+jv-1).gt.1).and.((jp+jv-1).lt.bsp_number)) then
                 !                                                            
                 !ao index...
                 ao_p = ao_to_bspline(mp,lp,(i +iv-1)-1) !e(1)
                 ao_q = ao_to_bspline(mq,lq,(ip+iv-1)-1) !e(1)
                 ao_t = ao_to_bspline(mt,lt,(j +jv-1)-1) !e(2)
                 ao_u = ao_to_bspline(mu,lu,(jp+jv-1)-1) !e(2)
                 !
                 bsp_vee_full(ao_p,ao_q,ao_t,ao_u) += c * angular
                 !
                 endif
                 endif 
                 endif
                 endif
                 
                end do
               end do
              end do
             end do
            end do
           end do !loop : radial
           
          end if !if : angular.ne.0
         end if !if : triangle

        end do !loop : k
           
       end do !loop : mu
      end do !loop : lu
     end do !loop : mt
    end do !loop : lt
   end do !loop : mq
  end do !loop : lq
 end do !loop : mp
end do !loop : lp
!
END_PROVIDER 
