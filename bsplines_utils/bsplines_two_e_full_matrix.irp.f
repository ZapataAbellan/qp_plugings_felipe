      BEGIN_PROVIDER [double precision, bsp_vee_full ,  (ao_dim,ao_dim,ao_dim,ao_dim)]
      BEGIN_DOC
      !Provides the two-electron matrix elements :                           
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
      double precision:: coef_1, coef_2, coef_3 
      
      !phase factor
      double precision:: fact

      !Slater radial integral coefficient R^k(p,q,t,u)
      double precision:: c

      !integres for loops in the assembling algorithm...
      integer::i,j,ip,jp,ii,jj,iv,jv,contador
      
      !Two-electron matrix elements in chemist notation...
      !electron    1 1 2 2           
      bsp_vee_full(:,:,:,:) = 0.d0

      do lp=0,bsp_lmax
       do lq=0,bsp_lmax
        do lt=0,bsp_lmax
         do lu=0,bsp_lmax
          !triangylar relation rules...
          kmin = max(abs(lp-lq),abs(lt-lu))
          kmax = min(abs(lp+lq),abs(lt+lu))
          do k=kmin,kmax   
           !triangular relations...
           if ((mod(lp+lq+k,2).eq.0).and.(mod(lt+lu+k,2).eq.0)) then
            !angular reduced elements...
            coef_1 = dsqrt(dfloat(2*lp+1)*dfloat(2*lq+1))*bsp_3j(lp,0,lq,0,k,0)
            coef_2 = dsqrt(dfloat(2*lt+1)*dfloat(2*lu+1))*bsp_3j(k,0,lt,0,lu,0)
            !magnetic numbers...
            do mp=-lp,lp
             do mq=-lq,lq
              do mt=-lt,lt
               do mu=-lu,lu
                !magnetic selection rules...
                fact = 0.d0
                if (((mp-mq)).eq.((mu-mt))) then
                 do mk=-k,k
                  fact += (-1.d0)**(-mk-mp-mt)*bsp_3j(lp,-mp,lq,mq,k,mk)*bsp_3j(k,-mk,lt,-mt,lu,mu)
                 enddo
                 !total angular element...
                 coef_3 = fact*coef_1*coef_2 
!                 write(55,'("k,lp,mp,lq,mq: (lp||lq)",5(1x,I2),1(1x,ES20.12E02))')k,lp,mp,lq,mq,coef_1
!                 write(55,'("k,lt,mt,lu,mu: (lt||lu)",5(1x,I2),1(1x,ES20.12E02))')k,lt,mt,lu,mu,coef_2
!                 write(55,'("k,lp,mp,lq,mq,lt,mt,lu,mu:",9(1x,I2),2(1x,ES20.12E02))')k,lp,mp,lq,mq,lt,mt,lu,mu,fact,coef_3
!                 write(55,*)' '
                  !radial integration...Rk(p,q|t,u)
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
                        bsp_vee_full(ao_p,ao_q,ao_t,ao_u) += c*coef_3
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
                  end do

                end if !magnetic selection rules
               end do !mu
              end do !mt
             end do !mq
            end do !mp

           end if !triangle relations...
          end do !loop: k      
         end do !loop:lu
        end do !loop:lt
       end do !loop:lq
      end do !loop:lp
      !
      END_PROVIDER 
