       BEGIN_PROVIDER [double precision, bsp_vee_full ,  (ao_dim,ao_dim,ao_dim,ao_dim)]
       BEGIN_DOC
       !Provides the two-electron matrix elements : 
       !
       !< p q | t u > = \sum_{k=0} R^k(p,q;t,u) \sum_{m_k=-k}^k delta_{m_k,m_p-m_t} delta_{m_k,m_q-m_u} 
       !              x (-1)^{m_k} c^k(l_p,m_p,l_t,m_t) c^k(l_q,m_q,l_u,m_u)
       !
       END_DOC
       implicit none       
       !atomic orbital angular momenta
       integer:: lp,lq,lt,lu
       integer:: mp,mq,mt,mu
       !multipolar expansion 
       integer:: k,kmax,kmin
       integer:: mk
       !atomic orbital index
       integer:: ao_p,ao_q,ao_t,ao_u
       !triangle relations
       integer:: ptk,quk
       !angular coefficient
       double precision:: angular
       !slater radial integral
       double precision:: c
       !integres for loops in the assembling algorithm...
       integer::i,j,ip,jp,ii,jj,iv,jv
       !
       !Gaunt coef. c^k(l_i,m_i,l_j,m_j)
       double precision::gaunt
       !
       !matrix elements
       bsp_vee_full(:,:,:,:) = 0.d0
       !
       !
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
           !limits in the multipolar expansion...
           kmin = max(abs(lp-lt),abs(lq-lu))
           kmax = min(abs(lp+lt),abs(lq+lu))
      
           !k-sum...
           do k=kmin,kmax
            
            !triangular relations between ao... 
            ptk = lp + lt + k                      
            quk = lq + lu + k
                                                                  
            if ((mod(ptk,2).eq.0).and.(mod(quk,2).eq.0)) then

             angular = 1.d0

             !print*,'lp,lq,lt,lu',lp,lq,lt,lu
             !print*,'k,kmin,kmax',k,kmin,kmax

             do mk=-k,k

              !delta functions...
              if ((mk.eq.(mp-mt)).and.(mk.eq.(mq-mu))) then
                
               angular = (-1.d0)**mk * gaunt(k,lp,mp,lt,mt) * gaunt(k,lq,mq,lu,mu) 

               !R^k

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
 
                     !we don't take the first and the last b-spline...
                     if (((i +iv-1).gt.1).and.((i +iv-1).lt.bsp_number)) then 
                     if (((j +jv-1).gt.1).and.((j +jv-1).lt.bsp_number)) then
                     if (((ip+iv-1).gt.1).and.((ip+iv-1).lt.bsp_number)) then
                     if (((jp+jv-1).gt.1).and.((jp+jv-1).lt.bsp_number)) then
 
                     ao_p = ao_to_bspline(mp,lp,(i +iv-1)-1)
                     ao_q = ao_to_bspline(mq,lq,(ip+iv-1)-1)
                     ao_t = ao_to_bspline(mt,lt,(j +jv-1)-1)
                     ao_u = ao_to_bspline(mu,lu,(jp+jv-1)-1)
                      
                     bsp_vee_full(ao_p,ao_q,ao_t,ao_u) += c * angular
                     
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

              end if !delta functions
             end do !loop: mk
            end if !triangular relations
           end do !loop: k
         
         end do
        end do
       end do
      end do
     end do
    end do
   end do
  end do


  END_PROVIDER 
