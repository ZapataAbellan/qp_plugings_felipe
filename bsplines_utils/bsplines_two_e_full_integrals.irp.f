       BEGIN_PROVIDER [double precision, bsp_v_abcd,(bsp_total_dim,bsp_total_dim,bsp_total_dim,bsp_total_dim)]
       BEGIN_DOC
       !computes the two-electron matrix elements : 
       !
       !< a b | c d > = \sum_k R^k(a,b,c,d)*c^k(la,ma,lc,mc)*c^k(lb,mb,ld,md) 
       !
       END_DOC
       implicit none       
       !angular momentum 
       integer:: la,lb,lc,ld
       !coef. harmonic expansion 
       integer:: k,kmax,kmin               
       !triangle relation
       integer::s2_ac,s2_bd
       !block
       integer:: b
       !constant
       real(8):: c
       !integres for loops...
       integer::i,j,ip,jp,ii,jj,iv,jv
       !Gaunt coef. (special case m=mp=0)
       real(8)::special_gaunt
       real(8)::angular
       !
       !computing the integrals...
       !
       bsp_v_abcd(:,:,:,:) = 0.d0
       !
       !loops over the angular momentum
       do la=0,bsp_lmax
        do lc=0,bsp_lmax
         do lb=0,bsp_lmax
          do ld=0,bsp_lmax
         
          !choosing the "k" in the expansion
           kmin = max(abs(la-lc),abs(lb-ld))
           kmax = min(abs(la+lc),abs(lb+ld))

           do k=kmin,kmax
            
            s2_ac = la + lc + k                      
            s2_bd = lb + ld + k
                                                                  
            if ((mod(s2_ac,2).eq.0).and.(mod(s2_bd,2).eq.0)) then
            
            angular = special_gaunt(la,k,lc)*special_gaunt(lb,k,ld)

            print*,'la,lb,lc,ld',la,lb,lc,ld
            print*,'k,kmin,kmax',k,kmin,kmax
            print*,'angular',angular
          

            !loop over the Slater (radial) integrals 
!            do jv=1,bsp_nv
!             jj=0
!             do j=1,bsp_order
!              do jp=1,bsp_order
!               jj = jj + 1
!               do iv=1,bsp_nv
!               ii=0
!                do i=1,bsp_order
!                 do ip=1,bsp_order
!                 ii = ii + 1
!                 if( iv < jv ) then
!                  c = roff1_kmax(ii,iv,k+1)*roff2_kmax(jj,jv,k+1) 
!                 else if( iv > jv ) then
!                  c = roff2_kmax(ii,iv,k+1)*roff1_kmax(jj,jv,k+1)
!                 else
!                  c = rdiag_kmax(ii,jj,iv,k+1)
!                 end if
!                 bsp_v_abcd(i+iv-1,j+jv-1,ip+iv-1,jp+jv-1) += &
!                &  angular&
!                & *s_coef_norm(i+iv-1) &
!                & *s_coef_norm(j+jv-1) &
!                & *s_coef_norm(ip+iv-1)&
!                & *s_coef_norm(jp+jv-1)*c 
!                 end do
!                end do
!               end do
!              end do
!             end do
!            end do 
!            
            end if

           end do !loop: k

          end do
         end do
        end do
       end do


   END_PROVIDER 
