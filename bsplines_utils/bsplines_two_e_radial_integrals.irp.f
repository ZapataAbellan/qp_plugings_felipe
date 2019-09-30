       BEGIN_PROVIDER [double precision, roff1_kmax, (bsp_order*bsp_order,bsp_nv, bsp_lmax+1)]
      &BEGIN_PROVIDER [double precision, roff2_kmax, (bsp_order*bsp_order,bsp_nv, bsp_lmax+1)]
      &BEGIN_PROVIDER [double precision, rdiag_kmax, (bsp_order*bsp_order,bsp_order*bsp_order,bsp_nv, bsp_lmax+1)]
       BEGIN_DOC 
       !This provides the radial Slater Integrals up tp "k_max" for the Coulomb potential
       END_DOC
       implicit none 
       !
       !Gauss-Legendre parameters:
       integer::v 
       real(8)::rvleft,rvright
       real(8),allocatable::ri(:),wi(:)
       real(8),allocatable::rii(:),wii(:) 
       !
       !Carl de Boor's parameters for b-splines:
       integer,parameter::nderiv=1
       integer::left,mflag,leftmki,leftmkj
       real(8),allocatable::valuesiv(:,:)
       real(8),allocatable::valuesjv(:,:)
       real(8),allocatable::work(:,:)
       !
       !Cell Algorithm required matrices:
       real(8)::c
       real(8),allocatable::rkoffdiv(:,:)
       real(8),allocatable::rkoffdjv(:,:)
       real(8),allocatable::rkdiag(:,:,:)
       real(8),allocatable::rkdiagtmp(:,:)
       real(8),allocatable::rkdiagtri(:,:)
       !
       !Coulomb Potential Function:
       integer:: k
       real(8)::pot_c_maxr
       real(8)::pot_c_minr
       !
       !
       !Integers for loops:
       integer::i,ip,j,jp
       integer::iv,jv
       integer::rglpiv,rglpjv
       integer::ii,jj
       integer::ih,jh,ihp,jhp
       !
       !Allocation of matrices:
       allocate(ri(bsp_glp))
       allocate(wi(bsp_glp))
       allocate(rii(bsp_glp))
       allocate(wii(bsp_glp))
       !
       allocate(work(bsp_order,bsp_order))
       !
       allocate(valuesiv(bsp_order,nderiv))
       allocate(valuesjv(bsp_order,nderiv))
       !
       allocate(rkoffdiv(bsp_order*bsp_order,bsp_nv))
       allocate(rkoffdjv(bsp_order*bsp_order,bsp_nv))
       allocate(rkdiag(bsp_order*bsp_order,bsp_order*bsp_order,bsp_nv))
       !
       allocate(rkdiagtmp(bsp_order,bsp_order))
       allocate(rkdiagtri(bsp_order*bsp_order,bsp_order*bsp_order))
       !
       !
       !
       !
       !Expansion term "K"
       do k=0,bsp_lmax
        !
        !Cell Integration...
        !
        do v=1,bsp_nv
          !
          rkoffdiv(:,v)=0.d0
          rkoffdjv(:,v)=0.d0
          rkdiag(:,:,v)=0.d0                                                                                      
          rkdiagtmp(:,:)=0.d0 
          rkdiagtri(:,:)=0.d0
          !
          rvleft=bsp_grid_knot(bsp_order+v-1)
          rvright=bsp_grid_knot(bsp_order+v)
          !
          call gauleg(rvleft,rvright,ri,wi,bsp_glp)
          do rglpiv=1,bsp_glp
             valuesiv(:,:)=0.d0
             call interv(bsp_grid_knot,bsp_number_knot,ri(rglpiv),left,mflag)
             leftmki = left - bsp_order
             call bsplvd(bsp_grid_knot,bsp_number_knot,bsp_order,ri(rglpiv),left,work,valuesiv,nderiv)
             !
             !Off-diagonal cells integrals...
             ii=0
             do i=1,bsp_order
                do ip=1,bsp_order
                   ii=ii+1
                   rkoffdiv(ii,v)= rkoffdiv(ii,v) + &
!                  &wi(rglpiv)*valuesiv(i,1)*(ri(rglpiv)**dfloat(k))*valuesiv(ip,1)
                   &wi(rglpiv)*valuesiv(i,1)*pot_c_minr(k,ri(rglpiv))*valuesiv(ip,1)
                   !
                   rkoffdjv(ii,v) =rkoffdjv(ii,v) + &
!                  &wi(rglpiv)*valuesiv(i,1)*(1.d0/ri(rglpiv)**(dfloat(k)+1))*valuesiv(ip,1)
                   &wi(rglpiv)*valuesiv(i,1)*pot_c_maxr(k,ri(rglpiv))*valuesiv(ip,1)
                end do
             end do
             !
             !Over-diagonal cells integrals...
             call gauleg(rvleft,ri(rglpiv),rii,wii,bsp_glp)
             do rglpjv=1,bsp_glp
                valuesjv(:,:)=0.d0
                call interv(bsp_grid_knot,bsp_number_knot,rii(rglpjv),left,mflag)
                leftmkj = left - bsp_order
                call bsplvd(bsp_grid_knot,bsp_number_knot,bsp_order,rii(rglpjv),left,work,valuesjv,nderiv)
                do j=1,bsp_order
                   do jp=1,bsp_order
                      rkdiagtmp(j,jp) = rkdiagtmp(j,jp) + &
!                     &wii(rglpjv)*valuesjv(j,1)*(rii(rglpjv)**dfloat(k))*valuesjv(jp,1)
                      &wii(rglpjv)*valuesjv(j,1)*pot_c_minr(k,rii(rglpjv))*valuesjv(jp,1)
                   end do
                end do
             enddo
             !Cell-triangle integrals... 
             ii=0
             do i=1,bsp_order
                do ip=1,bsp_order
                   ii=ii+1
                   jj=0
                   do j=1,bsp_order
                      do jp=1,bsp_order
                         jj=jj+1
                         rkdiagtri(ii,jj)=  rkdiagtri(ii,jj) + &
!                        &wi(rglpiv)*valuesiv(i,1)*(1.d0/ri(rglpiv)**(dfloat(k)+1))*valuesiv(ip,1)*rkdiagtmp(j,jp)
                         &wi(rglpiv)*valuesiv(i,1)*pot_c_maxr(k,ri(rglpiv))*valuesiv(ip,1)*rkdiagtmp(j,jp)
                      end do
                   end do
                end do
             end do
             !
             rkdiagtmp(:,:)=0.d0
             !
          enddo ! end loop : rglpiv
          !Sum up the two cell-triangle integrals...
          rkdiag(1:bsp_order*bsp_order,1:bsp_order*bsp_order,v)=rkdiagtri + transpose(rkdiagtri)
          
       enddo ! end loop : v
       !
       !
       roff1_kmax(:,:,k+1)  = rkoffdiv(:,:)     
       roff2_kmax(:,:,k+1)  = rkoffdjv(:,:)    
       rdiag_kmax(:,:,:,k+1)= rkdiag(:,:,:)    
       !
       enddo ! end loop : k

       END_PROVIDER

       !===============
       !COULOMB POTENTIAL 
       !===============
       function POT_C_MAXR(k,x)
       !===============
       implicit none
       real(8)::POT_C_MAXR
       integer::k
       real(8)::x
       POT_C_MAXR=x**(-(k+1))
       return
       end function POT_C_MAXR
       !
       !===============
       function POT_C_MINR(k,x)
       !===============
       implicit none
       real(8)::POT_C_MINR
       integer::k
       real(8)::x
       POT_C_MINR=x**k
       return
       end function POT_C_MINR
       !
       !
       !
