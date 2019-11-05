use map_module

subroutine give_all_two_e_ints_for_i(ao_p,bsp_vee_3_index)
 implicit none
use map_module
 integer, intent(in) :: ao_p
 double precision, intent(out) :: bsp_vee_3_index(ao_dim,ao_dim,ao_dim)
 BEGIN_DOC
 ! you enter with one ao_p index and you get out with 
 !
 ! the two-electron integrals : 
 !
 !bsp_vee_3_index(t,u,q) ===
 !
 !(t(2)u(2)|ao_p(1)q(1)) = sum_k R^k(ao_p(1)q(1);t(2)u(2)) sum_mk  <lpmp|C_mk^k(1)|lqmq> <ltmt|C_mk^k(2)|lumu>
 !
 END_DOC
      

 !atomic orbitals (ao)
 integer:: ao_q,ao_t,ao_u 
 !ao angular momentum
 integer:: lp,lq,lt,lu
 !ao magnetic momentum
 integer:: mp,mq,mt,mu

 !multipolar expansion 
 integer:: k,kmax,kmin,mk
 
 !Gaunt coefficient  <lm|C_mk^k|lpmp>
 double precision:: gaunt
 double precision:: angular

 !Slater radial integral coefficient R^k(p,q,t,u)
 double precision:: c

 !integres for loops in the assembling algorithm...
 integer::i,j,ip,jp,ii,jj,iv,jv
 !
 !Two-electron matrix elements in chemist notation...
 !electron    1 1 2 2           
 bsp_vee_3_index(:,:,:) = 0.d0

  !Loop over the angular space...
  !e(1)
   lp = ao_bspline_index_l(ao_p)
   mp = ao_bspline_index_m(ao_p)
    !e(1)
    do lq=0,bsp_lmax
     do mq=-lq,lq
      !e(2)
      do lt=0,bsp_lmax
       do mt=-lt,lt
        !e(2)
        do lu=0,bsp_lmax
         do mu=-lu,lu
         !
         !Multipolar expansion limits...
         kmin = max(abs(lp-lq),abs(lt-lu))
         kmax = min(abs(lp+lq),abs(lt+lu))
         ! 
         !Loop multipolar expansion...
         do k=kmin,kmax
          !
          !Loop over the radial grid...
          !Attention: 
          !R^k is stored in the chemist notation...
          !
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
                !we remove from the basis the first and the last b-spline...
                if (((i +iv-1).gt.1).and.((i +iv-1).lt.bsp_number)) then 
                if (((j +jv-1).gt.1).and.((j +jv-1).lt.bsp_number)) then
                if (((ip+iv-1).gt.1).and.((ip+iv-1).lt.bsp_number)) then
                if (((jp+jv-1).gt.1).and.((jp+jv-1).lt.bsp_number)) then
                 !
                 !ao index...
                 ao_q = ao_to_bspline(mq,lq,(ip+iv-1)-1) !e(1)
                 ao_t = ao_to_bspline(mt,lt,(j +jv-1)-1) !e(2)
                 ao_u = ao_to_bspline(mu,lu,(jp+jv-1)-1) !e(2)
                 !
                 !angular coefficient for a given k... 
                 angular = 0.d0
                 do mk=-k,k
                  angular +=  gaunt(lp,mp,k,mk,lq,mq) * gaunt(lt,mt,k,mk,lu,mu) 
                  !print*,'angular',angular
                 end do
                                
                 bsp_vee_3_index(ao_t,ao_u,ao_q) += c * angular
      
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
          !
         end do !loop: k
         !
        end do !loop:mu
       end do !loop:lu
       !
      end do !loop:mt
     end do !loop:lt
     !
    end do !loop:mq
   end do !loop:lq
      
      

 end

 subroutine fill_array_integrals_into_map(big_array,ao_p)
  use bitmasks
 implicit none
 use map_module
 integer, intent(in) :: ao_p
 double precision, intent(in) :: big_array(ao_dim,ao_dim,ao_dim)
 integer :: ao_q,ao_t,ao_u,index,n_integrals,size_buffer,i1,j1
 integer(key_kind),allocatable   :: buffer_i(:)
 real(integral_kind),allocatable :: buffer_value(:)

 size_buffer = min(ao_dim*ao_dim*ao_dim,16000000)

 n_integrals = 0

 allocate(buffer_i(size_buffer),buffer_value(size_buffer))
 do ao_q = ao_p, ao_dim
  j1 = ao_p+shiftr(ao_q*ao_q-ao_q,1)
  j1 = ao_q+shiftr(ao_p*ao_p-ao_p,1)
  do ao_u = 1, ao_dim           ! r1   
    i1 = shiftr(ao_u*ao_u-ao_u,1)
    if (i1 > j1) then
      exit
    endif
    do ao_t = 1, ao_u
      i1 += 1
      if (i1 > j1) then
        exit
      endif
      n_integrals += 1
      !
      buffer_value(n_integrals) = big_array(ao_t,ao_u,ao_q)
      call two_e_integrals_index(ao_p,ao_t,ao_q,ao_u,buffer_i(n_integrals))
      if (n_integrals == size_buffer) then
       call insert_into_ao_bsplines_integrals_map(n_integrals,buffer_i,buffer_value,&
           real(ao_integrals_threshold,integral_kind))
       n_integrals = 0
      endif
   enddo
  enddo
 enddo
 call insert_into_ao_bsplines_integrals_map(n_integrals,buffer_i,buffer_value,&
      real(ao_integrals_threshold,integral_kind))
 deallocate(buffer_i, buffer_value)

  
 end

 subroutine fill_all_integrals_into_ao_map
 implicit none
 integer :: ao_p
 double precision, allocatable :: bsp_vee_3_index(:,:,:)
 allocate(bsp_vee_3_index(ao_dim,ao_dim,ao_dim))
 do ao_p = 1, ao_dim
  call give_all_two_e_ints_for_i(ao_p,bsp_vee_3_index)
  call fill_array_integrals_into_map(bsp_vee_3_index,ao_p)
 enddo
!call map_merge(ao_bsplines_integrals_map)
 integer*8                      :: get_ao_bsplines_map_size, ao_bsplines_map_size
 ao_bsplines_map_size = get_ao_bsplines_map_size()
 print*,'ao_bsplines_map_size = ',ao_bsplines_map_size
 end
