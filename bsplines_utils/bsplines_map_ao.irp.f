use map_module

!! AO Map
!! ======

BEGIN_PROVIDER [ type(map_type), ao_bsplines_integrals_map ]
  implicit none
  BEGIN_DOC
  ! AO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call two_e_integrals_index(ao_dim,ao_dim,ao_dim,ao_dim,key_max)
  sze = key_max
  call map_init(ao_bsplines_integrals_map,sze)
  print*,  'AO map initialized : ', sze
END_PROVIDER

function get_ao_bsplines_map_size()
  implicit none
  integer (map_size_kind) :: get_ao_bsplines_map_size
  BEGIN_DOC
  ! Returns the number of elements in the AO map
  END_DOC
  get_ao_bsplines_map_size = ao_bsplines_integrals_map % n_elements
end

subroutine clear_ao_bsplines_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the MO map
  END_DOC
  call map_deinit(ao_bsplines_integrals_map)
  FREE ao_bsplines_integrals_map 
end


subroutine insert_into_ao_bsplines_integrals_map(n_integrals,                 &
      buffer_i, buffer_values, thr)
  use map_module
  implicit none

  BEGIN_DOC
  ! Create new entry into MO map, or accumulate in an existing entry
  END_DOC

  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  real(integral_kind), intent(in)    :: thr
  call map_update(ao_bsplines_integrals_map, buffer_i, buffer_values, n_integrals, thr)
end


double precision function get_ao_bsplines_two_e_integral(i,j,k,l,map) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one AO bi-electronic integral from the AO map
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  real(integral_kind)            :: tmp
  PROVIDE ao_bisplines_integrals_in_map
  !DIR$ FORCEINLINE
  call two_e_integrals_index(i,j,k,l,idx)
  !DIR$ FORCEINLINE
  call map_get(map,idx,tmp)
  result = tmp
end

BEGIN_PROVIDER [ logical, ao_bisplines_integrals_in_map ]
  use map_module
  implicit none
  BEGIN_DOC
  ! If True, the map of MO two-electron integrals is provided
  END_DOC
  provide ao_bsplines_integrals_map
  ao_bisplines_integrals_in_map = .True.
  call add_ao_bsplines_integrals_to_map

END_PROVIDER

BEGIN_PROVIDER [double precision,  count_array, (ao_dim,ao_dim,ao_dim,ao_dim)]
 implicit none
 count_array = 0.d0

END_PROVIDER 

subroutine add_ao_bsplines_integrals_to_map
  use bitmasks
  implicit none

  BEGIN_DOC
  ! Adds integrals to the AO map for the bspline basis
  END_DOC

  integer                         :: n_integrals
  integer                         :: size_buffer
  integer(key_kind),allocatable   :: buffer_i(:)
  real(integral_kind),allocatable :: buffer_value(:)

  size_buffer = min(ao_dim*ao_dim*ao_dim,16000000)

  n_integrals = 0
  allocate(buffer_i(size_buffer),buffer_value(size_buffer),key_ao_1(N_int),key_ao_2(N_int),key(N_int))

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
  double precision::gaunt
  integer(bit_kind), allocatable :: key_ao_1(:),key_ao_2(:),key(:)

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

        !mk sum...
        do mk=-k,k

         !delta functions...
         if ((mk.eq.(mp-mt)).and.(mk.eq.(mq-mu))) then
         
          angular = (-1.d0)**mk * gaunt(k,lp,mp,lt,mt) * gaunt(k,lq,mq,lu,mu) 

          !R^k(ao_p,ao_q,ao_t,ao_u)
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
 
                

                ! (pq|tu) 
                ! <pt|qu>
                ao_p = ao_to_bspline(mp,lp,(i +iv-1)-1) ! 1
                ao_q = ao_to_bspline(mq,lq,(ip+iv-1)-1) ! 1
                ao_t = ao_to_bspline(mt,lt,(j +jv-1)-1) ! 2
                ao_u = ao_to_bspline(mu,lu,(jp+jv-1)-1) ! 2

                key_ao_1 = 0_bit_kind
                key_ao_2 = 0_bit_kind
                key      = 0_bit_kind
                call set_bit_to_integer(ao_p,key_ao_1,N_int)
                call set_bit_to_integer(ao_q,key_ao_1,N_int)
                call set_bit_to_integer(ao_t,key_ao_2,N_int)
                call set_bit_to_integer(ao_u,key_ao_2,N_int)

                call set_bit_to_integer(ao_p,key,N_int)
                call set_bit_to_integer(ao_t,key,N_int)
                call set_bit_to_integer(ao_q,key,N_int)
                call set_bit_to_integer(ao_u,key,N_int)

                integer ::icount_1,icount_2,icount
                integer ::ik
                double precision :: factor
                icount = 0
                icount_1 = 0
                icount_2 = 0
                do ik = 1, N_int
                 icount_1 += popcnt(key_ao_1(ik))
                 icount_2 += popcnt(key_ao_2(ik))
                 icount += popcnt(key(ik))
                enddo
                
                factor = 1.d0
                if(icount_1 == 1 .and. icount_2 == 1 .and. icount == 2)then
                 factor = 0.5d0
                else if (icount_1 == 2 .and. icount_2 == 2 .and. icount == 3)then
                 factor = 1.d0/8.d0 
                else if (icount_1 == 1 .and. icount_2 == 2 .and. icount == 2)then
                 factor = 0.25d0
                else if (icount_1 == 2 .and. icount_2 == 1 .and. icount == 2)then
                 factor = 0.25d0
                else if (icount_1 == 2 .and. icount_2 == 2 .and. icount == 2)then
                 factor = 0.25d0
                else if (icount_1 == 1 .and. icount_2 == 2 .and. icount == 3)then
                 factor = 0.25d0
                else if (icount_1 == 2 .and. icount_2 == 1 .and. icount == 3)then
                 factor = 0.25d0
                else if (icount == 4)then
                 factor = 1.d0/8.d0
                endif

                ! (pq|tu) 

                n_integrals += 1
!               buffer_value(n_integrals) = c * angular * factor 
                buffer_value(n_integrals) = c * factor 
                call two_e_integrals_index(ao_p,ao_t,ao_q,ao_u,buffer_i(n_integrals))
                if (n_integrals == size_buffer) then
                  call insert_into_ao_bsplines_integrals_map(n_integrals,buffer_i,buffer_value,&
                      real(ao_integrals_threshold,integral_kind))
                  n_integrals = 0
                endif
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

     end do !loop:mu
    end do !loop:lu
        !
   end do !loop:mt
  end do !loop:lt
      !
 end do !loop:mq
 end do !loop:lq
    !
 end do !loop:mp
 end do !loop:lp


  call insert_into_ao_bsplines_integrals_map(n_integrals,buffer_i,buffer_value,&
      real(ao_integrals_threshold,integral_kind))
  deallocate(buffer_i, buffer_value)
  call map_merge(ao_bsplines_integrals_map)

  integer*8                      :: get_ao_bsplines_map_size, ao_bsplines_map_size
  ao_bsplines_map_size = get_ao_bsplines_map_size()
  print*,'ao_bsplines_map_size = ',ao_bsplines_map_size


end
