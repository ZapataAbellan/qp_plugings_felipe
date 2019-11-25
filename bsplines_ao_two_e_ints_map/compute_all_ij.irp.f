
subroutine compute_ao_integrals_bsplines_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC

  integer, intent(in)            :: j,l
  integer,intent(out)            :: n_integrals
  integer(key_kind),intent(out)  :: buffer_i(ao_dim*ao_dim)
  real(integral_kind),intent(out) :: buffer_value(ao_dim*ao_dim)

  integer                        :: i,k
  double precision               :: ao_two_e_integral_bsplines,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  double precision               :: thr
  integer                        :: kk, m, j1, i1
  double precision, allocatable  ;; jl_integrals(:,:)
  allocate(jl_integrals(ao_dim,ao_dim))

  thr = ao_integrals_threshold

  n_integrals = 0
  
  call give_all_bsplines_kl(j,l,jl_integrals)
  j1 = j+ishft(l*l-l,-1)
  do k = 1, ao_dim           ! r1
    i1 = ishft(k*k-k,-1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif
      !DIR$ FORCEINLINE
      integral = jl_integrals(i,k)  ! i,k : r1    j,l : r2
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral
    enddo
  enddo
end

