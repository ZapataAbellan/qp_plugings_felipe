program bsplines_utils
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  print *,'bsplines_n_order_max',bsplines_n_order_max
  integer :: i
! provide bsplines_grid_vector
!do i = 1, bsplines_n_interval
 i = 3
   print*,'bsplines_grid_vector(i) = ',bsplines_grid_vector(i)
!enddo
end
