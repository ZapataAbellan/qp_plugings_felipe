BEGIN_PROVIDER [ integer, bsplines_n_knot]
 implicit none
 bsplines_n_knot = bsplines_n_order_max + bsplines_number
END_PROVIDER 

BEGIN_PROVIDER [ integer, bsplines_n_interval]
 implicit none
 BEGIN_DOC
!  number of intervals for bspline grid
 END_DOC
 bsplines_n_interval = bsplines_number - bsplines_n_order_max + 1
END_PROVIDER 

BEGIN_PROVIDER [ double precision , bsplines_grid_vector, (bsplines_n_knot)]
 implicit none
 BEGIN_DOC
! Kuno's sequence ...
 END_DOC
 double precision  :: h                !initial step in the knot sequence for z*r
 integer           :: i
 !..determine h, the cobsplines_numbertant step size
 h=bsplines_box_size/dfloat(bsplines_n_interval)

 !..the multiple knots at the origin 
 bsplines_grid_vector(1:bsplines_n_order_max)=0.d0

 !.. the equally spaced points
 do i=bsplines_n_order_max+1,bsplines_number
    bsplines_grid_vector(i)=bsplines_grid_vector(i-1) + h
 end do

 !.. the multiple knots at the end of the grid
 bsplines_grid_vector(bsplines_number+1:bsplines_n_interval) = bsplines_grid_vector(bsplines_number)+ h


!do i=1,bsplines_number+bsplines_n_order_max
!   write(1,*)i,bsplines_grid_vector(i)
!end do

!print*,'bsplines_number      =',bsplines_number
!print*,'bsplines_n_order_max =',bsplines_n_order_max
!print*,'bsplines_n_interval  =',bsplines_n_interval
!print*,'h                    =',h
!print*,'bsplines_box_size    =',bsplines_box_size

END_PROVIDER 

subroutine pouet(a,b)
 implicit none
 double precision, intent(in) :: a
 double precision, intent(out) :: b(bsplines_n_interval)
 b(:) = 10.d0 * a
end

subroutine pouet_2(a,b)
 implicit none
 integer, intent(in)  :: a
 double precision, intent(out) :: b
 b = 10.d0 * dble(a)
end
