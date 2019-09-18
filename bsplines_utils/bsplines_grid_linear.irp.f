BEGIN_PROVIDER [ integer, bsp_number_knot]
 implicit none
BEGIN_DOC
! number of knots 
END_DOC
 bsp_number_knot = bsp_order + bsp_number
END_PROVIDER 

BEGIN_PROVIDER [ integer, bsp_nv]
 implicit none
 BEGIN_DOC
!  number of intervals in the grid
 END_DOC
 bsp_nv = bsp_number - bsp_order + 1
END_PROVIDER 

BEGIN_PROVIDER [ integer, bsp_dim]
 implicit none
 BEGIN_DOC
!  matrix dimensionm dimension of the used basis 
 END_DOC
 bsp_dim = bsp_number - 2
END_PROVIDER 

BEGIN_PROVIDER [ double precision , bsp_grid_knot, (bsp_number_knot)]
 implicit none
 BEGIN_DOC
! Vector containing the Knot sequence
 END_DOC
 double precision  :: h                !initial step in the knot sequence for z*r
 integer           :: i
 !..determine h, the cobsp_numbertant step size
 h=bsp_box_size/dfloat(bsp_nv)

 !..the multiple knots at the origin 
 bsp_grid_knot(1:bsp_order)=0.d0

 !.. the equally spaced points
 do i=bsp_order+1,bsp_number
    bsp_grid_knot(i)=bsp_grid_knot(i-1) + h
 end do

 !.. the multiple knots at the end of the grid
 bsp_grid_knot(bsp_number+1:bsp_number_knot) = bsp_grid_knot(bsp_number)+ h


!do i=1,bsp_number+bsp_order
!   write(1,*)i,bsp_grid_knot(i)
!end do

!print*,'bsp_number      =',bsp_number
!print*,'bsp_order =',bsp_order
!print*,'bsp_nv  =',bsp_nv
!print*,'h                    =',h
!print*,'bsp_box_size    =',bsp_box_size

END_PROVIDER 
