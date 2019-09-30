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
!  dimension of the b-spline basis set 
 END_DOC
 bsp_dim = bsp_number - 2
END_PROVIDER 

BEGIN_PROVIDER [integer, ao_dim]
 implicit none
 BEGIN_DOC
 ! dimension of the atomic orbital basis
 END_DOC
 integer :: l
 ao_dim = bsp_dim
 do l = 1,bsp_lmax
  ao_dim += bsp_dim*(2 * l+1)
 enddo
END_PROVIDER



