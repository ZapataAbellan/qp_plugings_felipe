BEGIN_PROVIDER [ integer, bsp_atomic_charge  ]
  implicit none
  BEGIN_DOC
! atomic charge
  END_DOC
  bsp_atomic_charge = int(nucl_charge(1))
END_PROVIDER 
