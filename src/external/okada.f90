MODULE okada
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
CONTAINS

  SUBROUTINE dc3d_wrapper(&
       & alpha, &
       & x, y, z, &
       & depth, dip, &
       & al1, al2, &
       & aw1, aw2, &
       & disl1, disl2, disl3, &
       & ux, uy, uz, &
       & uxx, uyx, uzx, &
       & uxy, uyy, uzy, &
       & uxz, uyz, uzz, &
       & iret) BIND(C, NAME='__dc3d__')

    REAL*8 :: &
         & alpha, &
         & x, y, z, &
         & depth, dip, &
         & al1, al2, &
         & aw1, aw2, &
         & disl1, disl2, disl3, &
         & ux, uy, uz, &
         & uxx, uyx, uzx, &
         & uxy, uyy, uzy, &
         & uxz, uyz, uzz

    INTEGER*8 :: iret

    CALL dc3d(&
         & alpha, &
         & x, y, z, &
         & depth, dip, &
         & al1, al2, &
         & aw1, aw2, &
         & disl1, disl2, disl3, &
         & ux, uy, uz, &
         & uxx, uyx, uzx, &
         & uxy, uyy, uzy, &
         & uxz, uyz, uzz, &
         & iret)

  END SUBROUTINE dc3d_wrapper

END MODULE okada
