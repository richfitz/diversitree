      subroutine dgcoov ( x, y, ia, ja, a, n, nz )
      implicit none
      double precision x(*), y(*)
*
*---  Computes y = A*x.  A is assumed here to be under the COOrdinates
*---  storage format, and is contained in arguments ia, ja, a
*
      integer n, nz, nzmax
      parameter( nzmax = 139301 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      integer i, j
 
      do j = 1,n
         y(j) = 0.0d0
      enddo
      do i = 1,nz
         y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
      enddo
      END
