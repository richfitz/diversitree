      subroutine DNCHBV( m, t, H,ldh, y, wsp )

      implicit none
      integer          m, ldh
      double precision t, H(ldh,m), y(m)
      complex*16       wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  DNCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is assumed to be upper-Hessenberg.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise. 
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the Hessenberg matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) upper Hessenberg matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine (thus avoiding an idle complex array elsewhere)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO
      integer ndeg, i, j, k, ip, ih, iy, iz
      parameter ( ndeg=7, ZERO=(0.0d0,0.0d0) )
      double precision alpha0
      complex*16 alpha(ndeg), theta(ndeg), tmpc

      intrinsic ABS,DBLE,MIN
      
*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*     
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      enddo
      do ip = 1,ndeg
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            wsp(iy+j-1) = wsp(iz+j-1)
            do i = 1,MIN( j+1,m )
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            enddo
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            do k = i,m
               wsp(ih+(j-1)*m+k-1) = ZERO
            enddo
         enddo
         do i = 1,m-1
*---        Get pivot and exchange rows ...
            if (ABS(wsp(ih+(i-1)*m+i-1)).lt.ABS(wsp(ih+(i-1)*m+i))) then
               call ZSWAP( m-i+1, wsp(ih+(i-1)*m+i-1),m, 
     .                     wsp(ih+(i-1)*m+i),m )
               call ZSWAP( 1, wsp(iy+i-1),1, wsp(iy+i),1 )
            endif
*---        Forward eliminiation ... 
            tmpc = wsp(ih+(i-1)*m+i) / wsp(ih+(i-1)*m+i-1)
            call ZAXPY( m-i, -tmpc, wsp(ih+i*m+i-1),m, wsp(ih+i*m+i),m )
            wsp(iy+i) = wsp(iy+i) - tmpc*wsp(iy+i-1)
         enddo
*---     Backward substitution ...    
         do i = m,1,-1
            tmpc = wsp(iy+i-1)
            do j = i+1,m
               tmpc = tmpc - wsp(ih+(j-1)*m+i-1)*wsp(iy+j-1)
            enddo
            wsp(iy+i-1) = tmpc / wsp(ih+(i-1)*m+i-1)
         enddo
*---     Accumulate the partial result in y ...     
         do j = 1,m
            y(j) = y(j) + DBLE( alpha(ip)*wsp(iy+j-1) )
         enddo
      enddo
      END
