*     This is a simple wrapper to a couple of the expokit functions from
*     within R.  The code here could easily have been written in C, and
*     that is the approach taken by mkn-expokit.c

*     I'll hard code in some upper bounds, I think.  Let nmax = 1024
*     (2^10), and nzmax = 102400 (100 nonzero elements per state).
*
*     Note that there is really no reason to call this from within
*     Fortran; this would be about as easy from C.
*     **** Sparse
      subroutine DSEXPMV(Q, n, ia, ja, nz, qnorm, v, t, tol, out, iflag)
      implicit none
      integer n, ia(nz), ja(nz), nz
      double precision Q(nz), qnorm, v(n), t

      integer nmax, nzmax, mmax
      parameter( nmax=1024, nzmax=102400, mmax=30 )
      integer lwsp, liwsp
*      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax+2 )
      parameter(lwsp=nmax*(mmax+1)+nmax+5*(mmax+2)**2+7,liwsp=nmax+2)

      integer m, itrace, iflag, iwsp(liwsp)
      double precision tol, anorm, wsp(lwsp), scal

      double precision out(n)

      if ( n .gt. nmax .or. nz .gt. nzmax ) then
         iflag = -1
         return
      endif

      itrace = 0
      iflag = 0
      scal = 1.0
*     Krylov basis.  This was the default in BiSSE/unresolved
      m = 15

      call DSEXPV(n, m, t, v, out, tol,
     .     qnorm, ia, ja, Q, nz, wsp,lwsp, iwsp,liwsp, itrace,
     .     iflag, scal)

      end

*     **** Dense
      subroutine DDEXPMV(Q, n, v, t, out, iflag)
      implicit none
      integer n, iflag
      double precision t, Q(n,n), v(n), out(n)

      double precision alpha, beta
      integer          iZERO, iONE
      integer ideg, iexp, ns, lwsp, nmax
      parameter( nmax=64 )
      parameter( ideg=6, lwsp=4*nmax*nmax+ideg+1 )
      parameter( alpha=1.0d0, beta=0.0d0, iZERO=0, iONE=1 )

      double precision wsp(lwsp)
      integer ipiv(n)

*---  compute E = exp(tQ) via Expokit's Pade algorithm
      call DGPADM(ideg, n, t, Q,n, wsp, lwsp, ipiv, iexp, ns, iflag)

*---  multiply out = E*v via BLAS (E is at wsp(iexp))
*     compute y = alpha*A*x + beta * y
*              TRANS  M x N  ALPHA A         LDA  X  INCX  BETA  Y    INCY
      call DGEMV('n', n,  n, alpha,wsp(iexp), n,  v, 1,    beta, out, 1)
      
      end

*     **** Full
      subroutine DEXPMF(Q, n, t, out, iflag)
      implicit none
      integer n, iflag
      double precision t, Q(n,n), v(n), out(n*n)

      double precision ZERO, ONE
      integer ideg, iexp, ns, lwsp, nmax
      parameter( nmax=64 )
      parameter( ideg=6, lwsp=4*nmax*nmax+ideg+1 )
      parameter( ZERO=0.0d0, ONE=1.0d0 )

      double precision wsp(lwsp)
      integer i, ipiv(n)

*---  compute E = exp(tQ) via Expokit's Pade algorithm
      call DGPADM(ideg, n, t, Q,n, wsp, lwsp, ipiv, iexp, ns, iflag)

      do i=1,(n*n)
         out(i) = wsp(iexp+i-1)
      enddo

      end

*     Now add checkpointing:
      subroutine DSEXPMVI(Q, n, ia, ja, nz, qnorm, v, t,lt, tol, 
     .     out, iflag)
      implicit none
      integer n, ia(nz), ja(nz), nz, lt
      double precision Q(nz), qnorm, v(n), t

      integer nmax, nzmax, mmax
      parameter( nmax=1024, nzmax=102400, mmax=30 )
      integer lwsp, liwsp
      parameter(lwsp=nmax*(mmax+1)+nmax+5*(mmax+2)**2+7,liwsp=nmax+2)

      integer m, itrace, iflag, iwsp(liwsp)
      double precision tol, anorm, wsp(lwsp), scal

      double precision out(n*lt)

      itrace = 0
      iflag = 0
      scal = 1.0
*     Krylov basis.  This was the default in BiSSE/unresolved
      m = 15

      call DSEXPVI(n, m, t,lt, v, out, tol,
     .     qnorm, ia, ja, Q, nz, wsp,lwsp, iwsp,liwsp, itrace,
     .     iflag, scal)

      end
