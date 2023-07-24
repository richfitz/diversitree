*----------------------------------------------------------------------|
      subroutine DMEXPV( n, m, t, v, w, tol, anorm, ia, ja, a, nz,
     .                   wsp,lwsp, iwsp,liwsp, cwsp,lcwsp, itrace,iflag)

      implicit none
      integer n, m, lwsp, liwsp, lcwsp, itrace, iflag, iwsp(liwsp),
     .     nz, nzmax
      parameter( nzmax=139301 )
      integer ia(nzmax), ja(nzmax)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp), a(nzmax)
      complex(kind=kind(0.0d0)) cwsp(lcwsp)

*-----Purpose----------------------------------------------------------|
*
*---  DMEXPV computes w = exp(t*A)*v - Customised for MARKOV CHAINS.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*     This is a customised version for Markov Chains. This means that a
*     check is done within this code to ensure that the resulting vector 
*     w is a probability vector, i.e., w must have all its components 
*     in [0,1], with sum equal to 1. This check is done at some expense
*     and the user may try DGEXPV which is cheaper since it ignores 
*     probability constraints.
*
*     IMPORTANT: The check assumes that the transition rate matrix Q
*                satisfies Qe = 0, where e=(1,...,1)'. Don't use DMEXPV
*                if this condition does not hold. Use DGEXPV instead.
*                DMEXPV/DGEXPV require the matrix-vector product 
*                y = A*x = Q'*x, i.e, the TRANSPOSE of Q times a vector.
*                Failure to remember this leads to wrong results.
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*                      
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested acurracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H     wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*              IMPORTANT: DMEXPV requires the product y = Ax = Q'x, i.e.
*              the TRANSPOSE of the transition rate matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = x_round, maximum among all roundoff errors (lower bound) 
*     wsp(4)  = s_round, sum of roundoff errors (lower bound)
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*     Markov chains are usually well-conditioned problems.
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 10000,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
* 
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hij, hump, SQR1,
     .                 roundoff, s_round, x_round

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2, DASUM

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) then
*         print *,'lwsp = ', lwsp
*         print *,'liwsp = ', liwsp
*         print *,'n = ', n
*         print *,'m = ', m
*         print *,'iflag = ', iflag
*         stop 'bad sizes (in input of DMEXPV)'
         return
      endif
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      sgn      = SIGN( 1.0d0,t )
      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      s_round  = 0.0d0
      x_error  = 0.0d0
      x_round  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

*     This is not initialised, causing compiler warnings.  For the time
C     where this is skipped, this should not matter (happy breakdown).
      avnorm   = 0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      call DCOPY( n, v,1, w,1 )
      beta = DNRM2( n, w,1 )
      vnorm = beta
      hump = beta
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )

      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call dgcoov( wsp(j1v-n), wsp(j1v), ia, ja, a, n, nz )
         do i = 1,j
            hij = DDOT( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
*            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call dgcoov( wsp(j1v-n), wsp(j1v), ia, ja, a, n, nz )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue

*
*---  compute w = beta*V*exp(t_step*H)*e1 ..
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         if ( iflag .lt. 0 ) then
*     This is an error that the calling function must recover from
            iflag = -42
            return
         endif
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         enddo
         wsp(iexph) = 1.0d0
         call DNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),cwsp(ifree+mx))
      endif

      continue
* 
*---  error estimate ...
* 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
*     Avoid compiler warning by referencing itrace (this is never
*     nonzero, but removing the argument would deviate from expokit
*     sources.
            itrace = 0
*            print*,'t_step =',t_old
*            print*,'err_loc =',err_loc
*            print*,'err_required =',delta*t_old*tol
*            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
c$$$            print*,"Failure in DMEXPV: ---"
c$$$            print*,"The requested tolerance is too high."
c$$$            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
      beta = DNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  Markov model constraints ...
*
      j = 0
      do i = 1,n
         if ( w(i).lt.0.0d0 ) then
            w(i) = 0.0d0
            j = j + 1
         endif
      enddo
      p1 = DASUM( n, w,1 )
      if ( j.gt.0 ) call DSCAL( n, 1.0d0/p1, w,1 )
      roundoff = DABS( 1.0d0-p1 ) / DBLE( n )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc, roundoff )
      err_loc = MAX( err_loc, rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
*      if ( itrace.ne.0 ) then
*         print*,'integration',nstep,'---------------------------------'
*         print*,'scale-square =',ns
*         print*,'wnorm     =',beta
*         print*,'step_size =',t_step
*         print*,'err_loc   =',err_loc
*         print*,'roundoff  =',roundoff
*         print*,'next_step =',t_new
*      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      s_round = s_round + roundoff
      x_error = MAX( x_error, err_loc )
      x_round = MAX( x_round, roundoff )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = x_round
      wsp(4)  = s_round
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      wsp(9)  = hump/vnorm
      wsp(10) = beta/vnorm

      END
