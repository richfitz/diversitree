*     These are the functions to use.  Where possible, they step
*     incrementally through a series of times.  Sometimes it is not
*     possible to do this without getting into a tangle, in which case
*     they drop down and use the more basic approach of computing the
*     exponential for each clade (see below).
*       BUCEXP - return complete state matrix
*       BUCEXPL - return likelihoods for a series of clades
*
*     This won't need to be used directly:
*       BLDMAT: construct infinitesimal rate matrix

*     The complete state space subroutines compute the likelihood vector
*     for all of state space over a series of times.  Results stored in
*     'w', a 2*n*lt vector, where n=nt*(nt+1)/2+1 (i.e. each possible
*     state up to nt species) In 'w', each set of n elements is the
*     probability of getting any state up to nt species ((0,0), (1,0),
*     (0,1), (2, 0), ..., nt).  The ith group of n elements corresponds
*     to the ith time [t(i)].  The first lt sets of n elements
*     correspond to the initial condition (1,0), and the second lt sets
*     of n elements correspond to the initial condition (0,1).
*     
*     Space is required as a function of the number of species, nt:
*       state space: nt(nt+1)/2 + 1
*       non-zeros: (7nt^2 - 7nt + 2)/2
*     if we absorb at 200 species, these are 20101 and 139301,
*     respectively.

*     The likelihood subroutines (BUCEXP*L) compute the likelihood of
*     observing clades in certain states:
*
*     This takes computes the statistic for 'lc' clades that originate
*     from 'lt' nodes, and therefore have only 'lt' distinct times.
*     Therefore, we compute the matrix exponential only for these 'lt'
*     distinct times, and with the resulting probability matrix compute
*     the likelihood of observing the given clade for the one or two
*     clades that originate from that node.
*
*     Parameters:
*       nt, la0, la1, mu0, mu1, q01, t, lt: as for BUCEXP
*       ti (int[lc]): vector of indices to 't'; the jth clade has time
*         t[ti[j]].  Conversely, the ith element of t is the time that
*         all clades, j, where ti[j] == i originate from.
*       Nc (int[lc]): Number of species in each clade (the ith clade has
*         Nc[i] species)
*       nsc (int[lc]): The number of species where the state is known
*         (0 < nsc[i] <= Nc[i])
*       k (int[lc]): The number of species known to be in state 1
*         (0 <= k[i] <= nsc[i])
*       ans (double[4*lc]): Likelihoods, in four groups of lc elements:
*         1: Likelihood of generating clade, starting from state 0
*         2: Same, but starting from state 1
*         3: Likelihood of extinction, starting from state 0
*         4: Same, but starting from state 1
*     See bucexplik() for a description of how Nc, nsc and k are used.
*     One important case is that where Nc[i] = nsc[i]; in this case,
*     there are N species in the clade, and we have perfect information
*     about the distribution of the species in the clade.  In this case,
*     the likelihood is just one of the numbers in the probability
*     matrix.

***   Raw state space (BUCEXP)
      subroutine BUCEXP(nt, la0, la1, mu0, mu1, q01, q10, t, lt, scal,
     .     tol, m, w, iflag)
      implicit none

      integer nt, lt
      double precision la0, la1, mu0, mu1, q01, q10, t(lt), scal, w(*)

      double precision infnorm

      integer nmax, nzmax, mmax
      parameter( nmax=20101, nzmax=139301, mmax=30 )
      integer lwsp, liwsp
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax )

      integer d, i, n, nz, m, itrace, iflag
      integer ia(nzmax), ja(nzmax), iwsp(liwsp)
      double precision tol, anorm, v(nmax), wsp(lwsp), a(nzmax)

      itrace = 0
      iflag = 0

      n  = nt*(nt + 1)/2 + 1
      nz = (7*nt*nt - 7*nt + 2)/2

      call BLDMAT(nt, la0, la1, mu0, mu1, q01, q10, ia, ja, a)

      anorm = infnorm(ia, ja, a, n, nz, wsp, lwsp)

      v(1) = 0.0d0
      v(2) = scal
      do i = 3,n
         v(i) = 0.0d0
      enddo

      do d=0,1
         call DSEXPVI(n, m, t, lt, v, w(d * lt * n + 1), tol,
     .        anorm, ia, ja, a, nz, wsp,lwsp, iwsp,liwsp, itrace,
     .        iflag, scal)
         if ( iflag .lt. 0 ) then
*     The trick here would be to return the number of successful
*     times done and just restart from there, but that's not a very
*     large optimization for the number of times this is used (I counted
*     70 times we dropped down to this in 45000 evaluations [1/660
*     evaluations]
            do i=1,lt
               call DSEXPV( n, m, t(i), v, w((lt*d + i-1)*n+1), tol,
     .              anorm,ia, ja, a, nz, wsp,lwsp, iwsp,liwsp, itrace,
     .              iflag, scal)
               if ( iflag .lt. 0 ) then
                  return
               endif
            enddo
         endif

         if ( d .eq. 0 ) then
            v(2) = 0.0d0
            v(3) = scal
         endif
      enddo

      end

***   Likelihood functions
      subroutine BUCEXPL(nt, la0, la1, mu0, mu1, q01, q10, t, lt,
     .     ti, Nc, nsc, k, lc, scal, tol, m, ans, iflag)
      implicit none

      integer nt, lt, lc, m, iflag
      integer ti(lc), Nc(lc), nsc(lc), k(lc)
      double precision la0, la1, mu0, mu1, q01, q10, t(lt), scal,
     .     tol, ans(4*lc)

      double precision bucexplik

      integer n, d, i, j
      double precision w( 2 * lt * (nt*(nt + 1)/2 + 1) )

      iflag = 0
      n = nt*(nt + 1)/2 + 1

      call BUCEXP(nt, la0, la1, mu0, mu1, q01, q10, t, lt, scal, tol,
     .     m, w, iflag)
      if ( iflag .lt. 0 ) then
         return
      endif

      do d = 0,1
         do i=1,lt
            do j=1,lc
               if ( ti(j) .eq. i ) then
                  ans(d*lc + j) =
     .                 bucexplik(Nc(j), nsc(j), k(j), 
     .                 w((lt*d + i-1)*n+1))
                  ans((2+d)*lc + j) = w((lt*d + i-1)*n+1)
               endif
            enddo
         enddo
      enddo

      end

***     Helper functions and subroutines:
*     Infinite norm of a COO stored matrix:
      double precision function infnorm(ia, ja, a, n, nz, wsp, lwsp)
      implicit none
      integer nzmax
      integer n, nz, lwsp
      parameter( nzmax = 139301 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax), wsp(lwsp)
      integer i

*     Purely to avoid compiler warning:
      i = ja(1)

      do i = 1,n
         wsp(i) = 0.0d0
      enddo
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      enddo
      infnorm = wsp(1)
      do i = 2,n
         if ( infnorm.lt.DBLE(wsp(i)) ) infnorm =  wsp(i)
      enddo

      return
      end

*     Likelihood calculation for a particular state space/clade
*     composition:
*     We need to know:
*       N:  The number of species in the clade
*       ns: The number of species of known state (0 < ns <= N)
*       k:  The number of species known to be in state '1' 
*           (0 <= k <= ns)
*       w:  The probability/state vector produced by bucexp()
      double precision function bucexplik(N, ns, k, w)
      implicit none
      integer N, ns, k
      double precision w(*)

      integer i, off
      double precision ans, hyperg

      off = N*(N + 1)/2 + 1
      ans = 0.0d0

      if ( ns .eq. N ) then
         ans = w(off + k)
      elseif ( ns .eq. 1 ) then
         if ( k .eq. 0 ) then
            do i = 0,(N-1)
               ans = ans + (N - i + 0.0d0)/N * w(off + i)
            enddo
         else
            do i = 1,N
               ans = ans + (i + 0.0d0)/N * w(off + i)
            enddo
         endif
      elseif ( ns .eq. 0 ) then
         do i = 0,N
            ans = ans + w(off + i)
         enddo
      else
         do i = k,(N - ns + k)
            ans = ans + hyperg(N, i, ns, k) * w(off + i)
         enddo
      endif
      
      bucexplik = ans
      end

*     Rate matrix
*       nt: the number of species at which we absorb.
*       la0/1: speciation rates in state 0 and 1
*       mu0/1: extinction rates in state 0 and 1
*       q01/q10: Character transition rates from 0->1 and 1->0
*       t: time at which solution is needed

*     WARNING: comment here looks out of date...
*     Speciation: Speciation to the current position is possible for
*     species in state 0 if there are more than one species (at least
*     two) currently in state 0 (if there is one species in state 0,
*     then speciation cannot account for this one species - instead, it
*     must have arisen from extinction or character state change).
*     Species in state '1' are considered before state '0', since that
*     keeps the column indices in order.
*
*     Character state transition.  Transition from 0->1 is possible
*     wherever '0' in this generation has at least one species,
*     irrespective of the number in state '1'.  "Nothing" always happens
*     at rate n0*p0 + n1*p1
*
*     Extinction
*
*     Change species number counters
*
*     The last speciation row is dealt with separately, since this is
*     transition into the "too many species" absorbing state
      subroutine BLDMAT(nt, la0, la1, mu0, mu1, q01, q10, ia, ja, a)
      implicit none
      
      integer nt
      double precision la0, la1, mu0, mu1, q01, q10
      
      integer n, nz, idx, i, n0, n1, ns, nzmax
      parameter( nzmax = 3000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      double precision p0, p1
      
      n = nt*(nt + 1)/2 + 1
      nz = (7*nt*nt - 7*nt + 2)/2

      n0 = 1
      n1 = 0
      ns = 1

      p0 = -(mu0 + la0 + q01)
      p1 = -(mu1 + la1 + q10)

      a(1) = mu0
      a(2) = mu1
      ia(1) = 1
      ia(2) = 1
      ja(1) = 2
      ja(2) = 3

      idx = 3

      DO i = 2,(n-1)
         IF ( n1 .gt. 1 ) THEN
            a(idx) = (n1 - 1) * la1
            ia(idx) = i
            ja(idx) = ns*(ns - 1)/2 + n1
            idx = idx + 1
         ENDIF
         IF ( n0 .gt. 1 ) THEN
            a(idx) = (n0 - 1) * la0
            ia(idx) = i
            ja(idx) = ns*(ns - 1)/2 + n1 + 1
            idx = idx + 1
         ENDIF

         IF ( n1 .gt. 0 ) THEN
            a(idx) = (n0 + 1) * q01
            ia(idx) = i
            ja(idx) = i - 1
            idx = idx + 1
         ENDIF

         a(idx) = n0 * p0 + n1 * p1
         ia(idx) = i
         ja(idx) = i
         idx = idx + 1

         IF ( n0 .gt. 0 ) THEN
            a(idx) = (n1 + 1) * q10
            ia(idx) = i
            ja(idx) = i + 1
            idx = idx + 1
         ENDIF
         
*     TODO: It is wasteful to compute (ns+1)(ns+2)/2 every time here
*     See above for similar cases, too [ns(ns-1)]
         IF ( ns < (nt - 1) ) THEN
            a(idx) = (n0 + 1) * mu0
            ia(idx) = i
            ja(idx) = (ns + 1)*(ns + 2)/2 + n1 + 1
            idx = idx + 1

            a(idx) = (n1 + 1) * mu1
            ia(idx) = i
            ja(idx) = (ns + 1)*(ns + 2)/2 + n1 + 2
            idx = idx + 1
         ENDIF

         IF ( n0 .gt. 0 ) THEN
            n0 = n0 - 1
            n1 = n1 + 1
         ELSE
            ns = ns + 1
            n0 = ns
            n1 = 0
         ENDIF
      ENDDO

      n0 = nt - 1
      n1 = 0
      DO i = (n - nt),(n - 1)
         a(idx) = n0 * la0 + n1 * la1
         ia(idx) = n
         ja(idx) = i
         idx = idx + 1
         n0 = n0 - 1
         n1 = n1 + 1
      ENDDO
      end
