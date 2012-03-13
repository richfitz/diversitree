*     A tidied version of the nucexp code, since it is getting a bit
*     ridiculous at the moment.

*     Generally, These are the functions to use.  Where possible, they
*     step incrementally through a series of times.  Sometimes it is not
*     possible to do this without getting into a tangle, in which case
*     they drop down and use the more basic approach of computing the
*     exponential for each clade (see below).
*       NUCEXP - return complete state matrix
*       NUCEXPL - return likelihoods for a series of clades
*
*     These functions are here mostly for testing; for each unique time
*     they compute the entire exponential.  These are much more robust
*     to errors.  The functions are analogous to NUCEXP and NUCEXPL,
*     taking identical arguments and retuning identical results to
*     rouding error.
*       NUCEXPSAFE
*       NUCEXPSAFEL
*
*     These ones are really old testing functions;
*       NUCEXP1
*       NUCEXP1L
*     they don't include the scal argument.
*
*     This won't need to be used directly:
*       NLDMAT: construct infinitesimal rate matrix

*     It might be nice if I can leverage a single mapping from the
*     output of a raw function (NUCEXP, NUCEXPG) to a likelihood one
*     (NUCEXPL, NUCEXPGL), but I don't want to keep around a huge pile
*     of extra results unnecessarily.  However, the intermediate
*     exponentiation *requires* that the entire matrix is returned, even
*     if we are only interested in a few elements.  Not having the
*     ability to call back makes it difficult to see a way around this.

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
*     Space is required as a function of the maximum number of species, nt:
*       state space: nt(nt+1)/2 + 1
*       non-zeros (bisse): (7nt^2 - 7nt + 2)/2
*       non-zeros (bisseness): (9nt^2 - 9nt - 2)/2
*     if we absorb at 200 species, these are 20101 and 139301 (bisse) or 179099 (bisseness),
*     respectively.

*     The likelihood subroutines (NUCEXP*L) compute the likelihood of
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
*       nt, mua, mub, laa, lab, qba, t, lt: as for BUCEXP
*       pac, paa, pbc, pba: for BiSSEness 
*       ti (int[lc]): vector of indices to 't'; the jth clade has time
*         t[ti[j]].  Conversely, the ith element of t is the time that
*         all clades, j, where ti[j] == i originate from.
*       Nc (int[lc]): Number of species in each clade (the ith clade has
*         Nc[i] species)
*       nsc (int[lc]): The number of species where the state is known
*         (0 < nsc[i] <= Nc[i])
*       k (int[lc]): The number of species known to be in state b
*         (0 <= k[i] <= nsc[i])
*       ans (double[4*lc]): Likelihoods, in four groups of lc elements:
*         1: Likelihood of generating clade, starting from state 0
*         2: Same, but starting from state 1
*         3: Likelihood of extinction, starting from state 0
*         4: Same, but starting from state 1
*     See nucexplik() for a description of how Nc, nsc and k are used.
*     One important case is that where Nc[i] = nsc[i]; in this case,
*     there are N species in the clade, and we have perfect information
*     about the distribution of the species in the clade.  In this case,
*     the likelihood is just one of the numbers in the probability
*     matrix.

***   Raw state space (NUCEXP/NUCEXPSAFE)
      subroutine NUCEXP(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, lt, scal,
     .     tol, m, w, iflag)
      implicit none

      integer nt, lt
      double precision mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t(lt), scal, w(*)

      double precision ninfnorm

      integer nmax, nzmax, mmax
      parameter( nmax=20101, nzmax=179099, mmax=30 )
      integer lwsp, liwsp
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax )

      integer d, i, n, nz, m, itrace, iflag
      integer ia(nzmax), ja(nzmax), iwsp(liwsp)
      double precision tol, anorm, v(nmax), wsp(lwsp), a(nzmax)

      itrace = 0
      iflag = 0

      n  = nt*(nt + 1)/2 + 1
*      nz = (7*nt*nt - 7*nt + 2)/2
      nz = (9*nt*nt - 9*nt - 2)/2

      call NLDMAT(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, ia, ja, a)

      anorm = ninfnorm(ia, ja, a, n, nz, wsp, lwsp)

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
*            print*,'[NUCEXP] WARNING: switching to manual calculation'
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
*                  print*,'[NUCEXP] WARNING: calculation failed'
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

      subroutine NUCEXPSAFE(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, lt,
     .     scal, tol, m, w, iflag)
      implicit none

      integer nt, lt
      double precision mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t(lt), scal, w(*)

      double precision ninfnorm

      integer nmax, nzmax, mmax
      parameter( nmax=20101, nzmax=179099, mmax=30 )
      integer lwsp, liwsp
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax )

      integer d, i, n, nz, m, itrace, iflag
      integer ia(nzmax), ja(nzmax), iwsp(liwsp)
      double precision tol, anorm, v(nmax), wsp(lwsp), a(nzmax)

      itrace = 0
      iflag = 0

      n  = nt*(nt + 1)/2 + 1
*      nz = (7*nt*nt - 7*nt + 2)/2
      nz = (9*nt*nt - 9*nt - 2)/2

      call NLDMAT(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, ia, ja, a)

      anorm = ninfnorm(ia, ja, a, n, nz, wsp, lwsp)

      v(1) = 0.0d0
      v(2) = scal
      do i = 3,n
         v(i) = 0.0d0
      enddo

      do d=0,1
         do i=1,lt
            call DSEXPV( n, m, t(i), v, w((lt*d + i-1)*n+1), tol,
     .           anorm,ia, ja, a, nz, wsp,lwsp, iwsp,liwsp, itrace,
     .           iflag, scal)
            if ( iflag .lt. 0 ) then
*               print*,'[NUCEXPSAFE] WARNING: calculation failed'
               return
            endif
         enddo

         if ( d .eq. 0 ) then
            v(2) = 0.0d0
            v(3) = scal
         endif
      enddo

      end

      subroutine NUCEXP1(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, lt, tol,
     .     m, w, iflag)
      implicit none

      integer nt, lt
      double precision mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t(lt), w(*)

      double precision ninfnorm

      integer nmax, nzmax, mmax
      parameter( nmax=20101, nzmax=179099, mmax=30 )
      integer lwsp, liwsp
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax )

      integer d, i, n, nz, m, itrace, iflag
      integer ia(nzmax), ja(nzmax), iwsp(liwsp)
      double precision tol, anorm, v(nmax), wsp(lwsp), a(nzmax)

      itrace = 0
      iflag = 0
      n  = nt*(nt + 1)/2 + 1
*      nz = (7*nt*nt - 7*nt + 2)/2
      nz = (9*nt*nt - 9*nt - 2)/2

      call NLDMAT(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, ia, ja, a)
      anorm = ninfnorm(ia, ja, a, n, nz, wsp, lwsp)

      v(1) = 0.0d0
      v(2) = 1.0d0
      do i = 3,n
         v(i) = 0.0d0
      enddo

      do d=0,1
         do i=1,lt
            call DMEXPV( n, m, t(i), v, w((lt*d + i-1)*n+1), tol, anorm,
     .           ia, ja, a, nz, wsp,lwsp, iwsp,liwsp, itrace, iflag)
            if ( iflag .lt. 0 ) then
*               print*,'[NUCEXP1] WARNING: calculation failed'
               return
            endif
         enddo
         if ( d .eq. 0 ) then
            v(2) = 0.0d0
            v(3) = 1.0d0
         endif
      enddo
      end

***   Likelihood functions
      subroutine NUCEXPL(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, lt,
     .     ti, Nc, nsc, k, lc, scal, tol, m, ans, iflag)
      implicit none

      integer nt, lt, lc, m, iflag
      integer ti(lc), Nc(lc), nsc(lc), k(lc)
      double precision mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t(lt), scal,
     .     tol, ans(4*lc)

      double precision NUCexplik

      integer n, d, i, j
      double precision w( 2 * lt * (nt*(nt + 1)/2 + 1) )

      iflag = 0
      n = nt*(nt + 1)/2 + 1

*     TODO: This is O(lt x lc), when it should be possible to do this in
*     O((lt+lc) * log(lt + lc)) (using an order()-style approach).
*     However, it's probably not that much of a time sink :)
      call NUCEXP(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, lt, scal, tol,
     .     m, w, iflag)
      if ( iflag .lt. 0 ) then
         print*,'[NUCEXPL] catching failure'
         return
      endif

      do d = 0,1
         do i=1,lt
            do j=1,lc
               if ( ti(j) .eq. i ) then
                  ans(d*lc + j) =
     .                 NUCexplik(Nc(j), nsc(j), k(j), 
     .                 w((lt*d + i-1)*n+1))
                  ans((2+d)*lc + j) = w((lt*d + i-1)*n+1)
               endif
            enddo
         enddo
      enddo

      end

      subroutine NUCEXPSAFEL(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, lt,
     .     ti, Nc, nsc, k, lc, scal, tol, m, ans, iflag)
      implicit none

      integer nt, lt, lc, m, iflag
      integer ti(lc), Nc(lc), nsc(lc), k(lc)
      double precision mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t(lt), scal,
     .     tol, ans(4*lc)

      double precision NUCexplik

      integer n, d, i, j
      double precision w( 2 * lt * (nt*(nt + 1)/2 + 1) )

      iflag = 0
      n = nt*(nt + 1)/2 + 1

      call NUCEXPSAFE(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, lt, scal,
     .     tol, m, w, iflag)
      if ( iflag .lt. 0 ) then
         print*,'[NUCEXPSAFE] catching failure'
         return
      endif
      do d = 0,1
         do i=1,lt
            do j=1,lc
               if ( ti(j) .eq. i ) then
                  ans(d*lc + j) =
     .                 NUCexplik(Nc(j), nsc(j), k(j), 
     .                 w((lt*d + i-1)*n+1))
                  ans((2+d)*lc + j) = w((lt*d + i-1)*n+1)
               endif
            enddo
         enddo
      enddo

      end

      subroutine NUCEXP1L(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, lt,
     .     ti, Nc, nsc, k, lc, tol, m, ans, iflag)
      implicit none

      integer nt, lt, lc, m, iflag
      integer ti(lc), Nc(lc), nsc(lc), k(lc)
      double precision mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t(lt), tol,
     .     ans(4*lc)

      double precision NUCexplik

      integer n, d, i, j
      double precision w( 2 * lt * (nt*(nt + 1)/2 + 1) )

      iflag = 0
      n  = nt*(nt + 1)/2 + 1

      call NUCEXP1(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, lt, tol, m, w, 
     .     iflag)
      if ( iflag .lt. 0 ) then
         print*,'[NUCEXP1L] catching failure'
         return
      endif
      do d = 0,1
         do i=1,lt
            do j=1,lc
               if ( ti(j) .eq. i ) then
                  ans(d*lc + j) =
     .                 NUCexplik(Nc(j), nsc(j), k(j), 
     .                 w((lt*d + i-1)*n+1))
                  ans((2+d)*lc + j) = w((lt*d + i-1)*n+1)
               endif
            enddo
         enddo
      enddo

      end

*     This should replace the inner loop so that we have
*      subroutine NUCEXP1L(nt, mua, mub, laa, lab, qba, qab, pac, paa,
*     .     pbc, pba, t, lt,
*     .     ti, Nc, nsc, k, lc, tol, m, ans)
*      implicit none
*
*      integer nt, lt, lc, m
*      integer ti(lc), Nc(lc), nsc(lc), k(lc)
*      double precision mua, mub, laa, lab, qba, qab, pac, paa,
*     .     pbc, pba, t(lt), tol,
*     .     ans(4*lc)
*
*      integer n
*      double precision w( 2 * lt * (nt*(nt + 1)/2 + 1) )
*
*      n  = nt*(nt + 1)/2 + 1
*      call NUCEXP1(nt, mua, mub, laa, lab, qba, qab, pac, paa,
*     .     pbc, pba, t, lt, tol, m, w)
*      call NUCEXPLIKELIHOODS(lt, ti, Nc, nsc, k, lc, w, n, ans)
*      end
      subroutine NUCEXPLIKELIHOODS(lt, ti, Nc, nsc, k, lc, w, n, ans)
      implicit none

      integer lt, lc, n
      integer ti(lc), Nc(lc), nsc(lc), k(lc)
      double precision w(2*lt*n), ans(4*lc)

      double precision NUCexplik
      integer d, i, j

      do d = 0,1
         do i=1,lt
            do j=1,lc
               if ( ti(j) .eq. i ) then
                  ans(d*lc + j) =
     .                 NUCexplik(Nc(j), nsc(j), k(j), 
     .                 w((lt*d + i-1)*n+1))
                  ans((2+d)*lc + j) = w((lt*d + i-1)*n+1)
               endif
            enddo
         enddo
      enddo
      
      end

***   Other
*     This is a historic function that is useful in debugging at the
*     Java end (since it is easy to make this return the same number as
*     the integrator does).
      subroutine NUCEXPONE(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, Nc, nsc,
     .     k, tol, m, ans, iflag)

      implicit none
      integer nt
      integer Nc(1), nsc(1), k(1), m, iflag
      double precision mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t(1), tol, ans(4)
      iflag = 0
*     the '1's are, in turn, lt, ti, [other args], lc
      call NUCEXP1L(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, t, 1, 1,
     .     Nc, nsc, k, 1, tol, m, ans, iflag)
      end

***     Helper functions and subroutines:
*     Infinite norm of a COO stored matrix:
      double precision function ninfnorm(ia, ja, a, n, nz, wsp, lwsp)
      implicit none
      integer nzmax
      integer n, nz, lwsp
      parameter( nzmax = 179099 )
*     SALLY:  Changed nzmax from 3000 to 179099
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax), wsp(lwsp)
      integer i

      do i = 1,n
         wsp(i) = 0.0d0
      enddo
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      enddo
      ninfnorm = wsp(1)
      do i = 2,n
         if ( ninfnorm.lt.DBLE(wsp(i)) ) ninfnorm =  wsp(i)
      enddo

      return
      end

*     Likelihood calculation for a particular state space/clade
*     composition:
*     We need to know:
*       N:  The number of species in the clade
*       ns: The number of species of known state (0 < ns <= N)
*       k:  The number of species known to be in state 'b' 
*           (0 <= k <= ns)
*       w:  The probability/state vector produced by NUCexp()
      double precision function NUCexplik(N, ns, k, w)
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
      
      NUCexplik = ans
      end

*     Rate matrix
*       nt: the number of species at which we absorb.
*       mua/b: extinction rates in state a and b
*       laa/b: speciation rates in state a and b
*       qba/qab: Character transition rates from a->b and b->a
*       t: time at which solution is needed

*     Speciation: Speciation to the current position is possible for
*     species in state a if there are more than one species (at least
*     two) currently in state a (if there is one species in state a,
*     then speciation cannot account for this one species - instead, it
*     must have arisen from extinction or character state change).
*     Species in state 'b' are considered before state 'a', since that
*     keeps the column indices in order.
*
*     Character state transition.  Transition from a->b is possible
*     wherever 'b' in this generation has at least one species,
*     irrespective of the number in state 'a'.  "Nothing" always happens
*     at rate na*pa + nb*pb
*
*     Extinction
*
*     Change species number counters
*
*     The last speciation row is dealt with separately, since this is
*     transition into the "too many species" absorbing state
      subroutine NLDMAT(nt, mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba, ia, ja, a)
      implicit none
      
      integer nt
      integer prodone, prodtwo
      double precision mua, mub, laa, lab, qba, qab, pac, paa,
     .     pbc, pba
      
      integer n, nz, idx, i, na, nb, ns, nzmax
      parameter( nzmax = 179099 )
*     SALLY:  Changed nzmax from 3000 to 179099
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      double precision pa, pb
      
      n = nt*(nt + 1)/2 + 1
*      nz = (7*nt*nt - 7*nt + 2)/2
      nz = (9*nt*nt - 9*nt - 2)/2

      na = 1
      nb = 0
      ns = 1

*     SALLY: Introduced for speed
      prodone = (ns + 1) * (ns + 2)/2
      prodtwo = ns * (ns - 1)/2

*     The chance of no change is unaffected by BiSSEness, which only alters what
*     happens at speciation.
      pa = -(mua + laa + qba)
      pb = -(mub + lab + qab)

*     SALLY: For all possible values of the index, idx, a(idx) gives the value of the matrix entry,
*     moving the system from ja(idx) to ia(idx). 
      a(1) = mua
      a(2) = mub
      ia(1) = 1
      ia(2) = 1
      ja(1) = 2
      ja(2) = 3

      idx = 3

      DO i = 2,(n-1)
*     BiSSE commands follow:
*         IF ( nb .gt. 1 ) THEN
*            a(idx) = (nb - 1) * lab
*            ia(idx) = i
*            ja(idx) = prodtwo + nb
*            idx = idx + 1
*         ENDIF
*         IF ( na .gt. 1 ) THEN
*            a(idx) = (na - 1) * laa
*            ia(idx) = i
*            ja(idx) = prodtwo + nb + 1
*            idx = idx + 1
*         ENDIF

*     BiSSEness commands follow:
         IF ( nb .gt. 1 ) THEN
            a(idx) = (nb - 1) * lab * (1 - pbc)
            IF ( na .gt. 0 ) THEN
               a(idx) = a(idx) + na * laa * pac * paa
            ENDIF
            ia(idx) = i
            ja(idx) = prodtwo + nb
            idx = idx + 1
            a(idx) = (na + 1) * laa * pac * (1 - paa)
            ia(idx) = i
            ja(idx) = prodtwo + nb - 1
            idx = idx + 1
         ENDIF
         IF ( nb .eq. 1 ) THEN
            IF ( na .gt. 0 ) THEN
               a(idx) = na * laa * pac * paa
            ia(idx) = i
            ja(idx) = prodtwo + nb
            idx = idx + 1
            ENDIF
         ENDIF
         IF ( na .gt. 1 ) THEN
            a(idx) = (na - 1) * laa * (1 - pac)
            IF ( nb .gt. 0 ) THEN
               a(idx) = a(idx) + nb * lab * pbc * pba
            ENDIF
            ia(idx) = i
            ja(idx) = prodtwo + nb + 1
            idx = idx + 1
            a(idx) = (nb + 1) * lab * pbc * (1 - pba)
            ia(idx) = i
            ja(idx) = prodtwo + nb + 2
            idx = idx + 1
         ENDIF
         IF ( na .eq. 1 ) THEN
            IF ( nb .gt. 0 ) THEN
               a(idx) = nb * lab * pbc * pba
            ia(idx) = i
            ja(idx) = prodtwo + nb + 1
            idx = idx + 1
            ENDIF
         ENDIF

*     If there are some species in state b, they could have come from state a:
*     SALLY:  Note that qab is actually q10 (it's the second of the two q elements)
*     Apparently, the notation used in the fortran code is that the last letter specifies the starting state.
         IF ( nb .gt. 0 ) THEN
            a(idx) = (na + 1) * qba
            ia(idx) = i
            ja(idx) = i - 1
            idx = idx + 1
         ENDIF

*     No change
         a(idx) = na * pa + nb * pb
         ia(idx) = i
         ja(idx) = i
         idx = idx + 1

*     If there are some species in state a, they could have come from state b:
         IF ( na .gt. 0 ) THEN
            a(idx) = (nb + 1) * qab
            ia(idx) = i
            ja(idx) = i + 1
            idx = idx + 1
         ENDIF
         
*     TODO: It is wasteful to compute (ns+1)(ns+2)/2 every time here
*     See above for similar cases, too [ns(ns-1)]
*     SALLY:  I fixed this a bit, defining prodone = (ns+1)(ns+2)/2 and prodtwo = ns(ns-1)/2
         IF ( ns < (nt - 1) ) THEN
            a(idx) = (na + 1) * mua
            ia(idx) = i
            ja(idx) = prodone + nb + 1
            idx = idx + 1

            a(idx) = (nb + 1) * mub
            ia(idx) = i
            ja(idx) = prodone + nb + 2
            idx = idx + 1
         ENDIF

         IF ( na .gt. 0 ) THEN
            na = na - 1
            nb = nb + 1
         ELSE
            ns = ns + 1
            prodone = (ns + 1) * (ns + 2)/2
            prodtwo = ns * (ns - 1)/2
            na = ns
            nb = 0
         ENDIF
      ENDDO

      na = nt - 1
      nb = 0
      DO i = (n - nt),(n - 1)
         a(idx) = na * laa + nb * lab
         ia(idx) = n
         ja(idx) = i
         idx = idx + 1
         na = na - 1
         nb = nb + 1
      ENDDO
      end
