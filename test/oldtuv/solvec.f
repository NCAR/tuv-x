      SUBROUTINE SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL,
     &     MAZIM, NN, NSTR, YLM0, YLMC, CC, 
     &     EVECC, EVAL, KK, GC, AAD, EVECCD, EVALD,
     &     WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0, ZJ, ZZ,
     &     OPRIM, LC, DITHER, mu2, glsave, dgl)

*bm  SOLVEC calls SOLEIG and UPBEAM; if UPBEAM reports a potenially 
*bm  unstable solution, the calculation is repeated with a slightly 
*bm  changed single scattering albedo; this process is iterates 
*bm  until a stable solution is found; as stable solutions may be 
*bm  reached either by increasing or by decreasing the single 
*bm  scattering albedo, both directions are explored ('upward' and
*bm  'downward' iteration); the solution which required the smaller 
*bm  change in the single scattering albedo is finally returned 
*bm  by SOLVEC.

cgy added glsave and dgl to call to allow adjustable dimensioning

c     .. Scalar Arguments ..

      INTEGER, intent(in) :: MAZIM, NN, NSTR, LC
      REAL, intent(in)    :: DELM0, FBEAM, PI, UMU0, OPRIM, DITHER
      REAL, intent(in)    :: mu2

c     ..
c     .. Array Arguments ..

      INTEGER, intent(inout) ::   IPVT(:)
      
      REAL, intent(inout) :: WK(:)
      REAL, intent(in) ::  
     $     CMU(:), CWT(:), YLM0(0:), YLMC(0:,:)

      REAL, intent(inout) :: AMB(:,:), APB(:,:), ARRAY(:,:)
      REAL, intent(inout) :: EVAL(:), EVECC(:,:), KK(:), GC(:,:)
      REAL, intent(inout) :: GLSAVE(0:), DGL(0:), CC(:,:), ZJ(:), ZZ(:)
      REAL, intent(out)   :: GL(0:)

      INTEGER :: K
      REAL(8) :: AAD(:,:), EVALD(:), EVECCD(:,:), WKD(:)

*bm   Variables for instability fix
      
      INTEGER :: UAGAIN, DAGAIN
      REAL    :: MINRCOND, ADD, UADD, DADD, SSA, DSSA, FACTOR
      
      LOGICAL ::  DONE, NOUP, NODN, DEBUG, INSTAB
      
*bm   reset parameters

      DONE = .FALSE.
      NOUP = .FALSE.
      NODN = .FALSE.


*bm   flag for printing debugging output      
*      DEBUG  = .TRUE.
      DEBUG  = .FALSE.

*bm   instability parameter; the solution is considered 
*bm   unstable, if the RCOND reported by SGECO is smaller 
*bm   than MINRCOND
      MINRCOND = 5000. * R1MACH(4)

*bm   if an instability is detected, the single scattering albedo
*bm   is iterated downwards in steps of DADD and upwards in steps 
*bm   of UADD; in practice, MINRCOND and -MINRCOND should 
*bm   be reasonable choices for these parameters
      DADD    = -MINRCOND
      UADD    = MINRCOND

      UAGAIN = 0
      DAGAIN = 0
      ADD   = DADD
      

*bm   save array GL( ) because it will be 
*bm   changed if an iteration should be neccessary
      DO K = MAZIM, NSTR - 1
         GLSAVE( K ) =  GL( K )
      ENDDO
      
      SSA = OPRIM


*bm   in case of an instability reported by UPBEAM (INSTAB)
*bm   the single scattering albedo will be changed by a small 
*bm   amount (ADD); this is indicated by DAGAIN or UAGAIN 
*bm   being larger than 0; a change in the single scattering 
*bm   albedo is equivalent to scaling the array GL( )

 666  IF ( DAGAIN > 0 .OR. UAGAIN > 0)  THEN
         FACTOR = (SSA + ADD) / SSA
         DO K = MAZIM, NSTR - 1
            GL( K ) =  GL( K ) * FACTOR
         ENDDO

         SSA = SSA + ADD
         
*bm   if the single scattering albedo is now smaller than 0
*bm   the downward iteration is stopped and upward iteration 
*bm   is forced instead

         IF( SSA < DITHER) THEN
            NODN = .TRUE.
            DAGAIN = -1
            goto 778
         ENDIF

*bm   if the single scattering albedo is now larger than its maximum 
*bm   allowed value (1.0 - DITHER), the upward iteration is 
*bm   stopped and downward iteration is forced instead

         IF( SSA > rONE - DITHER) THEN
            NOUP = .TRUE.
            UAGAIN = -1
            goto 888
         ENDIF
      ENDIF


c     ** Solve eigenfunction problem in Eq. STWJ(8B);
c        return eigenvalues and eigenvectors

 777     CALL SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL,
     &     MAZIM, NN, NSTR, YLMC, CC, EVECC, EVAL,
     &     KK, GC, AAD, EVECCD, EVALD, WKD )

c     ** Calculate particular solutions of
c        q.SS(18) for incident beam source

      IF ( FBEAM > rZERO ) THEN
         CALL  UPBEAM( mu2,
     $        ARRAY, CC, CMU, DELM0, FBEAM, GL,
     $        IPVT, MAZIM, NN, NSTR, PI, UMU0, WK,
     $        YLM0, YLMC, ZJ, ZZ, MINRCOND, INSTAB)
      ENDIF
      
c     ** Calculate particular solutions of
c        Eq. SS(15) for thermal emission source
c        (not available in psndo.f)
      
*bm   finished if the result is stable on the first try
      IF ( (.NOT. INSTAB) .AND. 
     $     (UAGAIN == 0) .AND. (DAGAIN == 0)) THEN
         goto 999
      ENDIF

*bm   downward iteration
      IF( INSTAB .AND. UAGAIN == 0 )  THEN
         DAGAIN = DAGAIN + 1
         GOTO 666
      ENDIF
      
*bm   upward iteration
      IF( INSTAB .AND. UAGAIN > 0 )  THEN
         UAGAIN = UAGAIN + 1
         GOTO 666
      ENDIF


*bm   ( DAGAIN .NE. 0 ) at this place means that the downward
*bm   iteration is finished 

 778  IF (DAGAIN /= 0 .AND. UAGAIN == 0) THEN
         
*bm   save downward iteration data for later use and 
*bm   restore original input data
         DO K = MAZIM, NSTR - 1
            DGL( K ) =  GL( K )
            GL( K ) =  GLSAVE( K )
         ENDDO

         DSSA = SSA
         SSA = OPRIM

*bm   start upward iteration
         ADD = UADD
         UAGAIN = UAGAIN + 1
         GOTO 666
      ENDIF

*bm   both iterations finished
 888  IF (DONE) THEN
         goto 998
      ENDIF


*bm  if neither upward nor downward iteration converged, the 
*bm  original conditions are restored and SOLEIG/UPBEAM 
*bm  is called for the last time 
         
      IF (NOUP .AND. NODN) THEN
         
         DO K = MAZIM, NSTR - 1
            GL( K ) =  GLSAVE( K )
         ENDDO
         
         SSA = OPRIM
         
         IF (DEBUG) THEN
            write (*,*) '! *** Neither upward nor downward iteration'
            write (*,*) '! *** converged; using original result.'
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ENDIF

*bm  if upward iteration did not converge, the stable downward conditions
*bm  are restored and SOLEIG/UPBEAM is called for the last time
      IF (NOUP) THEN
         DO K = MAZIM, NSTR - 1
            GL( K ) =  DGL( K )
         ENDDO
         
         SSA = DSSA
         
         IF (DEBUG) THEN
            write (*,*) '! *** The upward iteration did not converge.'
            write (*,*) '! *** Had to iterate ', DAGAIN,
     $           ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ENDIF

*bm  if downward iteration did not converge, we are done 
*bm  (the result of the upward iteration will be used)
      IF (NODN) THEN
         IF (DEBUG) THEN
            write (*,*) '! *** The downward iteration did not converge.'
            write (*,*) '! *** Had to iterate ', UAGAIN,
     $           ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF
         
         DONE = .TRUE.
         GOTO 998
      ENDIF

      
*bm   if both iterations converged, and if the upward iteration 
*bm   required more steps than the downward iteration, the stable 
*bm   downward conditions are restored and SOLEIG/UPBEAM is 
*bm   called for the last time 
         
      IF (UAGAIN > DAGAIN) THEN
         DO K = MAZIM, NSTR - 1
            GL( K ) =  DGL( K )
         ENDDO
         
         SSA = DSSA
         
         IF (DEBUG) THEN
            write (*,*) '! *** Both iterations converged;',
     $           ' using downward.'
            write (*,*) '! *** Had to iterate ', DAGAIN,
     $        ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ELSE
         
         IF (DEBUG) THEN
            write (*,*) '! *** Both iterations converged;',
     $           ' using upward.'
            write (*,*) '! *** Had to iterate ', UAGAIN,
     $        ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         goto 998
      ENDIF
      
*bm   finally restore original input data
 998  DO K = MAZIM, NSTR - 1
         GL( K ) =  GLSAVE( K )
      ENDDO
      
 999  CONTINUE

      END SUBROUTINE SOLVEC
