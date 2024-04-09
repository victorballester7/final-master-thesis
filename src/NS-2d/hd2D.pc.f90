!=================================================================
PROGRAM HD2D
!=================================================================
! HD2D code
!
! Numerically integrates the incompressible HD equations
! in 2 dimensions using the streamfunction formulation
! with an external force.
! A pseudo-spectral method is used to compute spatial
! derivatives, while variable order Runge-Kutta method
! is used to evolve the system in time domain.
! To compile, you need the FFTW library installed on
! your system. You should link with the FFTP subroutines
! and use the FFTPLANS and MPIVARS modules (see the file
! 'fftp_mod.f90').
!
! NOTATION: index 'i' is 'x'
!           index 'j' is 'y'
!=================================================================

!      USE mpivars
   USE fft
   USE ali
   USE var
   USE kes
   USE grid
   USE random
   IMPLICIT NONE

!
! Integration parameters
!     ord  : order of the Runge-Kutta method used

   INTEGER, PARAMETER :: ord = 4
   INTEGER :: ini
   INTEGER :: step
   INTEGER :: tstep
   INTEGER :: cstep
   INTEGER :: sstep

!
! streamfunction, vector potential, z component
! of the fields and external force matrixes

   DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: ps

   DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: fk

!
! Temporal data storing matrixes

   DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: C1
   DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: C2
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: R1
!
! Some auxiliary matrixes

   DOUBLE PRECISION :: ener, enerk
   DOUBLE PRECISION :: enst
   DOUBLE PRECISION :: dt,dt_new,cfl
   DOUBLE PRECISION :: kup,kdn
   DOUBLE PRECISION :: prm1,prm2
   DOUBLE PRECISION :: dump,tmp,tmp1
   DOUBLE PRECISION :: f0,u0
   DOUBLE PRECISION :: time
   DOUBLE PRECISION :: nu,hnu

   INTEGER :: mult,stat
   INTEGER :: t,o
   INTEGER :: i,j
   INTEGER :: ki
   INTEGER :: ic,id,iu
   INTEGER :: jc,jd,ju
   INTEGER :: timet,timec,times

   ! my variables
   INTEGER :: threshold, seed, iflow


   CHARACTER     :: c,d,u
   CHARACTER*3   :: node,ext
   CHARACTER*100 :: ldir
   CHARACTER*100 :: dir_data


!
   node = '000'

!
! Allocates memory

   ALLOCATE( R1(n,n) )
   ALLOCATE( C1(n/2+1,n) )
   ALLOCATE( C2(n/2+1,n) )
   ALLOCATE( ps(n/2+1,n) )
   ALLOCATE( fk(n/2+1,n) )
   ALLOCATE( ka(n), ka2(n/2+1,n) )

!
! Reads from the external file 'status.prm'
! the status of a previous run (if any)
!     stat: last output of a previous run
!     mult: time step multiplier

   OPEN(1,file='./status.prm',status='unknown')
   READ(1,*) stat
   READ(1,*) mult
   READ(1,*) time
   CLOSE(1)
!
! Reads from the external file 'parameter.txt' the
! parameters that will be used during the integration
!     dt   : time step size
!     step : total number of time steps to compute
!     tstep: number of steps between I/O writing
!     sstep: number of steps between power spectrum I/O
!     cstep: number of steps between information output
!     f0   : amplitude of the external kinetic force
!     u0   : amplitude of the initial streamfunction
!     kdn  : minimum wave number in the external force
!     kup  : maximum wave number in the external force
!     nu   : kinematic viscosity
!     hmu  : hypo viscosity
!     prm1 : free parameter 1
!     prm2 : free parameter 2
!     ldir : local directory for I/O
!     dir_data : global directory for I/O

   OPEN(1,file='./input.prm',status='unknown')
   READ(1,*) cfl         !  1
   READ(1,*) step        !  2
   READ(1,*) tstep       !  3
   READ(1,*) sstep       !  4
   READ(1,*) cstep       !  5
   READ(1,*) f0          !  6
   READ(1,*) u0          !  7
   READ(1,*) kdn         !  8
   READ(1,*) kup         !  9
   READ(1,*) nu          ! 10
   READ(1,*) hnu         ! 11
   READ(1,*) prm1        ! 12
   READ(1,*) prm2        ! 13
   READ(1,*) seed       ! 14
   READ(1,*) iflow       ! 15
   READ(1,'(a100)') ldir ! 16
   READ(1,'(a100)') dir_data
   CLOSE(1)
   cfl = cfl/dble(mult)
   step = step*mult
   tstep = tstep*mult
   sstep = sstep*mult
   cstep = cstep*mult

   print*, "cfl    =",cfl    !  1
   print*, "step   =",step   !  2
   print*, "tstep  =",tstep  !  3
   print*, "sstep  =",sstep  !  4
   print*, "cstep  =",cstep  !  5
   print*, "f0     =",f0     !  6
   print*, "u0     =",u0     !  7
   print*, "kdn    =",kdn    !  8
   print*, "kup    =",kup    !  9
   print*, "nu     =",nu     ! 10
   print*, "hnu    =",hnu    ! 11
   print*, "prm1   =",prm1   ! 12
   print*, "prm2   =",prm2   ! 13
   print*, "ldir   =",ldir   ! 14
   print*, "dir_data=",dir_data
!
! Some numerical constants

   ic = 48
   id = 48
   iu = 48
   jc = 48
   jd = 48
   ju = 48

!
! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

   kmax = (dble(n)/3.d0)**2
   tiny =  0.000001d0
   ener = 1.0d0
   enst =  1.0d0

!
! Builds the wave number and the square wave
! number matrixes

   DO i = 1,n/2
      ka(i) = dble(i-1)
      ka(i+n/2) = dble(i-n/2-1)
   END DO
   DO j = 1,n
      DO i = 1,n/2+1
         ka2(i,j) = ka(i)**2+ka(j)**2
      END DO
   END DO

!
! Initializes the FFT library
! Use FFTW_ESTIMATE in short runs and FFTW_MEASURE
! in long runs

   call rfftw2d_f77_create_plan(planrc,n,n,FFTW_REAL_TO_COMPLEX, &
      FFTW_ESTIMATE)
   call rfftw2d_f77_create_plan(plancr,n,n,FFTW_COMPLEX_TO_REAL, &
      FFTW_ESTIMATE)



!FORCING
   CALL forcing(iflow,f0,kup,kdn,seed,fk)

! INITIAL CONDITIONS
   IF (stat.eq.0) THEN
      ini = 1
      timet = tstep
      timec = cstep
      times = sstep

!STREAM FUNCTION R1 & VECTOR POTENTIAL R2
      DO j = 1,n
         DO i = 1,n
            R1(i,j) = 0.0d0
         END DO
      END DO
      DO ki=1,5
         DO j = 1,n
            DO i = 1,n
               !  IF ( (i.ge.n/2-threshold).and.(i.le.n/2+threshold) .and. (j.ge.n/2-threshold).and.(j.le.n/2+threshold) ) THEN
               !     R1(i,j) = 1.0d0
               !  ELSE
               !     R1(i,j) = 0.0d0
               !  ENDIF
               R1(i,j) = R1(i,j)                                     &
                  + 1.0d0*sin(2.0d0*ki *pi*(dble(i)-1)/dble(n)) &
                  + 1.0d0*sin(2.0d0*ki *pi*(dble(j)-1)/dble(n)) &
                  + 1.0d0*sin(2.0d0*ki *pi*(dble(i)-1)/dble(n)) &
                  * 1.0d0*sin(2.0d0*ki *pi*(dble(j)-1)/dble(n))
            END DO
         END DO
      ENDDO
      CALL rfftwnd_f77_one_real_to_complex(planrc, R1, ps)
      CALL energy(ps,ener,1)
      tmp1=u0/sqrt(ener)
      DO j = 1,n
         DO i = 1,n/2+1
            IF ((ka2(i,j).le.kmax).and.(ka2(i,j).ge.tiny)) THEN
               ps(i,j) = tmp1*ps(i,j)
            ELSE
               ps(i,j) = 0.0d0
            ENDIF
         END DO
      END DO
   ELSE
      print*,'READING...',stat
      ini = int((stat-1)*tstep)
      dump = dble(ini)/dble(sstep)+1
      times = 0
      timet = 0
      timec = 0
      jc = 48+int(dump/100)
      jd = 48+int(dump/10)-int(dump/100)*10
      ju = 48+int(dump)-int(dump/10)*10
      ic = 48+int(float(stat)/100)
      id = 48+int(float(stat)/10)-int(float(stat)/100)*10
      iu = 48+int(stat)-int(float(stat)/10)*10
      c = char(ic)
      d = char(id)
      u = char(iu)

      OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
         // c // d // u //'.out',form='unformatted')
      READ(1) R1
      CALL rfftwnd_f77_one_real_to_complex(planrc, R1, ps)
   ENDIF



!
! Time integration scheme starts here
! Uses Runge-Kutta of order 'ord'
!#################### MAIN LOOP ######################
   RK : DO t = ini,step

!!!!!!!  RANDOM FORCING  !!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL forcing(iflow,f0,kup,kdn,seed,fk)  !! set fk=0
      ! CALL energy(fk,enerk,1)
      ! tmp1=f0/sqrt(0.5*enerk*dt)*randu(seed)
      ! fk  = tmp1*fk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Every 'cstep' steps, generates external files
! to check consistency and convergency. See the
! mhdcheck subroutine for details.


      dt_new = 0.0000002d0
      IF (timec.eq.cstep) THEN
         timec = 0
         CALL hdcheck(ps,fk,time,ener,enst,dir_data)
         if ( dt_new.ge. cfl/(dble(n)*sqrt(ener)) ) THEN
            dt_new = cfl/(dble(n)*sqrt(ener))
         endif
         if ( dt_new .le. 0.5*cfl/(dble(n)*sqrt(ener)) ) THEN
            dt_new = 0.5d0*cfl/(dble(n)*sqrt(ener))
         endif
         if ( dt_new .ge. 0.5d0 ) THEN
            dt_new = 0.5d0
         endif
         dt=dt_new/mult
         ! if it's the first time, print header
         IF (t.eq.ini) THEN
            print*,"DBG", "           t", "   dt", "                        time","                     energy","        enstrophy"
         ENDIF
         print*,"DBG",t,dt,time,ener,enst
      ENDIF

! Every 'sstep' steps, generates external files
! with the power spectrum

      IF (times.eq.sstep) THEN
         times = 0
         ju = ju+1
         IF (ju.eq.58) THEN
            ju = 48
            jd = jd+1
         ENDIF
         IF (jd.eq.58) THEN
            jd = 48
            jc = jc+1
         ENDIF
         c = char(jc)
         d = char(jd)
         u = char(ju)
         ext = c // d // u
         CALL spectrum(ps,ext,1,dir_data)
         CALL Eprof(ps,ext,dir_data)
         CALL laplak2(ps,C1)     ! make W

         !  CALL ring_vorticity(C1,ext,dir_data)

         CALL vectrans(ps,ps,C1,'euu',ext,dir_data)
         CALL vectrans(C1,ps,C1,'vrt',ext,dir_data)
      ENDIF

! Every 'tstep' steps, stores the results of the integration

      IF (timet.eq.tstep) THEN
         timet = 0
         iu = iu+1
         IF (iu.eq.58) THEN
            iu = 48
            id = id+1
         ENDIF
         IF (id.eq.58) THEN
            id = 48
            ic = ic+1
         ENDIF
         c = char(ic)
         d = char(id)
         u = char(iu)

         DO j = 1,n
            DO i = 1,n/2+1
               C1(i,j) = ps(i,j)/dble(n)**2 ! normalizes for the FFT
            END DO
         END DO
         CALL rfftwnd_f77_one_complex_to_real(plancr, C1, R1)
         OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
            // c // d // u // '.out',form='unformatted')
         WRITE(1) R1
         CLOSE(1)

         DO j = 1,n
            DO i = 1,n/2+1
               C1(i,j) = ps(i,j)*ka2(i,j)/dble(n)**2
            END DO
         END DO
         CALL rfftwnd_f77_one_complex_to_real(plancr, C1, R1)
         OPEN(1,file=trim(ldir) // '/hd2Dww.' // node // '.' &
            // c // d // u // '.out',form='unformatted')
         WRITE(1) R1
         CLOSE(1)

      ENDIF

      timet = timet+1
      times = times+1
      timec = timec+1
      time = time+dt


! Runge-Kutta step 1
! Copies the streamfunction into the auxiliary matrix C1

      DO j = 1,n
         DO i = 1,n/2+1
            C1(i,j) = ps(i,j)
         END DO
      END DO

! Runge-Kutta step 2

      DO o = ord,1,-1

         CALL laplak2(C1,C2)     ! make W
         CALL poisson(C1,C2,C1)  ! u grad w


         DO j = 1,n
            DO i = 1,n/2+1
               IF ((ka2(i,j).le.kmax).and.(ka2(i,j).ge.tiny)) THEN
                  C1(i,j) = (ps(i,j)+ dt*( - C1(i,j)/ka2(i,j)  &
                     + fk(i,j))/dble(o)) &
                     /(1.0d0 + (hnu/ka2(i,j)**2+nu*ka2(i,j))*dt/dble(o))
               ELSE
                  C1(i,j) = 0.0d0
               ENDIF

            END DO
         END DO

      END DO

! Runge-Kutta step 3
! Copies the result from the auxiliary matrixes into ps, az


      DO j = 1,n
         DO i = 1,n/2+1
            IF ((ka2(i,j).le.kmax).and.(ka2(i,j).ge.tiny)) THEN
               ps(i,j) = C1(i,j)
            ELSE
               ps(i,j) = 0.0d0
            ENDIF
         END DO
      END DO

   END DO RK
!##############  END OF MAIN LOOP ###################

!
! End of Runge-Kutta

   CALL rfftwnd_f77_destroy_plan(plancr)
   CALL rfftwnd_f77_destroy_plan(planrc)
   DEALLOCATE( R1 )
   DEALLOCATE( ps,fk )
   DEALLOCATE( C1,C2 )
   DEALLOCATE( ka,ka2)

END PROGRAM HD2D
