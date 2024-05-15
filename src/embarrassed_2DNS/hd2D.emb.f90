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

   USE mpivars
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
   DOUBLE PRECISION :: Re_nu
   DOUBLE PRECISION :: nu,mu
   DOUBLE PRECISION :: phase1,phase2,phase3,phase4

   INTEGER :: mult,stat
   INTEGER :: t,o
   INTEGER :: i,j
   INTEGER :: ki
   INTEGER :: ic,id,iu
   INTEGER :: jc,jd,ju
   INTEGER :: timet,timec,times
   INTEGER :: inu,imu

   ! my variables
   INTEGER :: threshold, seed, iflow,myseed
   INTEGER :: start_time, end_time, rate


   CHARACTER     :: c,d,u
   CHARACTER*3   :: node,ext
   CHARACTER*100 :: ldir



!
! Initializes the MPI library
   CALL MPI_INIT(ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
   ic = 48+int(myrank/100) ! id for the node (cent)
   id = 48+int(myrank/10)-int(myrank/100)*10 ! id for the node (dec)
   iu = 48+int(myrank)-int(myrank/10)*10 ! id for the node (unit)
   c = char(ic)
   d = char(id)
   u = char(iu)
   node = c // d // u

   ! count execution time
   IF (myrank.eq.0) CALL system_clock(start_time, rate)


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

   IF (myrank.eq.0) THEN
      OPEN(1,file='./status.prm',status='unknown')
      READ(1,*) stat
      READ(1,*) time
      CLOSE(1)
   ENDIF
   CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

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
!     inu  : power of the squared wave number in the viscosity, e.g. nu* laplacian^inu
!     mu   : hipo viscosity
!     imu  : power of the squared wave number in the hyper viscosity, e.g. mu* laplacian^(-imu)
!     prm1 : free parameter 1
!     prm2 : free parameter 2
!     ldir : local directory for I/O

   IF (myrank.eq.0) THEN
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
      READ(1,*) Re_nu                       ! 10
      READ(1,*) inu                      ! 11
      READ(1,*) mu                       ! 12
      READ(1,*) imu                      ! 13
      READ(1,*) iflow                    ! 14
      READ(1,*) seed                     ! 15
      READ(1,*) prm1                     ! 16
      READ(1,*) prm2                     ! 17
      READ(1,'(a100)') ldir              ! 18
      CLOSE(1)

      OPEN(1,file=trim(ldir) // '/dim.txt')
      WRITE(1,*) n, nprocs
      CLOSE(1)

      ! formula for the Reynolds number associated with the kinematic viscosity
      nu = (f0**2 * 4 * kdn**2/pi)**(1.0/3.0) * kup**(-4.0/3.0) / Re_nu

      print*, "dim    =",n      !  0
      print*, "nprocs =",nprocs !  0
      print*, "cfl    =",cfl    !  1
      print*, "step   =",step   !  2
      print*, "tstep  =",tstep  !  3
      print*, "sstep  =",sstep  !  4
      print*, "cstep  =",cstep  !  5
      print*, "f0     =",f0     !  6
      print*, "u0     =",u0     !  7
      print*, "kdn    =",kdn    !  8
      print*, "kup    =",kup    !  9
      print*, "Re_nu  =",Re_nu     ! 10
      print*, "nu     =",nu     ! 10
      print*, "inu    =",inu    ! 11
      print*, "mu     =",mu    ! 12
      print*, "imu    =",imu    ! 13
      print*, "iflow  =",iflow  ! 14
      print*, "seed   =",seed   ! 15
      print*, "prm1   =",prm1   ! 16
      print*, "prm2   =",prm2   ! 17
      print*, "ldir   =",trim(ldir)   ! 18
   ENDIF

   CALL MPI_BCAST(  CFL,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 1
   CALL MPI_BCAST( step,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 2
   CALL MPI_BCAST(tstep,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 3
   CALL MPI_BCAST(sstep,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 4
   CALL MPI_BCAST(cstep,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 5
   CALL MPI_BCAST(   f0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 6
   CALL MPI_BCAST(   u0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 7
   CALL MPI_BCAST(  kdn,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 8
   CALL MPI_BCAST(  kup,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 9
   CALL MPI_BCAST(   nu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !10
   CALL MPI_BCAST(  inu,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) !11
   CALL MPI_BCAST(   mu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !12
   CALL MPI_BCAST(  imu,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) !13
   CALL MPI_BCAST(iflow,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) !14
   CALL MPI_BCAST( seed,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) !15
   CALL MPI_BCAST( prm1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !16
   CALL MPI_BCAST( prm2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !17
   CALL MPI_BCAST( ldir,100,MPI_CHARACTER,     0,MPI_COMM_WORLD,ierr) !19
   myseed=seed+myrank


!
! Some numerical constants

   ic = 48
   id = 48
   iu = 48
   jc = 48
   jd = 48
   ju = 48-1

!
! Some constants for the FFT
!     kmax2: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

   kmax2 = (dble(n)/3.d0)**2
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


! Sets the initial conditions.

! INITIAL CONDITIONS
   IF (stat.eq.0) THEN ! we start from scratch
      ini = 1
      timet = 0
      timec = cstep
      times = sstep

!STREAM FUNCTION R1
      CALL initialcond(iflow,u0,kup,kdn,seed,myseed,ps)
      CALL energy(ps,ener,1)
      tmp1=u0/sqrt(ener)
      DO j = 1,n
         DO i = 1,n/2+1
            IF ((ka2(i,j).le.kmax2).and.(ka2(i,j).ge.tiny)) THEN
               ps(i,j) = tmp1*ps(i,j)
            ELSE
               ps(i,j) = 0.0d0
            ENDIF
         END DO
      END DO
   ELSE ! we start from a previous run
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

      OPEN(1,file=trim(ldir) // '/output/hd2Dps.' // node // '.' &
         // c // d // u //'.out',form='unformatted')
      READ(1) R1
      CLOSE(1)
      print*,"READING DONE!"
      CALL rfftwnd_f77_one_real_to_complex(planrc, R1, ps)
      print*,"READING & FFT DONE!"
   ENDIF



   if (myrank.eq.0) print*,"START RK",step
!
! Time integration scheme starts here
! Uses Runge-Kutta of order 'ord'
!#################### MAIN LOOP ######################
   ini = 1
   RK : DO t = ini,step
      ! update the time step
      ! CALL CFL_condition(CFL,ps,inu,nu,dt)

      ! keep dt constant
      dt = 2.0d-5

!!!!!!!  RANDOM FORCING  !!!!!!!!!!!
      CALL forcing(iflow,f0,kup,kdn,seed,myseed,fk)
      CALL energy(fk,enerk,1)
      !IF (enerk.le.tiny) THEN ! no forcing
      !   tmp1=1.0d0
      !ELSE
      tmp1=f0/sqrt(0.5*enerk*dt) ! we normalize the energy injection rate, not the forcing amplitude
      !ENDIF
      fk  = tmp1*fk
      ! DO j = 1,n
      !    DO i = 1,n/2+1
      !       IF ((ka2(i,j).le.kmax2).and.(ka2(i,j).ge.tiny)) THEN
      !          fk(i,j) = tmp1*fk(i,j)
      !       ELSE
      !          fk(i,j) = 0.0d0
      !       ENDIF
      !    END DO
      ! END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Every 'cstep' steps, generates external files
! to check consistency and convergency. See the
! hdcheck subroutine for details.

      IF (timec.eq.cstep) THEN
         timec = 0
         CALL hdcheck(ps,fk,time,inu,nu,imu,mu,ener,enst,node,ldir)
         ! if it's the first time, print header
         IF (t.eq.ini .AND. myrank.eq.0) THEN
            print*,"DBG", "           t", "   dt", "                        time","                     energy","        enstrophy"
         ENDIF
         if (myrank.eq.0) print*,"DBG",t,dt,time,ener,enst
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
         CALL spectrum(ps,ext,node,ldir)
         CALL EnergyEnstropy_profiles(ps,ext,node,ldir)
         CALL laplak2(ps,C1)     ! make W
         CALL vectrans(ps,ps,C1,'euu',ext,node,ldir)
         CALL vectrans(C1,ps,C1,'vrt',ext,node,ldir)
         IF (myrank.eq.0) THEN
            OPEN(1,file=trim(ldir)//'/spectra_times.txt',position='append')
            WRITE(1,13) ext,time
13          FORMAT( A3,    F12.6)
            CLOSE(1)
         ENDIF
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
         ext = c // d // u
         CALL outputfields(ps,fk,ext,node,ldir)

         IF (myrank.eq.0) THEN
            OPEN(1,file=trim(ldir)//'/field_times.txt',position='append')
            WRITE(1,12) c//d//u,time
12          FORMAT( A3,    F12.6)
            CLOSE(1)
            OPEN(1,file='status.prm')
            WRITE(1,*) c//d//u,'          % stat'
            WRITE(1,*) time,'  % time'
            CLOSE(1)
         ENDIF
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
               IF ((ka2(i,j).le.kmax2).and.(ka2(i,j).ge.tiny)) THEN
                  ! this reformulation can be deduced from the general formula of y'=Ay+b: y_{n+1}=e^(A*dt)[y_n + int_0^dt e^(A(-s))b(s)ds]
                  ! explicit Euler is equivalent to approximating the latter integral by e^(A*0)b(0)dt= b(0)dt
                  ! ==> y_{n+1}=e^(A*dt)[y_n + b(y_n)dt]
                  C1(i,j) = (ps(i,j)+ dt*( - C1(i,j)/ka2(i,j)  &
                     + fk(i,j))/dble(o))
                  tmp    = exp(-(mu/ka2(i,j)**imu+nu*ka2(i,j)**inu)*dt/dble(o))
                  !! tmp     = 1.0D0 /(1.0d0 + (-mu/ka2(j,i)**imu+nu*ka2(j,i)**inu)*dt/dble(o))
                  C1(i,j) = C1(i,j) *tmp
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
            IF ((ka2(i,j).le.kmax2).and.(ka2(i,j).ge.tiny)) THEN
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
   CALL MPI_FINALIZE(ierr)
   CALL rfftwnd_f77_destroy_plan(plancr)
   CALL rfftwnd_f77_destroy_plan(planrc)
   DEALLOCATE( R1 )
   DEALLOCATE( ps,fk )
   DEALLOCATE( C1,C2 )
   DEALLOCATE( ka,ka2)

   IF (myrank.eq.0) THEN
      CALL system_clock(end_time)
      print*, "Execution time (seconds): ", (end_time - start_time)/real(rate)
   ENDIF

END PROGRAM HD2D
