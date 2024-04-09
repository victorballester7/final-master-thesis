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
!
! 2004 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
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
   DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: f1

!
! Temporal data storing matrixes

   DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:) :: C1
   DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:) :: C2
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: R1
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: R2
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: xv,yv
   INTEGER, ALLOCATABLE, DIMENSION (:)    :: fractal
!
! Some auxiliary matrixes

   DOUBLE PRECISION :: ener
   DOUBLE PRECISION :: enerk,enst
   DOUBLE PRECISION :: dt,dt_new,CFL
   DOUBLE PRECISION :: kup,kdn,kr
   DOUBLE PRECISION :: prm1,prm2
   DOUBLE PRECISION :: dump,tmp
   DOUBLE PRECISION :: tmp1,tmp2,tmp3,tmp4,tmp5
   DOUBLE PRECISION :: f0,u0
   DOUBLE PRECISION :: time
   DOUBLE PRECISION :: nu,mu
   DOUBLE PRECISION :: phase1,phase2
   DOUBLE PRECISION :: phase3,phase4
   DOUBLE PRECISION :: advfc,xfpos,yfpos

   INTEGER :: stat
   INTEGER :: t,o
   INTEGER :: i,j,ir,jr
   INTEGER :: ki,kj,jj,ii
   INTEGER :: ic,id,iu
   INTEGER :: jc,jd,ju,jt
   INTEGER :: timet,timec,times
   INTEGER :: seed,seed1,myseed
   INTEGER :: iflow,inu,imu,itmp,jtmp,nfat,nthn,ifat,ithn

   CHARACTER     :: c,d,u,th
   CHARACTER*3   :: node,ext
   CHARACTER*4   :: ext4
   CHARACTER*100 :: ldir

!
! Initializes the MPI library

   CALL MPI_INIT(ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
   ic = 48+int(myrank/100)
   id = 48+int(myrank/10)-int(myrank/100)*10
   iu = 48+int(myrank)-int(myrank/10)*10
   c = char(ic)
   d = char(id)
   u = char(iu)
   node = c // d // u

!
! Allocates memory for distributed blocks

   CALL range(1,n/2+1,nprocs,myrank,ista,iend)
   CALL range(1,n,nprocs,myrank,jsta,jend)

   ALLOCATE( R1(n,jsta:jend) )
   ALLOCATE( R2(n,jsta:jend) )
   ALLOCATE( C1(n,ista:iend) )
   ALLOCATE( C2(n,ista:iend) )
   ALLOCATE( ps(n,ista:iend) )
   ALLOCATE( fk(n,ista:iend) )
   ALLOCATE( f1(n,ista:iend) )
   ALLOCATE( ka(n), ka2(n,ista:iend) )
   ALLOCATE( xv(8),yv(8) )
   ALLOCATE( fractal(1024) )

!
! Reads from the external file 'status.inp'
! the status of a previous run (if any)
!     stat: last output of a previous run
!     mult: time step multiplier

   IF (myrank.eq.0) THEN
      OPEN(1,file='status.inp',status='unknown')
      READ(1,*) stat
      READ(1,*) time
      CLOSE(1)
      OPEN(1,file='old_status.inp')
      WRITE(1,*) stat,'  % stat'
      WRITE(1,*) time,'  % time'
      CLOSE(1)
   ENDIF
   CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!      IF (myrank.eq.0) THEN
!         OPEN(1,file='xyfpos.inp',status='unknown')
!         READ(1,*) xfpos
!         READ(1,*) yfpos
!         CLOSE(1)
!      ENDIF
!      CALL MPI_BCAST(xfpos,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!      CALL MPI_BCAST(yfpos,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

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
!     seed : seed for random function
!     prm1 : free parameter 1
!     prm2 : free parameter 2
!     ldir : local directory for I/O

   IF (myrank.eq.0) THEN
      OPEN(1,file='parameter.inp',status='unknown')
      READ(1,*) CFL                      ! 1
      READ(1,*) step                     ! 2
      READ(1,*) tstep                    ! 3
      READ(1,*) sstep                    ! 4
      READ(1,*) cstep                    ! 5
      READ(1,*) f0                       ! 6
      READ(1,*) u0                       ! 7
      READ(1,*) kdn                      ! 8
      READ(1,*) kup                      ! 9
      READ(1,*) nu                       ! 10
      READ(1,*) inu                      ! 11
      READ(1,*) mu                       ! 12
      READ(1,*) imu                      ! 13
      READ(1,*) iflow                    ! 14
      READ(1,*) seed                     ! 15
      READ(1,*) prm1                     ! 16
      READ(1,*) prm2                     ! 17
      READ(1,'(a100)') ldir              ! 18
      CLOSE(1)
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
   CALL MPI_BCAST( seed,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) !14
   CALL MPI_BCAST(iflow,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) !15
   CALL MPI_BCAST( prm1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !16
   CALL MPI_BCAST( prm2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !17
   CALL MPI_BCAST( ldir,100,MPI_CHARACTER,     0,MPI_COMM_WORLD,ierr) !18
   myseed = seed+myrank
!
! Some numerical constants

   ic = 48
   id = 48
   iu = 48
   jt = 48
   jc = 48
   jd = 48
   ju = 48

!
! Some constants for the FFT
!     kmax: maximum truncation for dealiMR=           0 ii,jj         944          13 vx,vy -0.664422238994306asing
!     tiny: minimum truncation for dealiasing

   kmax = (dble(n)/3.d0)
   kmax2= kmax**2
   tiny =  0.000001d0
   kmin =  tiny
   kmin2=  kmin**2
   enerk = 1.0d0
   enst =  1.0d0

! Builds the wave number and the square wave
! number matrixes

   DO i = 1,n/2
      ka(i) = dble(i-1)
      ka(i+n/2) = dble(i-n/2-1)
   END DO
   DO i = ista,iend
      DO j = 1,n
         ka2(j,i) = ka(i)**2+ka(j)**2
      END DO
   END DO

! Initializes the FFT library
! Use FFTW_ESTIMATE in short runs and FFTW_MEASURE
! in long runs

   CALL fftp2d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
      FFTW_MEASURE)
   CALL fftp2d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
      FFTW_MEASURE)

!tmp1=randu(seed)/sqrt(dt)
! Sets the initial conditions.
!!!!!! INITIAL CONDITIONS !!!!!!!!!!!!!!!
   IF (stat.eq.0) THEN
      ini = 1
      timet = 0 !!tstep
      times = sstep
      timec = cstep

!        STREAM FUNCTION R1
      DO i = ista,iend
         DO j = 1,n
            !! IF ((ka2(j,i).le.kup*kup).and.(ka2(j,i).ge.kdn*kdn)) THEN
            IF ((ka2(j,i).le.kup*kup).and.(ka2(j,i).ge.kmin2  )) THEN
               tmp = ka2(j,i)/(kdn+1)**2
               if (tmp.gt.80.0d0) tmp=80.0d0
               phase1 = 2*pi*randu(myseed)
               ps(j,i) = exp(im*phase1)*exp(-tmp)
            ELSE
               ps(j,i) = 0.
            ENDIF
         END DO
      END DO
      IF (myrank.eq.0) THEN
         DO j = 1,n
            jj=n-j+2
            if (jj.gt.n) jj=1
            ps(j,1) = conjg(ps(jj,1))
         END DO
      ENDIF
      CALL energy(ps,enerk,1)
      CALL MPI_BCAST(enerk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      tmp=u0/sqrt(enerk)
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
               ps(j,i) = tmp*ps(j,i)
            ELSE
               ps(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
      !! if (myrank.eq.0) print*,"PRCV 1"
      !! CALL parseval(ps,ps)

   ELSE
      print*,'READING...',stat
      ini = int((stat-1)*tstep)
      dump = dble(ini)/dble(sstep)+1
      times = 0
      timet = 0
      timec = 0

      jt = 48+int(dump/1000)
      jc = 48+int(dump/100)-int(dump/1000)*10
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
      if (myrank.eq.0) print*,"READING DONE!"
      CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)
      if (myrank.eq.0) print*,"READING & FFT DONE!"
      CALL parseval(ps,ps)
   ENDIF

!!!!!!! STEADY oo RANDOM FORCING  !!!!!!!!!!!
   IF (iflow.le.3) THEN
      CALL CFL_condition(CFL,ps,inu,nu,dt)
      CALL forcing(iflow,f0,kup,kdn,dt,myseed,fk)
      !print*,myrank,iflow,f0,kup,kdn,dt,myseed
   ENDIF
   ext4='ffff'
   CALL spectrum(fk,ext4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   DO i = ista,iend
      DO j = 1,n
         C1(j,i) = ps(j,i)/dble(n)**2
      END DO
   END DO
   CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
   OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
      // '000' // '.out',form='unformatted')
   WRITE(1) R1
   CLOSE(1)
   DO i = ista,iend
      DO j = 1,n
         C1(j,i) = ps(j,i)*ka2(j,i)/dble(n)**2
      END DO
   END DO
   CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
   OPEN(1,file=trim(ldir) // '/hd2Dww.' // node // '.' &
      // '000' // '.out',form='unformatted')
   WRITE(1) R1
   CLOSE(1)
   CALL forcing(iflow,f0,kup,kdn,dt,myseed,fk)  !! set fk=0
   CALL energy(fk,enerk,1)
   tmp1=f0/sqrt(0.5*enerk)
   CALL MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   fk  = tmp1*fk
   DO i = ista,iend
      DO j = 1,n
         C1(j,i) = fk(j,i)*ka2(j,i)/dble(n)**2
      END DO
   END DO
   CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
   OPEN(1,file=trim(ldir) // '/hd2Dfw.' // node // '.' &
      // '000' // '.out',form='unformatted')
   WRITE(1) R1
   CLOSE(1)
   DO i = ista,iend
      DO j = 1,n
         C1(j,i) = fk(j,i)/dble(n)**2
      END DO
   END DO
   CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
   OPEN(1,file=trim(ldir) // '/hd2Dfp.' // node // '.' &
      // '000' // '.out',form='unformatted')
   WRITE(1) R1
   CLOSE(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (myrank.eq.0) print*,"START RK",step
! Time integration scheme starts here
! Uses Runge-Kutta of order 'ord'
!#################### MAIN LOOP ######################
   ini=1
   RK : DO t = ini,step
      !if (myrank.eq.0) print*,"DBG RK",t
      CALL CFL_condition(CFL,ps,inu,nu,dt)
      !if (myrank.eq.0) print*,"DBG RK",t,dt



!!!!!!!  RANDOM FORCING  !!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL forcing(iflow,f0,kup,kdn,dt,myseed,fk)  !! set fk=0
      CALL energy(fk,enerk,1)
      tmp1=f0/sqrt(0.5*enerk*dt)*randu(seed)
      CALL MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      fk  = tmp1*fk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Every 'cstep' steps, generates external files
! to check consistency and convergency. See the
! hdcheck subroutine for details.
      IF (timec.eq.cstep) THEN
         timec = 0
         CALL hdcheck(ps,fk,time,inu,nu,imu,mu,ener,enst)
         if (myrank.eq.0) print*,time,ener,enst
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
         IF (jc.eq.58) THEN
            jc = 48
            jt = jt+1
         ENDIF
         th= char(jt)
         c = char(jc)
         d = char(jd)
         u = char(ju)
         ext4 = th // c // d // u
         CALL spectrum(ps,ext4)
         CALL laplak2(ps,C1)     ! make W
         CALL vectrans(ps,ps,C1,'euu',ext4)
         CALL vectrans(C1,ps,C1,'vrt',ext4)
         !CALL structure(ps,ext4)
         IF (myrank.eq.0) THEN
            OPEN(1,file='spectra_times.txt',position='append')
            WRITE(1,13) ext4,time
13          FORMAT( A4,    F12.6)
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
         DO i = ista,iend
            DO j = 1,n
               C1(j,i) = ps(j,i)/dble(n)**2
            END DO
         END DO
         CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
         OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
            // c // d // u // '.out',form='unformatted')
         WRITE(1) R1
         CLOSE(1)
         DO i = ista,iend
            DO j = 1,n
               C1(j,i) = ps(j,i)*ka2(j,i)/dble(n)**2
            END DO
         END DO
         CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
         OPEN(1,file=trim(ldir) // '/hd2Dww.' // node // '.' &
            // c // d // u // '.out',form='unformatted')
         WRITE(1) R1
         CLOSE(1)
!            DO i = ista,iend
!               DO j = 1,n
!                  C1(j,i) = fk(j,i)*ka2(j,i)/dble(n)**2
!               END DO
!            END DO
!            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
!            OPEN(1,file=trim(ldir) // '/hd2Dfw.' // node // '.' &
!                 // c // d // u // '.out',form='unformatted')
!            WRITE(1) R1
!            CLOSE(1)
!            DO i = ista,iend
!               DO j = 1,n
!                  C1(j,i) = fk(j,i)/dble(n)**2
!               END DO
!            END DO
!            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
!            OPEN(1,file=trim(ldir) // '/hd2Dfp.' // node // '.' &
!                 // c // d // u // '.out',form='unformatted')
!            WRITE(1) R1
!            CLOSE(1)

         IF (myrank.eq.0) THEN
            OPEN(1,file='field_times.txt',position='append')
            WRITE(1,12) c//d//u,time
12          FORMAT( A3,    F12.6)
            CLOSE(1)
            !OPEN(1,file='status.inp')
            !WRITE(1,*) c//d//u,'          % stat'
            !WRITE(1,*) time,'  % time'
            !CLOSE(1)
         ENDIF
      ENDIF

      timet = timet+1
      times = times+1
      timec = timec+1
      time = time+dt

! Runge-Kutta step 1
! Copies the streamfunction into the auxiliary matrix C1

      DO i = ista,iend
         DO j = 1,n
            C1(j,i) = ps(j,i)
         END DO
      END DO

! Runge-Kutta step 2

      DO o = ord,1,-1

         CALL laplak2(C1,C2)     ! make W
         CALL poisson(C1,C2,C1)  ! u grad w

         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
                  C1(j,i) = ps(j,i)+ dt*(-C1(j,i)/ka2(j,i) + fk(j,i))/dble(o)
                  tmp     = exp(-(mu/ka2(j,i)**imu + nu*ka2(j,i)**inu)*dt/dble(o))
                  !! tmp     = 1.0D0 /(1.0d0 + (mu/ka2(j,i)**imu+nu*ka2(j,i)**inu)*dt/dble(o))
                  C1(j,i) = C1(j,i) *tmp
               ELSE
                  C1(j,i) = 0.0d0
               ENDIF
            END DO
         END DO

      END DO

! Runge-Kutta step 3
! Copies the result from the auxiliary matrixes into ps

      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
               ps(j,i) = C1(j,i)
            ELSE
               ps(j,i) = 0.0d0
            ENDIF
         END DO
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   END DO RK
!##############  END OF MAIN LOOP ###################

!
! End of Runge-Kutta

   CALL MPI_FINALIZE(ierr)
   CALL fftp2d_destroy_plan(plancr)
   CALL fftp2d_destroy_plan(planrc)
   DEALLOCATE( R1 )
   DEALLOCATE( ps,fk )
   DEALLOCATE( C1,C2 )
   DEALLOCATE( ka,ka2 )

END PROGRAM HD2D
