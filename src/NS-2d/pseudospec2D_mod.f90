!=================================================================
! MODULES for 2D codes
!
! 2003 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

!=================================================================

MODULE grid
!
! n: number of points in the spatial grid
! choose: 128   256  512  1024  2048  4096  8192  16384
   INTEGER :: n = 256
   SAVE

END MODULE grid
!=================================================================

MODULE fft
!
!      USE fftplans
!      TYPE (FFTPLAN) :: planrc, plancr
!      SAVE
   INCLUDE 'fftw_f77.i'
   INTEGER, PARAMETER  :: ikind = 8
   INTEGER(kind=ikind) :: planrc,plancr
   SAVE


END MODULE fft
!=================================================================

MODULE ali
   DOUBLE PRECISION :: kmax
   DOUBLE PRECISION :: tiny
   SAVE

END MODULE ali
!=================================================================

MODULE var
   DOUBLE PRECISION    :: pi = 3.1415926535897932384d0
   DOUBLE COMPLEX :: im = (0.,1.)
   SAVE

END MODULE var
!=================================================================

MODULE kes
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: ka
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: ka2
   SAVE

END MODULE kes
!=================================================================

MODULE random
CONTAINS
   REAL FUNCTION randu(idum)
!
! Uniform distributed random numbers between -1 and
! 1. The seed idum must be between 0 and the value
! of mask

      INTEGER, PARAMETER :: iq=127773,ir=2836,mask=123459876
      INTEGER, PARAMETER :: ia=16807,im=2147483647
      INTEGER            :: k,idum
      REAL, PARAMETER    :: am=1./im

      idum = ieor(idum,mask)
      k = idum/iq
      idum = ia*(idum-k*iq)-ir*k
      IF (idum.lt.0) idum = idum+im
      randu = am*idum
      randu = (randu-.5)*2
      idum = ieor(idum,mask)

   END FUNCTION randu

END MODULE random
!=================================================================
MODULE fftplans
!
! Set the variable ikind to:  4 in 32 bits machines
!                             8 in 64 bits machines
!                             7 when using GNU compilers
! Set the variable csize to:  8 if L1 cache is <= 64 kb
!                            16 if L1 cache is 128 kb
   INCLUDE 'fftw_f77.i'
   INTEGER, PARAMETER  :: ikind = 8
   INTEGER, PARAMETER  :: csize = 16
   TYPE FFTPLAN
      INTEGER(kind=ikind) :: planr,planc
      INTEGER :: n
      INTEGER, DIMENSION (:), POINTER :: itype1, itype2
   END TYPE FFTPLAN
   SAVE

END MODULE fftplans
!=================================================================

!  MODULE mpivars
!      INCLUDE 'mpif.h'
!      INTEGER, SAVE :: ista,iend
!      INTEGER, SAVE :: jsta,jend
!      INTEGER, SAVE :: ksta,kend
!      INTEGER, SAVE :: nprocs,myrank
!      INTEGER, SAVE :: ierr
!
!  END MODULE mpivars
!=================================================================
