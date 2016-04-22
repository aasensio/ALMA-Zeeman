module dustOpacity
use iso_c_binding, only : c_double, c_int
implicit none
contains
       subroutine kappa_dust(freq,iddust,adust,kdust) bind(c)

!      Input: freq                       [Hz]
!             type of dust 'iddust'
!             dust grains radius 'adust' [um]
!      Returns dust absortion coefficient 'kdust' [cm2/g(dust)]
!              and 'dustok' if called from command interpreter 

       implicit none

	   real(c_double), intent(in) :: freq
	   real(c_double), intent(in) :: adust
	   integer(c_int), intent(in) :: iddust
	   real(c_double), intent(out) :: kdust   ! dust absortion coefficient [cm2/g(dust)]
	   
       real*8 hplanck,kboltz,clight,pi
       parameter(hplanck=6.62606896d-27) ! Planck constant    [erg s] 
       parameter(kboltz=1.3806504d-16)   ! Boltzmann constant [erg/K] 
       parameter(clight=2.99792458d10)   ! speed of light     [cm/s] 
       parameter(pi=3.1415926536d0)      ! pi number 

       integer nmaxiddust,nmaxwl
       parameter(nmaxiddust=1)           ! Number of dust types = 1
       parameter(nmaxwl=100)

       integer i,j,jl,ndata(nmaxiddust)
       real*8 x0,n0,k0,wl0,qabs,slope
       real*8 wl(nmaxiddust,nmaxwl)      ! wavelength [um]
       real*8 n(nmaxiddust,nmaxwl)       ! n index
       real*8 k(nmaxiddust,nmaxwl)       ! k index
       real*8 rhodust(nmaxiddust)        ! dust mass density  [g cm-3]       
       complex*16 ri0

       logical dustok
       logical wlincr                    ! data ordered by increasing(decreasing) wavelentgh
       
       dustok = .TRUE.

!      iddust=1: AC (amorphous carbon); Suh 2000, MNRAS, 315, 740
!      complex index of refraction (n,k) as a function of wavelength [um]
       ndata(1)=98
       rhodust(1)=2.0d0       ! g cm-3
       data (wl(1,i),i=1,98)/& ! wavelength [um]
      .3600E+05,.2000E+05,.1300E+05,.4000E+04,.1300E+04,.7000E+03,&
      .5000E+03,.4000E+03,.3000E+03,.2500E+03,.2000E+03,.1500E+03,&
      .1400E+03,.1300E+03,.1200E+03,.1100E+03,.1050E+03,.1000E+03,&
      .9500E+02,.9000E+02,.8500E+02,.8000E+02,.7500E+02,.7000E+02,&
      .6500E+02,.6000E+02,.5500E+02,.5000E+02,.4500E+02,.4000E+02,&
      .3500E+02,.3000E+02,.2800E+02,.2700E+02,.2600E+02,.2500E+02,&
      .2400E+02,.2300E+02,.2200E+02,.2000E+02,.1900E+02,.1800E+02,&
      .1700E+02,.1600E+02,.1500E+02,.1450E+02,.1400E+02,.1350E+02,&
      .1300E+02,.1250E+02,.1200E+02,.1160E+02,.1130E+02,.1100E+02,&
      .1050E+02,.1000E+02,.9850E+01,.9700E+01,.9550E+01,.9400E+01,&
      .9000E+01,.8500E+01,.8000E+01,.7500E+01,.7000E+01,.6500E+01,&
      .6000E+01,.5500E+01,.5000E+01,.4500E+01,.4000E+01,.3500E+01,&
      .3000E+01,.2700E+01,.2200E+01,.2000E+01,.1700E+01,.1300E+01,&
      .1150E+01,.1000E+01,.8500E+00,.7000E+00,.5500E+00,.4400E+00,&
      .3600E+00,.3000E+00,.2500E+00,.2000E+00,.1500E+00,.1200E+00,&
      .1000E+00,.8000E-01,.6000E-01,.5000E-01,.4000E-01,.3000E-01,&
      .2000E-01,.1000E-01/
       data (n(1,i),i=1,98)/& ! complex index of refraction n
      .2805E+01,.2761E+01,.2774E+01,.2729E+01,.2697E+01,.2656E+01,&
      .2635E+01,.2619E+01,.2597E+01,.2581E+01,.2561E+01,.2534E+01,&
      .2525E+01,.2517E+01,.2508E+01,.2498E+01,.2493E+01,.2487E+01,&
      .2481E+01,.2474E+01,.2467E+01,.2458E+01,.2449E+01,.2439E+01,&
      .2428E+01,.2415E+01,.2400E+01,.2384E+01,.2364E+01,.2340E+01,&
      .2310E+01,.2257E+01,.2232E+01,.2219E+01,.2207E+01,.2194E+01,&
      .2183E+01,.2172E+01,.2159E+01,.2138E+01,.2126E+01,.2114E+01,&
      .2101E+01,.2088E+01,.2075E+01,.2067E+01,.2060E+01,.2050E+01,&
      .2042E+01,.2032E+01,.2024E+01,.2015E+01,.2010E+01,.2002E+01,&
      .1991E+01,.1981E+01,.1978E+01,.1972E+01,.1969E+01,.1963E+01,&
      .1952E+01,.1935E+01,.1921E+01,.1906E+01,.1891E+01,.1875E+01,&
      .1858E+01,.1840E+01,.1821E+01,.1801E+01,.1780E+01,.1757E+01,&
      .1733E+01,.1715E+01,.1685E+01,.1671E+01,.1650E+01,.1622E+01,&
      .1612E+01,.1605E+01,.1605E+01,.1624E+01,.1574E+01,.1535E+01,&
      .1481E+01,.1428E+01,.1371E+01,.1299E+01,.1240E+01,.1184E+01,&
      .1016E+01,.9526E+00,.9146E+00,.9249E+00,.9307E+00,.9526E+00,&
      .9826E+00,.9896E+00/
       data (k(1,i),i=1,98)/& ! complex index of refraction k
      .2670E-01,.3301E-01,.3796E-01,.5727E-01,.8446E-01,.1055E+00,&
      .1191E+00,.1291E+00,.1433E+00,.1533E+00,.1665E+00,.1853E+00,&
      .1903E+00,.1957E+00,.2017E+00,.2086E+00,.2122E+00,.2164E+00,&
      .2205E+00,.2251E+00,.2303E+00,.2357E+00,.2417E+00,.2485E+00,&
      .2560E+00,.2642E+00,.2737E+00,.2846E+00,.2972E+00,.3124E+00,&
      .3310E+00,.3566E+00,.3542E+00,.3526E+00,.3511E+00,.3496E+00,&
      .3474E+00,.3454E+00,.3432E+00,.3379E+00,.3375E+00,.3371E+00,&
      .3365E+00,.3359E+00,.3352E+00,.3350E+00,.3345E+00,.3346E+00,&
      .3342E+00,.3341E+00,.3338E+00,.3334E+00,.3333E+00,.3334E+00,&
      .3332E+00,.3327E+00,.3327E+00,.3328E+00,.3327E+00,.3332E+00,&
      .3330E+00,.3295E+00,.3254E+00,.3208E+00,.3160E+00,.3109E+00,&
      .3054E+00,.2994E+00,.2932E+00,.2863E+00,.2782E+00,.2696E+00,&
      .2596E+00,.2532E+00,.2406E+00,.2350E+00,.2253E+00,.2095E+00,&
      .2024E+00,.1940E+00,.1890E+00,.2222E+00,.2700E+00,.3064E+00,&
      .3338E+00,.3502E+00,.3601E+00,.3597E+00,.3298E+00,.4224E+00,&
      .3937E+00,.2782E+00,.1590E+00,.1127E+00,.7059E-01,.3488E-01,&
      .1096E-01,.1010E-02/

       if(.not.dustok)then
          if(iddust.gt.nmaxiddust)then
             write(*,*)' E- IDdust outside valid range:'
             write(*,*)'   =0: dust not considered'
             write(*,*)'   =1: AC (Suh 2000, MNRAS, 315, 740)'
             write(*,*)
             dustok=.false.
          else
             dustok=.true.
          endif       
          if(adust.lt.1.0d-3.or.adust.gt.100.0d0)then
             write(*,*)' E- aDUST outside valid range: 0.001-100 um'
             write(*,*)
             dustok=.false.
          else
             dustok=.true.
          endif       
          return
       endif

       wl0=(clight/freq)*1.0d4       ! wavelength [um]       
       x0=2.0d0*pi*adust/wl0         ! dust size parameter = 2pi*a_dust/wavelength

!      Interpolate/extrapolate complex index of refraction n,k in log units of wavelength 
       if(wl(iddust,2).gt.wl(iddust,1))wlincr=.true.
       if(wl(iddust,2).lt.wl(iddust,1))wlincr=.false.
       jl=0
       if(wlincr)then ! data ordered by increasing wavelength
          do j=1,ndata(iddust)-1
             if(log10(wl0).ge.log10(wl(iddust,j)).and.log10(wl0).le.log10(wl(iddust,j+1)))then
                jl=j             
             endif
          enddo
          if(log10(wl0).lt.log10(wl(iddust,1)))jl=1
          if(log10(wl0).gt.log10(wl(iddust,ndata(iddust))))  jl=ndata(iddust)-1
       else           ! data ordered by decreasing wavelength
          do j=1,ndata(iddust)-1
             if(log10(wl0).le.log10(wl(iddust,j)).and.log10(wl0).ge.log10(wl(iddust,j+1)))then
                jl=j             
             endif
          enddo
          if(log10(wl0).gt.log10(wl(iddust,1)))jl=1
          if(log10(wl0).lt.log10(wl(iddust,ndata(iddust)))) jl=ndata(iddust)-1
       endif
       slope=(n(iddust,jl+1)-n(iddust,jl))/(log10(wl(iddust,jl+1))-log10(wl(iddust,jl)))
       n0=n(iddust,jl)+slope*(log10(wl0)-log10(wl(iddust,jl)))
       slope=(k(iddust,jl+1)-k(iddust,jl))/(log10(wl(iddust,jl+1))-log10(wl(iddust,jl)))
       k0=k(iddust,jl)+slope*(log10(wl0)-log10(wl(iddust,jl)))

       ri0=dcmplx(n0,-k0)
       call mie(x0,ri0,qabs)

!      Compute kdust: k_nu = (3/4)Qabs/(a_dust*rho_dust)
       kdust=qabs*3.0d0/4.0d0/(adust*1.0d-4)/rhodust(iddust)

       return
       end subroutine kappa_dust



!_______________________________________________________________________

!*********** Light Scattering by Spherical Particles *************
!                                                                *
!   Calculations of extinction, scattering, absorption, etc      *
!  efficiency factors for homogeneous spheres (the Mie theory).  *
!................................................................*
!  Input data:          filename:   mie.dat                      *
!   ri = n-k*i: complex index of refraction                      *
!           Nx: number of sizes                                  *
!           x1: size 1                                           *
!           x2: size 2                                           *
!          ...: ...                                              *
!................................................................*
!  Output data:         filename:   mie.out                      *
!   Qext : extinction factors                                    *
!   Qsca : scattering factors                                    *
!   Qabs : absorption factors                                    *
!   Qbk  : backscattering factors                                *
!   Qpr  : radiation pressure factors                            *
!  albedo: particle's albedo                                     *
!   g    : asymmetry factor                                      *
!................................................................*
! NB! In order to treat very large particles,                    *
!     one needs to enlarge the parameter NTERMS.                 *
!................................................................*
! created by N.V. Voshchinnikov                                  *
! with a support of the Volkswagen Foundation (Germany)          *
! (c) 1989/98 Astronomical Institute, St.Petersburg University   *
!*****************************************************************
!
!      \_______________________   modified for LVG  _________________________/ 
!
       subroutine mie(x,ri,qabs)
!
!      \______________________________  end  ________________________________/ 
	  IMPLICIT REAL*8(A-H,O-Z)                                          
	  integer :: nterms, nx, i
      parameter(nterms=300)
      
      COMPLEX*16 ri
      DIMENSION xx(50)
!   20 FORMAT(2D14.10)                                                 
!   21 FORMAT(3X,'SPHERES: homogeneous')                                 
!   22 FORMAT(3X,'THEORY:  exact')                                        
!  223 FORMAT(I4)                                          
!   24 FORMAT(4X,'x',7X,'Qext',4X,'Qsca',4X,'Qabs',4X,'Qbk',5x,'Qpr',
!     *       4X,'Albedo',3x,'g')
!   25 FORMAT(13X,F8.4,9X,I3,6X,2D15.4,F10.2,5X,D12.3,F11.3)             
!   29 FORMAT(I4/(F8.4))                                                 
!   30 FORMAT(3x,'m = '2F14.10,'*i')
!   31 FORMAT(F8.4)                                                      
!   32 FORMAT(1H  )                                                      
!   33 FORMAT(1X,8(F8.4))
!  233 FORMAT(4X,i4,2F10.3,f11.1)                                             
!   34 FORMAT( 8X,'x = ',F7.3,5x,'Qext = ',f8.4)          
!   36 FORMAT(1X,64 ('='))
!   37 FORMAT(1X,60 ('*')//)                                             
!   38 FORMAT(1X,10I4)                                                   
!   41 FORMAT(1X,64('.'))
!   42 FORMAT(1X,64('-'))
!   44 FORMAT(1X,A8/(1X,2F8.4))                                          
!   48 FORMAT( 5X,F8.4,2X,2(F10.3,3X),3X,F8.2)                           
!
!*                  INPUT
!
! LVG        print *,'start Mie'
! LVG      open (unit=05,file='mie.dat',status='old',access='sequential')
! LVG      open (unit=07,file='mie.out',status='unknown',access='append')
! LVG*      PI = 4d0 * datan(1d0)
! LVG      READ (5,20) ri
! LVG      READ (5,29) nx,(xx(I),I=1,Nx)
!      \_______________________   modified for LVG  _________________________/ 
!
       nx=1
       xx(1)=x
!
!      \______________________________  end  ________________________________/ 
!
!-------------------------------------------------------------------
!***              * Spheres *
!-------------------------------------------------------------------
!*                  Efficiency factors
!
      qext=0d0
      qsca=0d0
      qabs=0d0
      qbk =0d0
      qpr =0d0
      alb =0d0
      g   =0d0
!
!*                  Exact theory
!
! LVG      write (7,22)
! LVG      write (7,30) ri
! LVG      write (7,41)
! LVG      write (7,24)
! LVG      write (7,42)
! LVG      write (*,42)

        do 276 i = 1, nx
        x=xx(i)
        if(x.le.1d-6) go to 444
        call shexq(ri,x,qext,qsca,qabs,qbk,qpr,alb,g)
  444   continue
! LVG        write (7,33) x,qext,qsca,qabs,qbk,qpr,alb,g
! LVG        write (*,*) x, qext
  276   continue
!
! LVG 1000 continue
! LVG      write (7,42)
! LVG      write (*,42)
!      \_______________________   modified for LVG  _________________________/ 
!
! LVG      STOP
       return
!
!      \______________________________  end  ________________________________/ 
      END subroutine mie
!--------------------------------------------------------------------
! **********   shexq - Sphers: homogeneous
!                      Theory: exact
!                      Results: efficiency factors
!--------------------------------------------------------------------
      SUBROUTINE shexq(RI,X,QEXT,QSCA,qabs,qbk,qpr,alb,g)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)                      
      integer :: nterms, num, num2, num1
      parameter(nterms=300)
      
      DIMENSION RU(2*nterms), BESJ(nterms), BESY(nterms), RA(nterms), RB(nterms)
      AX=1.0D0/X                                                        
      NUM=NM(X)                                                         
!      \_______________________   modified for LVG  _________________________/ 
!
!      NUM2=NM(cDSQRT(RI*DCONJG(RI))*X)
       ax2=cdsqrt(ri*dconjg(ri))*x
       num2=nm(ax2)
!
!      \______________________________  end  ________________________________/ 

        if(num.gt.nterms) then
        write(*,*) 'nterms, num=', nterms, num
        stop
        end if
        if(num2.gt.2*nterms) then
        write(*,*) '2*nterms, num2=', 2*nterms, num2
        stop
        end if


      CALL AA(AX,RI,NUM2,RU)                                            
      CALL BESSEL(AX,NUM,BESJ,BESY)                                     
      CALL AB(AX,RI,NUM,NUM1,RU,BESJ,BESY,RA,RB)
      CALL QQ1(AX,NUM1,QEXT,QSCA,qbk,qpr,RA,RB)                             
      qabs=qext-qsca
      alb=qsca/qext
      g=(qext-qpr)/qsca
      RETURN                                                            
      END SUBROUTINE shexq                                                  
!--------------------------------------------------------------------
! NM-auxiliary function for AA & BESSEL
!    (number NM is calculated using X)
! see: Trudy Astronom. Observ. LGU V.28,P.14,1971
!    for X>1 value of NM was raised
! August 1989, AO LGU
!--------------------------------------------------------------------
      FUNCTION NM(X)      
      real*8 x
      integer :: nm
      IF(X.LT.1) GO TO 11
      IF(X.GT.100) GO TO 12
      NM=INT(1.25*X+15.5)
      RETURN
   11 NM=INT(7.5*X+9.0)
      RETURN
   12 NM=INT(1.0625*X+28.5)
      RETURN
      END FUNCTION NM
!--------------------------------------------------------------------
! AA-subroutine for calculations of the ratio of the derivative
!    to the function for Bessel functions of half order with
!    the complex argument: J'(N)/J(N).
!    The calculations are given by the recursive expression
!    ``from top to bottom'' beginning from N=NUM.
!    RU-array of results.
!    A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
!    RI - complex refractive index.
! August 1989, AO LGU
!--------------------------------------------------------------------
      SUBROUTINE AA(A,RI,NUM,RU)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      integer :: i, j, num1, num, i1
      DIMENSION RU(NUM)
      S=A/RI
      RU(NUM)=(NUM+1.0D0)*S
      NUM1=NUM-1
      DO 13 J=1,NUM1
      I=NUM-J
      I1=I+1
      S1=I1*S
   13 RU(I)=S1-1.0D0/(RU(I1)+S1)
      RETURN
      END SUBROUTINE AA
!--------------------------------------------------------------------
! BESSEL-subroutine for calculations of the Bessel functions
!    of half order and first (J(N)) second (Y(N)) kinds with
!    the real argument X (A=1/X).
!    The calculations are started from N=NUM.
!    Desription of method see:
!    V.M.Loskutov, Trudy Astronom. Observ. LGU V.28,P.14,1971
!    BESJ-array of functions J(N), BESY-array of functions Y(N)
! August 1989, AO LGU
!--------------------------------------------------------------------
      SUBROUTINE BESSEL(A,NUM,BESJ,BESY)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z)
      integer :: num, i, n, i2, num1, i1
      DIMENSION BESJ(NUM+1),BESY(NUM+1) 
      BESJ(NUM+1)=0.0D0
      BESJ(NUM)=0.1D-11
      N=2*NUM+1
      NUM1=NUM-1
      DO 11 I=1,NUM1
      N=N-2
      I1=NUM-I
   11 BESJ(I1)=N*A*BESJ(I1+1)-BESJ(I1+2)
      B1=A*BESJ(1)-BESJ(2)
      B2=-A*B1-BESJ(1)
      N=2*(NUM/2)
      B=1.2533141373155002D0*DSQRT(A)
      C=B*BESJ(1)
      DO 12 I=3,N,2
      B=B*(I-0.5D0)*(I-2.0D0)/(I-2.5D0)/(I-1.0D0)
   12 C=C+B*BESJ(I)
      C=1.0D0/C
      DO 13 I=1,NUM
   13 BESJ(I)=C*BESJ(I)
      BESY(1)=-C*B1
      BESY(2)=C*B2
      DO 14 I=3,NUM
      I2=I-2
   14 BESY(I)=(2.0D0*I2+1.0D0)*A*BESY(I-1)-BESY(I2)
      RETURN
      END SUBROUTINE BESSEL
!--------------------------------------------------------------------
! AB-subroutine for calculations of the complex coefficients
!    A(N), B(N) for spherical particles.
!    A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
!    RI - complex refractive index.
!    The coefficients are calculated up to the number NUM1.LE.NUM,
!    for which |A(N)**2+B(N)**2|.LE.10**(-60)
!    RA-array of coefficients A(N), RB-array of coefficients B(N)
! August 1989, AO LGU
!--------------------------------------------------------------------
      SUBROUTINE AB(A,RI,NUM,NUM1,RU,BESJ,BESY,RA,RB)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      integer :: num, num1, i, n
      DIMENSION RU(NUM),BESJ(NUM),BESY(NUM),RA(NUM),RB(NUM)
      S3=(0.0D0,1.0D0)
      DO 11 I=1,NUM-1
      N=I+1
      S=RU(I)/RI+I*A
      S1=S*BESJ(N)-BESJ(I)
      S2=S*BESY(N)-BESY(I)
      RA(I)=S1/(S1-S3*S2)
      S=RU(I)*RI+I*A
      S1=S*BESJ(N)-BESJ(I)
      S2=S*BESY(N)-BESY(I)
      RB(I)=S1/(S1-S3*S2)
      P=RA(I)*DCONJG(RA(I))+RB(I)*DCONJG(RB(I))
      IF (P.LE.1D-60) GO TO 12
   11 CONTINUE
   12 NUM1=I
      RETURN
      END SUBROUTINE AB
!--------------------------------------------------------------------
! QQ1-subroutine for calculations of the efficiency factors for
!     extinction (QEXT), scattering (QSCA), backscattering (QBK)
!     and radiation pressure (QPR) for spherical particles.
! August 1989, AO LGU
!--------------------------------------------------------------------
      SUBROUTINE QQ1(A,NUM,QEXT,QSCA,qbk,qpr,RA,RB)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      DIMENSION RA(NUM),RB(NUM)
      integer :: num, i, n
      B=2.0D0*A*A
      C=0.0D0
      D=0.0D0
      s=(0d0,0d0)
      r=(0d0,0d0)
      N=1
      DO 11 I=1,NUM-1
      N=N+2      
      r=r+(i+0.5d0)*(-1)**i*(ra(i)-rb(i))      
      s=s+i*(i+2d0)/(i+1d0)*(ra(i)*dconjg(ra(i+1))+rb(i)*dconjg(rb(i+1)))+n/i/(i+1d0)*(ra(i)*dconjg(rb(i)))    
      C=C+N*(RA(I)+RB(I))
   11 D=D+N*(RA(I)*DCONJG(RA(I))+RB(I)*DCONJG(RB(I)))
      QEXT=B*C
      QSCA=B*D                          
      qbk=2d0*b*r*dconjg(r)
      qpr=qext-2d0*b*s
      RETURN
      END SUBROUTINE QQ1
!=== eof ===
end module dustOpacity