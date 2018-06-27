!======================================================================!
!     Title  : nscat.f90                                               !
!     Author : Yusa Shusaku                                            !
!     Date   : 2008-6-20-Fri ~ 2008-6-21-Sat                           !
!     Last modified : 2008-6-22-Sun                                    !
!                                                                      !
!     A program which computes angular distribution of differential    !
!     cross section for elastic neutron-nucleus scattering.            !
!     From the theory of scattering, the differential cross section    !
!     is equal to the square of scattering amplitude. Since the        !
!     scattering amplitude is expressed in terms of phase shifts,      !
!     we can obtain differential cross section by computing them.      !
!                                                                      !
!     *** Structure of this program ***                                !
!                                                                      !
!     module      Com_var                                              !
!     program     main                                                 !
!     subroutine  Inn_Sol                                              !
!     subroutine  Numerov                                              !
!     subroutine  Phase_Shift                                          !
!     function    dsigma                                               !
!     subroutine  Tot_Cros_Sec                                         !
!     function    kk                                                   !
!     function    V                                                    !
!     function    jL                                                   !
!     function    nL                                                   !
!     function    PL                                                   !
!     subroutine  Title                                                !
!     subroutine  Display                                              !
!     subroutine  Gnuplot                                              !
!                                                                      !
!======================================================================!
      module Com_var
!----------------------------------------------------------------------!
!     Definitions of various global constants :                        !
!     Lc    ... Cutoff angular momentum.                               !
!     rmax  ... Range of nuclear potential.                            !
!     E     ... Collision energy in CM system.                         !
!     Amass ... Mass number of target nucleus.                         !
!     rmass ... Reduced mass.                                          !
!     k     ... Wave number.                                           !
!     Rn    ... Radius of target nucleus.                              !
!----------------------------------------------------------------------!
      implicit none
      integer, parameter :: Lc=14
      real*8,  parameter, private :: r0=1.2d0
      real*8,  parameter :: PI=3.1415926535897932d0, dx=0.0007d0
      real*8,  parameter :: rmax=15.0d0, Amass=40.0d0
      real*8,  parameter :: Mn=938.0d0, hbarc=197.329d0
      real*8,  parameter :: E=10.0d0
      real*8,  parameter :: V0=50.0d0, ALPHA=0.63d0
      real*8,  parameter :: rmass=Mn * Amass / (1.0d0 + Amass) 
      real*8,  parameter :: k=sqrt(2.0d0 * E * rmass) / hbarc
      real*8,  parameter :: Rn=r0 * Amass ** (1.0d0/3.0d0)
      end module
!======================================================================!
      program main
!----------------------------------------------------------------------!
!     Main program.                                                    !
!----------------------------------------------------------------------!
      use Com_var, only : Lc
      implicit none
      integer :: L
      real*8  :: delta(0:Lc)
      
      call Title
      do L=0, Lc
         call Phase_Shift(delta)
      end do
      call Display(delta)
      !call Gnuplot(delta, 7, 'gnu')

      stop
      end program
!======================================================================!
      subroutine  Inn_Sol(L, a, b, s)
!----------------------------------------------------------------------!
!     A subroutine which integrates schroedinger equation from         !
!     origin to rmax using Numerov method.                             !
!     After that we compute 's' which is necessary for calculating     !
!     phase shifts.                                                    !
!----------------------------------------------------------------------!
      use Com_var, only : rmax, dx
      implicit none
      integer, intent(in) :: L
      real*8, intent(out) :: a, b, s
      real*8 :: y(3), x, ya

      x = 1.0d-10
      y(3) = x ** (L + 1)
      y(2) = y(3) + dble(L + 1) * x ** L * dx

      x = x + 2.0d0 * dx
      do
         call Numerov(L, x, y)
         if (x > rmax) exit
         x = x + dx
         y(3) = y(2)
         y(2) = y(1)
      end do

      a = x - 1.0d0 * dx
      ya = y(2)
      b = x
      s = b * ya / (a * y(1))

      return
      end subroutine
!======================================================================!
      subroutine  Numerov(L, x, y)
!----------------------------------------------------------------------!
!     A subroutine which computes 'y(x+dx)' from the values 'y(x)'     !
!     and 'y(x-dx).                                                    !
!----------------------------------------------------------------------!
      use Com_var, only : dx
      implicit none
      integer, intent(in) :: L
      real*8, intent(in) :: x
      real*8, intent(inout) :: y(3)
      real*8, external :: kk
      real*8, parameter :: wh=dx*dx/12.0d0
      real*8  :: ynm2, ynm1

      ynm1 = (2.0d0 - 10.0d0 * wh * kk(L,x-dx)) * y(2)
      ynm2 = (1.0d0 + wh * kk(L,x-2.0d0*dx)) * y(3)
      y(1) = (ynm1 - ynm2) / (1.0d0 + wh * kk(L,x))

      return
      end subroutine
!======================================================================!
      subroutine  Phase_Shift(delta)
!----------------------------------------------------------------------!
!     A subroutine which computes phase shift with angular momentum    !
!     'L'.                                                             !
!----------------------------------------------------------------------!
      use Com_var, only : Lc, k
      implicit none
      integer :: L
      real*8, intent(out) :: delta(0:Lc) 
      real*8, external :: jL, nL
      real*8  :: s, a, b, arg

      do L=0, Lc
         call Inn_Sol(L, a, b, s)
         arg = (s * jL(L,k*b) - jL(L,k*a)) / (s * nL(L,k*b) - nL(L,k*a))
         delta(L) = atan(arg)
      end do

      return
      end subroutine
!======================================================================!
      function  dsigma(delta, t)  result(ds)
!----------------------------------------------------------------------!
!     A function which returns differential cross section at angle     !
!     't'. At first, we compute scattering amplitude from the phase    !
!     shifts and compute differential cross section by squaring it.    !
!----------------------------------------------------------------------!
      use Com_var, only : Lc, k
      implicit none
      integer :: L
      real*8, intent(in) :: delta(0:Lc), t
      real*8, external :: PL
      real*8 :: ds
      complex*16, parameter :: i=(0.0d0, 1.0d0)
      complex*16 :: f

      f = (0.0d0, 0.0d0)
      do L=Lc, 0, -1
         f = f + dble(2 * L + 1) * exp(i*delta(L)) * sin(delta(L)) &
     &         * PL(L,cos(t))
      end do
      f = f / k 
      ds = abs(f) ** 2

      return
      end function
!======================================================================!
      subroutine  Tot_Cros_Sec(delta, tsig)
!----------------------------------------------------------------------!
!     A subroutine which computes total cross section from the         !
!     phase shifts.                                                    !
!----------------------------------------------------------------------!
      use Com_var, only : Lc, PI, k
      implicit none
      integer :: L
      real*8, intent(in)  :: delta(0:Lc)
      real*8, intent(out) :: tsig

      tsig = 0.0d0
      do L=Lc, 0, -1
         tsig = tsig + dble(2 * L + 1) * sin(delta(L)) ** 2
      end do
      tsig = tsig * 4.0d0 * PI / (k * k)

      return
      end subroutine
!======================================================================!
      function kk(L, x)  result(k)
!----------------------------------------------------------------------!
!     Definition of function 'kk' which appears in the Schroedinger    !
!     equation when we write it in the following form :                !
!      d^2y                                                            !
!     ------ + kk * y = 0 .                                            !
!      dx^2                                                            !
!----------------------------------------------------------------------!
      use Com_var, only : rmass, hbarc, E
      implicit none
      integer, intent(in) :: L
      real*8, intent(in) :: x
      real*8, external :: V
      real*8 :: k
      
      k = 2.0d0 * rmass * (E - V(x)) / (hbarc * hbarc)  &
     &        - dble((L * (L + 1))) / (x * x)

      end function
!======================================================================!
      function  V(x)  result(WS) 
!----------------------------------------------------------------------!
!     Definition of Woods-Saxon potential.                             !
!----------------------------------------------------------------------!
      use Com_var, only : V0, ALPHA, Rn
      implicit none
      real*8, intent(in) :: x
      real*8 :: WS

      WS = - V0 / (1.0d0 + exp((x - Rn) / ALPHA))

      return
      end function
!======================================================================!
      function  jL(L, x)  result(sphB)
!----------------------------------------------------------------------!
!     Definition of sphrical Bessel function.                          !
!     We use recursion relation in computing it.                       !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: L
      integer :: i
      real*8, intent(in) :: x
      real*8 :: sphB, j0, j1

      j0 = sin(x) / x
      j1 = (sin(x) - x * cos(x)) / (x * x)

      if (L >= 2) then
         do i=2, L
            sphB = dble(2 * i - 1) / x * j1 - j0
            j0 = j1
            j1 = sphB
         end do
      else if (L == 0) then
         sphB = j0
      else if (L == 1) then
         sphB = j1
      end if

      return
      end function
!======================================================================!
      function  nL(L, x)  result(sphN)
!----------------------------------------------------------------------!
!     Definition of sphrical Neumann function.                         !
!     We use recursion relation in computing it.                       !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: L
      integer :: i
      real*8, intent(in) :: x
      real*8 :: sphN, n0, n1

      n0 = - cos(x) / x
      n1 = - (cos(x) + x * sin(x)) / (x * x)

      if (L >= 2) then
         do i=2, L
            sphN = dble(2 * i - 1) / x * n1 - n0
            n0 = n1
            n1 = sphN
         end do
      else if (L == 0) then
         sphN = n0
      else if (L == 1) then
         sphN = n1
      end if

      return
      end function
!======================================================================!
      function  PL(L, x)  result(P)
!----------------------------------------------------------------------!
!     Definition of Legendre polynomials.                              !
!     We use recursion relation in computing it.                       !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: L
      integer :: i
      real*8, intent(in) :: x
      real*8 :: P0, P1, P

      P0 = 1.0d0
      P1 = x

      if (L >= 2) then
         do i=2, L
            P = (dble(2*i - 1) * x * P1 - dble(i - 1) * P0) / dble(i)
            P0 = P1
            P1 = P
         end do
      else if (L == 0) then
         P = P0
      else if (L == 1) then
         P = P1
      end if

      return
      end function
!======================================================================!
      subroutine  Title

      write(6,*)
      write(6,*) '    ********************************'
      write(6,*) '       Neutron-Nucleus Scattering   '
      write(6,*) '     ~Computation of cross section~ '
      write(6,*) '    ********************************'
      write(6,*)

      return
      end subroutine
!======================================================================!
      subroutine  Display(delta)
!----------------------------------------------------------------------!
!     A subroutine which computes angular distribution of the          !
!     differential cross section and the total cross section and       !
!     displays them on the screen.                                     !
!----------------------------------------------------------------------!
      use Com_var, only : PI, Lc
      implicit none
      integer :: i
      real*8 :: t, tsig
      real*8, intent(in) :: delta(0:Lc)
      real*8, external :: dsigma
      character(len=25), parameter :: FM1='(2x,i4,3x,a,5x,f12.7,a)'
      character(len=15), parameter :: FM2='(2x,a,f8.3,a)'

      write(6,*) '  Angle | Diff. Cross Section'
      write(6,*) '---------------------------------'
      do i=0, 180, 10
         t = dble(i) * PI / 180.0d0 
         write(6,FM1) i,'|',dsigma(delta, t), '   fm^2'
      end do
      write(6,*)
      write(6,*) 'Absolute error < 10^(-7)'
      write(6,*)
      write(6,*) '*** Total Cross Section ***'
      call Tot_Cros_Sec(delta, tsig)
      write(6,FM2) 'sigma = ', tsig ,'  fm^2'
      write(6,*)
      write(6,*) 'Type following commands :'
      write(6,*) '%> gnuplot < gnu'
      write(6,*) '%> gv Diff_CS.eps'
      write(6,*)

      return
      end subroutine
!======================================================================!
      subroutine  Gnuplot(delta, Fnum, Fname)
!----------------------------------------------------------------------!
!     A subroutine which makes a file for gnuplot. If we redirect      !
!     it into the gnuplot, we can obtain an eps file. This file is     !
!     a graph of differential cross section as a function of angle.    !
!----------------------------------------------------------------------!
      use Com_var, only : Lc, PI, E
      implicit none
      integer, intent(in) :: Fnum
      integer :: i
      real*8, intent(in) :: delta(0:Lc)
      real*8, external :: dsigma
      real*8 :: t
      character(len=*), intent(in) :: Fname
      character(len=17), parameter :: FM='(1x,f8.3,f13.8)'
      character(len=15), parameter :: Q='(1x,a,f5.1,a)'
      
      open(Fnum,file=Fname)
      write(Fnum,*) 'set term postscript eps enhanced color'
      write(Fnum,*) "set output 'Diff_CS.eps"
      write(Fnum,*) "set title 'Differential Cross Section'"
      write(Fnum,*) "set xlabel '{/Symbol q}_{cm} (degree)'"
      write(Fnum,*) "set ylabel 'd{/Symbol s} / d{/Symbol W} (mb)'"
      write(Fnum,*) "set label 'n + ^{40}Ca' at 80,0.4"
      write(Fnum,Q) "set label 'E_{cm} =", E, " MeV'at 80,0.15"
      write(Fnum,*) 'unset key'
      write(Fnum,*) 'set logscale y'
      write(Fnum,*) 'set format y "10^{%L}"'
      write(Fnum,*) 'set size 0.5, 0.5'
      write(Fnum,*) "plot '-' with line"
      do i=0, 900
         t = dble(i) * PI / 900.0d0 
         write(Fnum,FM)  0.2d0 * dble(i), dsigma(delta,t)*0.001d0
      end do
      close(Fnum)
      
      return
      end subroutine
