module params
! Defines model parameters.
implicit none
real, parameter :: one    = 1.0
real, parameter :: pi     = 3.14159
real, parameter :: dt     = 1.0                  ! Timestep                                         (days)
real, parameter :: k_Ca   = 0.00563              ! For future CO2                                      (-)
real, parameter :: a      = 0.1502               ! Allometry constant for mass from diameter (kg[DM] cm-1)
real, parameter :: b      = 2.476                ! Allometry constant for mass from diamter            (-)
real, parameter :: KS_def = 0.0288               ! Neutral CO2 concentraiton               (g[Su] g[DM]-1)
real, parameter :: n      = 4                    ! Hill coefficient                                    (-)
real, parameter :: alpha  = 2.0 * 0.0089 / 365.0 ! Source capacity            (kg[Su] cm[sapwood]-2 day-1)
real, parameter :: beta   = 2.0 * 0.329 / 365.0  ! Sink capacity                                (cm day-1)
integer, parameter :: nyr_co2 = 2101             ! Size of CO2 concentration array                     (n)
integer, parameter :: nt = nint (1.0 / dt)       ! No. timesteps in day                                (n)
integer, parameter :: sim = 5                    ! Simulation version                                  (n)
end module params

module vars
! Declare model variables.
use params
implicit none
real, dimension (nyr_co2) :: ca, ca_ssp   ! Mean annual atmospheric CO2 mixing ratio  (ppm)
real, dimension (nyr_co2*365*nt) :: ca_it ! Mean atmospheric CO2 mixing ratio at model timestep (ppm)
real :: ca0, ca1, ac, bc, D, M, S, dM, gG, fa, U, gA, fi, Sc, dS
real :: dS1, dS2, dS3, dS4, dM1, dM2, dM3, dM4, calib_G, calib_A
real :: KS
integer :: iyr, kyr, kday, it, i, j
end module vars

program sugar
use params
use vars
implicit none
write (*,*) alpha*365.0,beta*365.0
!stop
! Read annual mean CO2 mixing ratios, 1700-2023 CE (ppm)
open (10,file='global_co2_ann_1700_2023.txt',status='old')
do kyr = 1700, 2023
  read (10,*) iyr, ca (kyr)
end do
close (10)
! Read annual mean CO2 mixing ratios, 1850-2100 CE (ppm)
open (10,file='ssp585.txt',status='old')
do kyr = 1850, 2100
  read (10,*) iyr, ca_ssp (kyr)
end do
close (10)
! Use SSP585 CO2 for 2000-2100 CE.
ca (2000:2100) = ca_ssp (2000:2100)
ca (2101) = ca (2100) + (ca (2100) - ca (2099))
!do kyr = 2024, nyr_co2
!  ca (kyr) = ca (2023)* exp (k_Ca * (kyr-2023))
!end do
! For diagnostics of CO2.
open (20,file='ca.txt',status='unknown')
do kyr = 1700, nyr_co2
  write (20,*) kyr, ca (kyr)
end do
close (20)
! Assume Ca is on 1 jan.
j = 1699*365*nt + 1 !+ 182*nt
! Compute CO2 on model timestep.
do kyr = 1700, 2100
  ca0 = ca (kyr)
  ca1 = ca (kyr+1)
  bc = (ca1 - ca0) / (365*nt)
  ac = ca0 - bc * 1.0
  do kday = 1, 365
    do it = 1, nt
      i = (kday-1)*nt + it
      ca_it (j) = ac + bc * float (i)
      !ca_it (j) = ca (kyr)
      !write (*,'(5i8,2f12.4)') kyr,kday,it,i,j,ca(kyr),ca_it(j)
      j = j + 1
    end do
  end do
  !stop
end do
! For diagostic of CO2.
open (20,file='ca_it.txt',status='unknown')
! '5' is default simulation with full mode.
if (sim == 5) then
  open (21,file='o5.txt',status='unknown')
  calib_G = one
  calib_A = one
endif
if (sim == 6) then
  open (21,file='o6.txt',status='unknown')
  calib_G = 0.4 / (2.0 * 0.329)
  calib_A = one
end if
if (sim == 7) then
  open (21,file='o7.txt',status='unknown')
  calib_G = one
  calib_A = one
  ca_it (:) = ca (1700)
end if
if (sim == 8) then
  open (21,file='o8.txt',status='unknown')
  calib_G = one
  calib_A = 1.36
end if
if (sim == 9) then
  open (21,file='o9.txt',status='unknown')
  calib_G = one
  calib_A = 1.36
  ca_it (:) = ca (1700)
end if
KS = KS_def
if (sim == 10) then
  open (21,file='o10.txt',status='unknown')
  calib_G = 1.75
  calib_A = 0.84
  KS = 1.2 * KS_def
end if
if (sim == 11) then
  open (21,file='o11.txt',status='unknown')
  calib_G = 1.75
  calib_A = 0.84
  KS = 1.2 * KS_def
  ca_it (:) = ca (1700)
end if
if (sim == 12) then
  open (21,file='o12.txt',status='unknown')
  calib_G = one
  calib_A = one
end if
if (sim == 13) then
  open (21,file='o13.txt',status='unknown')
  calib_G = one
  calib_A = one
end if
D = 0.1
M = a * D ** b
S = KS * M
gG = one
fi = 0.5
fa = 0.5
j = 1699*365*nt + 1
do kyr = 1700, 2100
  !****adf
  if (kyr == 1801) then
    D = 0.1
    M = a * D ** b
    S = KS * M
  end if
  if (kyr == 2000) ca_it = ca_it + 165.0 !****adf
  !****adf
  do kday = 1, 365
    do it = 1, nt
      if ((sim == 12) .and. (kyr >= 2000)) then
        gA = fC (ca_it (j)) / 2.0
      else
        gA = fC (ca_it (j))
      end if
      if ((sim == 13) .and. (kyr >= 2000)) then
        gG = one / 2.0
      else
        gG = gG
      end if
      call grow (M        , S        , dS1,dM1)
      call grow (M+dM1/2.0, S+dS1/2.0, dS2,dM2)
      call grow (M+dM2/2.0, S+dS2/2.0, dS3,dM3)
      call grow (M+dM3    , S+dS3    , dS4,dM4)
      dM = (1.0 / 6.0) * (dM1 + 2.0 * dM2 + 2.0 * dM3 + dM4)
      dS = (1.0 / 6.0) * (dS1 + 2.0 * dS2 + 2.0 * dS3 + dS4)
      M = M + dt * dM
      S = S + dt * dS
      if ((sim == 6) .or. (sim == 8) .or. (sim == 9)) then
        Sc = KS
      else
        Sc = S / M
      end if
      fi = KS ** n / (KS ** n + Sc ** n)
      fa = Sc ** n / (KS ** n + Sc ** n)
      D = (M / a) ** (1.0 / b)
      U = gA * fi * calib_A * alpha * pi * (D / 2.0) ** 2
      if ((sim == 8) .or. (sim == 9)) dM = U / 0.474
      write (20,*) j, ca_it (j)
      write (21,'(8f16.7)') float(kyr)+float(kday)/365.0+float(it/(365*nt)), &
                            1.0e3*dM,1.0e3*U,M,D,S/M,fi,fa
      j = j + 1
    end do ! it
  end do ! kday
  !write (*,*) kyr, M, D, KS
  if (kyr==1999) write (*,*) Ca_it(j),D,M,U,dM
  if (kyr==2010) write (*,*) Ca_it(j),D,M,U,dM
end do ! kyr
close (20)
close (21)

contains
function fC (Ca)
implicit none
real, parameter :: KC = 400.0
real :: fC, Ca, Ci
Ci = 0.7 * Ca
fC = Ci / (KC + Ci)
fC = fC / 0.392
end function fC

end program sugar

subroutine grow (Mg, Sg, dSg, dMg)
use params
use vars
implicit none
real :: Mg, Sg, dSg, dMg, Ug, Dg
Dg = (Mg / a) ** (1.0 / b)
if ((sim == 6) .or. (sim == 8) .or. (sim == 9)) then
  Sc = KS
else
  Sc = Sg / Mg
end if
fi = KS ** n / (KS ** n + Sc ** n)
fa = Sc ** n / (KS ** n + Sc ** n)
Ug = gA * fi * calib_A * alpha * pi * (Dg / 2.0) ** 2
dMg = gG * fa * calib_G * beta * (a * b) * Dg ** (b - one)
if ((sim == 8).or. (sim == 9)) dMg = Ug / 0.474
dSg = (Ug - 0.474 * dMg) / 0.421
end subroutine grow
