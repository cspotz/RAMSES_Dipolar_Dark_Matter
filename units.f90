subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use constants, only: Mpc2cm, mH, kB, rhoc
  use cooling_module, only: X
  implicit none

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  ! scale_l converts distance from user units into cm
  ! 1 kpc = 3.08567758 x 10^21 cm
  !il faut entrer la taille du domaine et mettre dans les CI entre 0 et 1
  scale_l = 6.171999999999999D+023

  ! scale_d converts mass density from user units into g/cc
  ! From Msun/kpc^3 to g/cm^3
  ! 1 Msun/kpc^3 = (1.98855e+33 g) / (3.08567758e+21 cm)^3 = 6.76838229e-32 g/cm^3
  !il faut mettre la masse en g et cette routine les converti en Msol
  scale_d = 1.6332441716447662D-025! (3.2858295805171536D42 / (scale_l)**3)

  ! scale_t converts time from user units into seconds
  ! G = 6.674D-8 cm^3 g^-1 s^-2
  ! G = 4.302D-6 kpc Msun^-1 (km/s)^2
  scale_t = 1.0/sqrt(6.674D-8 * scale_d)   ! 1.85e+18

  ! scale_v converts velocity in user units into cm/s
  !1 cm/s= 1.0D+5 km/s
  scale_v = scale_l/scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  !scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  !scale_nH = X/mH * scale_d



end subroutine units
