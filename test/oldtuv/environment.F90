! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \todo rewrite micm_environment module to follow musica style conventions
module micm_environment

  use musica_constants,                only : musica_dk

  type environment_t
    !> Number density of air [# cm-3]
    real(musica_dk) :: number_density_air
    !> Temperature [K]
    real(musica_dk) :: temperature
    !> Pressure [Pa]
    real(musica_dk) :: pressure
    !> Aerosol surface area density \bug what are the units for surface area density?
    real(musica_dk) :: aerosol_surface_area_density(4)
    !> Aerosol diameter \bug what are the units for aerosol diameter?
    real(musica_dk) ::  aerosol_diameter(4)
    !> Water vapor volume mixing ratio [mol mol-1]
    real(musica_dk) :: h2ovmr
    !> Molecular oxygen number density [# cm-3]
    real(musica_dk) :: o2_number_density
  end type environment_t

end module micm_environment
