! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Time repeated columns; optionally write all fields and rates for comparison
! with another build. Run from a directory containing data/.
program column_sweep
  use iso_fortran_env, only : int32, int64
  use musica_config, only : config_t
  use musica_constants, only : dk => musica_dk
  use musica_mpi, only : musica_mpi_init, musica_mpi_finalize
  use tuvx_core, only : core_t
  use tuvx_grid, only : grid_t
  use tuvx_solver, only : radiation_field_t
  use tuvx_profile_warehouse, only : profile_warehouse_t
  use tuvx_profile, only : profile_t
  implicit none
  class(core_t), pointer :: model
  class(grid_t), pointer :: grid
  type(profile_warehouse_t), pointer :: profiles
  type(config_t) :: config
  type(radiation_field_t) :: field
  type profile_snapshot_t
    class(profile_t), pointer :: value => null()
  end type
  type(profile_snapshot_t) :: snapshots(4)
  character(len=11), parameter :: names(*) = [character(len=11) :: 'temperature','air','O2','O3']
  character(len=1024) :: path, arg, output
  integer :: columns, c, nz, nw, np, nh, nd, unit, i, j
  integer(int64) :: start, finish, rate, ticks
  real(dk), allocatable :: photo(:,:), heat(:,:), dose(:,:)
  real(dk) :: zenith, factor, offset

  call musica_mpi_init()
  call get_command_argument(1,path)
  if( len_trim(path) == 0 ) then
    print *, 'Usage: tuvx_column_sweep CONFIG [COLUMNS] [OUTPUT]'
    stop 2
  end if
  columns = 300
  call get_command_argument(2,arg)
  if( len_trim(arg) > 0 ) read(arg,*) columns
  if( columns < 1 ) stop 2
  call get_command_argument(3,output)
  call config%from_file(trim(path))
  model => core_t(config)
  grid => model%get_grid('height','km')
  nz = grid%size()+1
  deallocate(grid)
  grid => model%get_grid('wavelength','nm')
  nw = grid%size()
  deallocate(grid)
  np = model%number_of_photolysis_reactions()
  nh = model%number_of_heating_rates()
  nd = model%number_of_dose_rates()
  allocate(photo(nz,np),heat(nz,nh),dose(nz,nd))
  ! Snapshot the model profiles, then vary their values between runs. This
  ! exercises column-dependent spectra without changing the configuration.
  profiles => model%get_profile_warehouse()
  do i = 1, size(names)
    if( i == 1 ) then
      snapshots(i)%value => profiles%get_profile(trim(names(i)),'K')
    else if( profiles%exists(trim(names(i)),'molecule cm-3') ) then
      snapshots(i)%value => profiles%get_profile(trim(names(i)),'molecule cm-3')
    end if
  end do
  if( len_trim(output) > 0 ) then
    open(newunit=unit,file=trim(output),access='stream',form='unformatted',status='replace')
    write(unit) int([columns,nz,nw,np,nh,nd],int32)
  end if
  ticks = 0
  call model%run(0.0_dk,1.0_dk,photo,dose,heat)
  do c = 0, columns-1
    do i = 1, size(names)
      if( .not. associated(snapshots(i)%value) ) cycle
      do j = 1, size(profiles%profiles_)
        associate(p => profiles%profiles_(j)%val_, base => snapshots(i)%value)
          if( p%handle_ /= names(i) ) cycle
          factor = 1.0_dk + 0.05_dk*sin(real(c,dk))
          if( i == 4 ) factor = 1.0_dk + 0.3_dk*sin(real(c,dk))
          offset = 10.0_dk*sin(real(c,dk))
          if( i == 1 ) then
            p%edge_val_ = base%edge_val_+offset
            p%mid_val_ = base%mid_val_+offset
          else
            if( allocated(base%edge_val_) ) p%edge_val_ = base%edge_val_*factor
            if( allocated(base%mid_val_) ) p%mid_val_ = base%mid_val_*factor
            if( allocated(base%layer_dens_) ) p%layer_dens_ = base%layer_dens_*factor
            if( allocated(base%exo_layer_dens_) ) p%exo_layer_dens_ = base%exo_layer_dens_*factor
            if( allocated(base%burden_dens_) ) p%burden_dens_ = base%burden_dens_*factor
          end if
        end associate
      end do
    end do
    zenith = 100.0_dk*real(c,dk)/real(max(1,columns-1),dk)
    call system_clock(start,rate)
    call model%run(zenith,1.0_dk+0.017_dk*sin(real(c,dk)),photo,dose,heat)
    call system_clock(finish)
    ticks = ticks+finish-start
    if( len_trim(output) > 0 ) then
      field = model%get_radiation_field()
      write(unit) field%fdr_,field%fdn_,field%fup_,field%edr_,field%edn_,field%eup_,photo,heat,dose
    end if
  end do
  print '(a,f12.2)', 'microseconds per column: ',1.e6_dk*real(ticks,dk)/real(rate,dk)/columns
  if( len_trim(output) > 0 ) close(unit)
  do i = 1, size(names)
    if( associated(snapshots(i)%value) ) deallocate(snapshots(i)%value)
  end do
  deallocate(model)
  call musica_mpi_finalize()
end program column_sweep
