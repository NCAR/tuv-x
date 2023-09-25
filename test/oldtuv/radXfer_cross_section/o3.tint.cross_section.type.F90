! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This o3_tint_cross_section module

!> The o3 temperature interpolation cross_section type and related functions
module micm_radXfer_o3_tint_cross_section_type

  use micm_radXfer_tint_cross_section_type,    only : tint_cross_section_t
  use musica_constants,                        only : musica_dk, musica_ik

  implicit none

  private
  public :: o3_tint_cross_section_t

  integer(musica_ik), parameter :: iONE = 1_musica_ik
  integer(musica_ik), parameter :: iTWO = 2_musica_ik
  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for o3_tint_cross_section
  type, extends(tint_cross_section_t) :: o3_tint_cross_section_t
    real(musica_dk) :: v185(1), v195(1), v345(1)
  contains
    !> Initialize the cross section
    procedure :: initialize
    !> Calculate the cross section
    procedure :: calculate => run
    !> refraction
    procedure :: refraction
  end type o3_tint_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize o3_tint_cross_section_t object
  subroutine initialize( this, config, gridWareHouse, ProfileWareHouse )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use netcdf_util,                     only : netcdf_t
    use photo_utils,                     only : inter2
    use musica_assert,                   only : die_msg
    use micm_grid_warehouse,             only : grid_warehouse_t
    use micm_1d_grid,                    only : base_grid_t
    use micm_Profile_warehouse,     only : Profile_warehouse_t

    !> o3 tint cross section type
    class(o3_tint_cross_section_t), intent(inout) :: this
    !> cross section configuration object
    type(config_t), intent(inout) :: config
    !> The warehouses
    type(grid_warehouse_t), intent(inout)         :: gridWareHouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse

    !> Local variables
    real(musica_dk), parameter  :: refracDensity = 2.45e19_musica_dk
    real(musica_dk), parameter  :: w185 = 185._musica_dk
    real(musica_dk), parameter  :: w195 = 195._musica_dk
    real(musica_dk), parameter  :: w345 = 345._musica_dk

    character(len=*), parameter :: Iam = 'o3 tint cross section initialize: '
    character(len=*), parameter :: Hdr = 'cross_section_'

    integer(musica_ik) :: retcode, nmdlLambda
    integer(musica_ik) :: parmNdx, fileNdx, Ndxl, Ndxu
    integer(musica_ik) :: nParms, nTemps
    real(musica_dk)              :: tmp
    real(musica_dk), allocatable :: refracNdx(:)
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_parameter(:)
    logical :: found, monopos
    character(len=:), allocatable :: msg
    type(netcdf_t), allocatable   :: netcdf_obj
    type(string_t), allocatable   :: netcdfFiles(:)
    class(base_grid_t), pointer :: lambdaGrid
    type(string_t)     :: Handle

    write(*,*) Iam,'entering'

    !> Get model wavelength grids
    Handle = 'Photolysis, wavelength'
    lambdaGrid => gridWareHouse%get_grid( Handle )

    !> get cross section netcdf filespec
    call config%get( 'netcdf files', netcdfFiles, Iam, found=found )

has_netcdf_file: &
    if( found ) then
      allocate( this%cross_section_parms(size(netcdfFiles)) )
file_loop: &
      do fileNdx = iONE,size(this%cross_section_parms)
        allocate( netcdf_obj )
    !> read netcdf cross section parameters
        call netcdf_obj%read_netcdf_file( filespec=netcdfFiles(fileNdx)%to_char(), Hdr=Hdr )
        nParms = size(netcdf_obj%parameters,dim=2)
    !> must have at least one parameter
        if( nParms < iONE ) then
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  array must have 1 or more parameters'
          call die_msg( 400000002, msg )
        endif
    !> refraction index
        refracNdx = this%refraction( netcdf_obj%wavelength, refracDensity )
        netcdf_obj%wavelength = refracNdx * netcdf_obj%wavelength

        associate( Xsection => this%cross_section_parms(fileNdx) )
    !> interpolation temperatures must be in netcdf file
        if( allocated(netcdf_obj%temperature) ) then
          Xsection%temperature = netcdf_obj%temperature
          nTemps = size( Xsection%temperature )
          if( nTemps > iONE ) then
            Xsection%deltaT = Xsection%temperature(iTWO:nParms) - Xsection%temperature(iONE:nParms-iONE)
            monopos = all( Xsection%deltaT > rZERO )        
            if( .not. monopos ) then
              if( any( Xsection%deltaT > rZERO ) ) then
                write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  temperature array not monotonic'
                call die_msg( 400000005, msg )
              endif
              do Ndxl = iONE,nParms/iTWO
                Ndxu = nParms - Ndxl + iONE
                tmp = Xsection%temperature(Ndxl)
                Xsection%temperature(Ndxl) = Xsection%temperature(Ndxu)
                Xsection%temperature(Ndxu) = tmp
                data_parameter = netcdf_obj%parameters(:,Ndxl)
                netcdf_obj%parameters(:,Ndxl) = netcdf_obj%parameters(:,Ndxu)
                netcdf_obj%parameters(:,Ndxu) = data_parameter
              enddo
              Xsection%deltaT = Xsection%temperature(iTWO:nParms) - Xsection%temperature(iONE:nParms-iONE)
            endif
          endif
        else
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),' does not have interpolation temperatures'
          call die_msg( 400000006, msg )
        endif

    !> interpolate from data to model wavelength grid
        if( allocated(netcdf_obj%wavelength) ) then
          if( .not. allocated(Xsection%array) ) then
            allocate(Xsection%array(lambdaGrid%ncells_,nParms))
          endif
          do parmNdx = iONE,nParms
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters(:,parmNdx)
            call this%addpnts( config, data_lambda, data_parameter )
            call inter2(xto=lambdaGrid%edge_, &
                        yto=Xsection%array(:,parmNdx), &
                        xfrom=data_lambda, &
                        yfrom=data_parameter,ierr=retcode)
          enddo
        else
          Xsection%array = netcdf_obj%parameters
        endif
        end associate
        deallocate( netcdf_obj )
      enddo file_loop
      this%v185 = this%refraction( (/ w185 /), refracDensity ) * w185
      this%v195 = this%refraction( (/ w195 /), refracDensity ) * w195
      this%v345 = this%refraction( (/ w345 /), refracDensity ) * w345
      this%cross_section_parms(2)%array(:,4) = this%cross_section_parms(1)%array(:,1)
    else has_netcdf_file
      write(msg,*) Iam//'must have at least one netcdf input file'
      call die_msg( 400000008, msg )
    endif has_netcdf_file

    write(*,*) Iam,'exiting'

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, gridWareHouse, ProfileWareHouse ) result( cross_section )

    use micm_grid_warehouse,             only : grid_warehouse_t
    use micm_1d_grid,                    only : base_grid_t
    use micm_Profile_warehouse,     only : Profile_warehouse_t
    use micm_Profile,               only : base_profile_t
    use musica_string,                   only : string_t

    !> Arguments
    !> o3 tint cross section
    class(o3_tint_cross_section_t), intent(in) :: this
    !> The warehouses
    type(grid_warehouse_t), intent(inout)         :: gridWareHouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse
    !> Calculated cross section
    real(kind=musica_dk), allocatable             :: cross_section(:,:)

    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    character(len=*), parameter :: Iam = 'o3 tint cross section calculate: '

    integer(musica_ik) :: layer
    integer(musica_ik) :: nTemp
    integer(musica_ik) :: fileNdx, tNdx, wNdx
    real(musica_dk)    :: Tadj, Tstar
    real(musica_dk), allocatable  :: wrkCrossSection(:,:)
    class(base_grid_t), pointer :: zGrid
    class(base_grid_t), pointer :: lambdaGrid
    class(base_profile_t), pointer :: mdlTemperature
    type(string_t)     :: Handle

    write(*,*) Iam,'entering'

    Handle = 'Vertical Z'
    zGrid => gridWareHouse%get_grid( Handle )
    Handle = 'Photolysis, wavelength'
    lambdaGrid => gridWareHouse%get_grid( Handle )
    Handle = 'Temperature'
    mdlTemperature => ProfileWareHouse%get_Profile( Handle )

    allocate( wrkCrossSection(lambdaGrid%ncells_,zGrid%ncells_) )
    wrkCrossSection = rZERO

layer_loop: &
    do layer = 1,zGrid%ncells_
lambda_loop: &
    do wNdx = iONE,lambdaGrid%ncells_
!     associate( lambda => lambdaGrid%mid_(wNdx) )
      associate( lambda => lambdaGrid%edge_(wNdx) )
        if( lambda < this%v185(1) ) then
          fileNdx = 3_musica_ik
        elseif( this%v185(1) <= lambda .and. lambda < this%v195(1) ) then
          fileNdx = 4_musica_ik
        elseif( this%v195(1) <= lambda .and. lambda < this%v345(1) ) then
          fileNdx = 2_musica_ik
        else
          wrkCrossSection(wNdx,layer) = wrkCrossSection(wNdx,layer) + this%cross_section_parms(1)%array(wNdx,1)
          cycle lambda_loop
        endif
      end associate
      associate( Temp => this%cross_section_parms(fileNdx)%temperature, wrkXsect => this%cross_section_parms(fileNdx) )
      nTemp = size( Temp )
      Tadj  = min( max( mdlTemperature%mid_val_(layer),Temp(iONE) ),Temp(nTemp) )
      do tNdx = iTWO,nTemp 
        if( Tadj <= Temp(tNdx) ) then
          exit
        endif
      enddo
      tNdx = tNdx - iONE
      Tstar = (Tadj - Temp(tNdx))/wrkXsect%deltaT(tNdx)
      wrkCrossSection(wNdx,layer) = wrkCrossSection(wNdx,layer) &
                            + wrkXsect%array(wNdx,tNdx) &
                            + Tstar * (wrkXsect%array(wNdx,tNdx+1) - wrkXsect%array(wNdx,tNdx))
      end associate
    enddo lambda_loop
    enddo layer_loop

    cross_section = transpose( wrkCrossSection )

    write(*,*) Iam,'exiting'

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function refraction( this, wavelength, atmDensity ) result( refrac )

    !> o3 tint cross section
    class(o3_tint_cross_section_t), intent(in) :: this
    !> input vacuum wavelength, nm and air density, molecule cm-3
    real(musica_dk), intent(in)                :: atmDensity
    real(musica_dk), intent(in)                :: wavelength(:)

    !> result
    real(musica_dk)                            :: refrac(size(wavelength))

    real(musica_dk), parameter :: rONE = 1.0_musica_dk
    real(musica_dk), parameter :: divisor = 2.69e19_musica_dk * 273.15_musica_dk/288.15_musica_dk

    integer(musica_ik) :: wNdx
! output refractive index for standard air
! (dry air at 15 deg. C, 101.325 kPa, 0.03% CO2)
    real(musica_dk) :: wrk, sig, sigsq

! from CRC Handbook, originally from Edlen, B., Metrologia, 2, 71, 1966.
! valid from 200 nm to 2000 nm
! beyond this range, use constant value

    do wNdx = iONE,size(wavelength)
      wrk = min( max( wavelength(wNdx),200._musica_dk ),2000._musica_dk ) 
      sig = 1.e3_musica_dk/wrk
      sigsq = sig * sig

      wrk = 8342.13_musica_dk + 2406030._musica_dk/(130._musica_dk - sigsq) &
          + 15997._musica_dk/(38.9_musica_dk - sigsq)

! adjust to local air density

      wrk = wrk * atmDensity/divisor

! index of refraction:

      refrac(wNdx) = rONE + 1.e-8_musica_dk * wrk
    enddo

    end function refraction

end module micm_radXfer_o3_tint_cross_section_type
