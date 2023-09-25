      module la_srb_type

      use musica_constants, only : ik => musica_ik, dk => musica_dk, lk => musica_lk

      implicit none

      private
      public :: la_srb_t

      integer(ik), parameter :: nPoly = 20_ik

      integer(ik), parameter :: kla = 2_ik
      integer(ik), parameter :: nla = kla - 1_ik
      real(dk), parameter :: wlla(kla) = (/ 121.4_dk, 121.9_dk /)

      integer(ik), parameter :: ksrb = 18_ik
      integer(ik), parameter :: nsrb = ksrb - 1_ik
      real(dk), parameter :: wlsrb(ksrb) = &
         (/175.4_dk, 177.0_dk, 178.6_dk, 180.2_dk, 181.8_dk, 183.5_dk, 185.2_dk, 186.9_dk, &
           188.7_dk, 190.5_dk, 192.3_dk, 194.2_dk, 196.1_dk, 198.0_dk, 200.0_dk, 202.0_dk, &
           204.1_dk, 206.2_dk/)

      integer(ik), parameter :: iONE = 1_ik
      real(dk), parameter :: rZERO = 0.0_dk
      real(dk), parameter :: rONE  = 1.0_dk
      real(dk), parameter :: rTWO  = 2.0_dk
      real(dk), parameter :: rTEN  = 10.0_dk
      real(dk), parameter :: precis  = 1.e-7_dk
      real(dk), parameter :: largest = 1.e36_dk

      type :: la_srb_t
        integer(ik) :: ila, isrb
        logical(lk) :: has_la, has_srb, has_la_srb
        real(dk)    :: AC(nPoly,nsrb)
        real(dk)    :: BC(nPoly,nsrb) ! Chebyshev polynomial coeffs
      contains
        procedure :: initialize  
        procedure :: calculate_OD => la_srb_OD
        procedure :: calculate_xs => la_srb_xs
        procedure, private :: lymana_OD, lymana_xs
        procedure, private :: schum_OD, schum_xs
        procedure, private :: init_srb_xs
        procedure, private :: effxs
        procedure, private :: calc_params
        procedure, private :: chebev
      end type la_srb_t

      contains

      subroutine initialize( this, gridWareHouse )
!----------------------------------------------------------------------
! check that the user wavelength grid, WL(IW), is compatible 
! with the wavelengths for the parameterizations of the Lyman-alpha and SRB.
! Also compute and save corresponding grid indices (ILA, ISRB)
!----------------------------------------------------------------------
      use micm_grid_warehouse,  only : grid_warehouse_t
      use micm_1d_grid,         only : base_grid_t
      use musica_assert,        only : die_msg
      use musica_string,        only : string_t

      !> Arguments
      class(la_srb_t), intent(inout)        :: this
      type(grid_warehouse_t), intent(inout) :: gridWareHouse

      !> Local variables
      character(len=*), parameter :: Iam = 'la_srb initialize: '

      integer(ik) :: iw, nw
      type(string_t)                     :: Handle
      class(base_grid_t), pointer      :: lambdaGrid

      write(*,*) ' '
      write(*,*) Iam // 'entering'

      Handle = 'Photolysis, wavelength' ; lambdaGrid => gridWareHouse%get_grid( Handle )

      !> Are la and srb grids fully "inside" the model grid?
      this%has_la  = lambdaGrid%edge_(iONE) <= wlla(iONE) .and. lambdaGrid%edge_(lambdaGrid%ncells_+iONE) >= wlla(kla)
      this%has_srb = lambdaGrid%edge_(iONE) <= wlsrb(iONE) .and. lambdaGrid%edge_(lambdaGrid%ncells_+iONE) >= wlsrb(ksrb)
      this%has_la_srb = this%has_la .or. this%has_srb
has_la_srb: &
      if( this%has_la_srb ) then
        nw = lambdaGrid%ncells_ + iONE
        if( this%has_la ) then
! locate Lyman-alpha wavelengths on grid
          this%ila = 0
          do iw = iONE, nw
            if(ABS(lambdaGrid%edge_(iw) - wlla(iONE)) < rTEN*precis) then
              this%ila = iw
              EXIT
            endif
          enddo
! check Lyman-alpha wavelength grid
          if(this%ila == 0) then
            write(*,*) 'For wavelengths below 205.8 nm, only the'
            write(*,*) 'pre-specified wavelength grid is permitted'
            write(*,*) 'Use nwint=-156, or edit subroutine gridw.f'
            call die_msg( 20001,' Lyman alpha grid mis-match - 1' )
          endif
          do iw = 2_ik, nla + iONE
            if(ABS(lambdaGrid%edge_(this%ila + iw - iONE) - wlla(iw)) > rTEN*precis) then
              call die_msg( 20002,' Lyman alpha grid mis-match - 2' )
            endif
          enddo
        endif
        if( this%has_srb ) then
! locate Schumann-Runge wavelengths on grid
          this%isrb = 0
          do iw = iONE, nw
            if(ABS(lambdaGrid%edge_(iw) - wlsrb(iONE)) < rTEN*precis) then
              this%isrb = iw
              EXIT
            endif
          enddo
! check Schumann-Runge wavelength grid
          if(this%isrb == 0) then
            write(*,*) 'For wavelengths below 205.8 nm, only the'
            write(*,*) 'pre-specified wavelength grid is permitted'
            write(*,*) 'Use nwint=-156, or edit subroutine gridw.f'
            call die_msg( 20003,' SRB grid mis-match - 1' )
          endif
          do iw = 2_ik, nsrb + iONE
            if(ABS(lambdaGrid%edge_(this%isrb + iw - iONE) - wlsrb(iw)) > rTEN* precis) then
              call die_msg( 20002,' SRB grid mismatch - w' )
            endif
          enddo
!> Loads Chebyshev polynomial Coeff.
          call this%init_srb_xs
        endif
      endif has_la_srb

      write(*,*) ' '
      write(*,*) Iam // 'exiting'

      end subroutine initialize

      subroutine la_srb_OD( this,gridWareHouse,ProfileWareHouse,Airvcol,Airscol,dto2 )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Compute equivalent optical depths for O2 absorption, and O2 effective    =*
!=  absorption cross sections, parameterized in the Lyman-alpha and SR bands =*
!-----------------------------------------------------------------------------* 
!=  parameterS:                                                              =*
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
!=            grid                                                           =*
!=  Z       - REAL, specified altitude working grid (km)                  (I)=*
!=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
!=            wavelength grid                                                =*
!=  WL      - REAL, vector of lxower limits of wavelength intervals in    (I)=*
!=            working wavelength grid                                        =*
!=  CZ      - REAL, number of air molecules per cm^2 at each specified    (I)=*
!=            altitude layer                                                 =*
!=  ZEN     - REAL, solar zenith angle                                    (I)=*
!=                                                                           =*
!=  O2XS1   - REAL, O2 cross section from rdo2xs                          (I)=*
!=                                                                           =*
!=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
!=            vertical layer at each specified wavelength                    =*
!=  O2XS    - REAL, molecular absorption cross section in SR bands at     (O)=*
!=            each specified altitude and wavelength.  Includes Herzberg     =*
!=            continuum.                                                     =*
!-----------------------------------------------------------------------------*

    use micm_Profile_warehouse,   only : Profile_warehouse_t
    use micm_Profile,             only : base_profile_t
    use micm_grid_warehouse,      only : grid_warehouse_t
    use micm_1d_grid,             only : base_grid_t
    use musica_string,            only : string_t

    !> Arguments
    class(la_srb_t), intent(inout)  :: this
    !> Vertical Profile warehouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse
    !> grid warehouse
    type(grid_warehouse_t), intent(inout)    :: gridWareHouse

    real(dk), intent(in)    :: Airvcol(:), Airscol(:)
    real(dk), intent(inout) :: dto2(:,:)

    !> Local variables
    real(dk), parameter         :: o2Vmr = .2095_dk
    character(len=*), parameter :: Iam = 'la_srb OD: '

    integer(ik) :: nz, nzm1, nw, i, iz, iw
    real(dk)    :: secchi(size(Airscol))
    real(dk)    :: o2scol(size(Airscol))
    class(base_grid_t), pointer :: zGrid
    class(base_grid_t), pointer :: lambdaGrid
    class(base_profile_t), pointer :: temperature
    type(string_t) :: Handle

!----------------------------------------------------------------------
! Lyman-alpha variables
! O2 optical depth and equivalent cross section in the Lyman-alpha region
! Wavelengths for Lyman alpha and SRB parameterizations:
!----------------------------------------------------------------------
    real(dk)    :: dto2la(size(Airvcol),nla)
!----------------------------------------------------------------------
! grid on which Koppers' parameterization is defined
! O2 optical depth and equivalent cross section on Koppers' grid
!----------------------------------------------------------------------
    real(dk)    :: dto2k(size(Airvcol),nsrb)

    write(*,*) ' '
    write(*,*) Iam // 'entering'

has_la_srb: &
    if( this%has_la_srb ) then
!-----------------------------------------------------------------------------
!> get specific grids and vertical profiles
!-----------------------------------------------------------------------------
      Handle = 'Vertical Z'
      zGrid => gridWareHouse%get_grid( Handle )
      Handle = 'Photolysis, wavelength'
      lambdaGrid => gridWareHouse%get_grid( Handle )

      Handle = 'Temperature'
      temperature => ProfileWareHouse%get_Profile( Handle )

      nw   = lambdaGrid%ncells_ + iONE
      nzm1 = zGrid%ncells_
      nz   = nzm1 +  iONE
      !>----------------------------------------------------------------------
      !> O2 slant column
      !>----------------------------------------------------------------------
      o2scol(:) = o2Vmr * Airscol(:)
      !>----------------------------------------------------------------------
      !> Effective secant of solar zenith angle.  
      !> Use 2.0 if no direct sun (value for isotropic radiation)
      !> For nz, use value at nz-1
      !>----------------------------------------------------------------------
      where( Airscol(1:nzm1) <= .1_dk*largest )
        secchi(1:nzm1) = Airscol(1:nzm1)/Airvcol(1:nzm1)
      elsewhere
        secchi(1:nzm1) = rTWO
      endwhere
      secchi(nz) = secchi(nzm1)
!---------------------------------------------------------------------
! Lyman-Alpha parameterization, output values of O2 optical depth
! and O2 effective (equivalent) cross section
!----------------------------------------------------------------------
      if( this%has_la ) then
        call this%lymana_OD( o2scol,secchi,dto2la )
        do iw = this%ila, this%ila + nla - iONE
          dto2(:,iw) = dto2la(:,iw-this%ila+iONE)
        enddo
      endif
!------------------------------------------------------------------------------
! Koppers' parameterization of the SR bands, output values of O2
! optical depth and O2 equivalent cross section 
!------------------------------------------------------------------------------
      if( this%has_srb ) then
        call this%schum_OD( o2scol,temperature%edge_val_,secchi,dto2k )
        do iw = this%isrb, this%isrb + nsrb - iONE
          dto2(:,iw) = dto2k(:,iw-this%isrb+iONE)
        enddo
      endif
    endif has_la_srb

    write(*,*) ' '
    write(*,*) Iam // 'exiting'

    end subroutine la_srb_OD

    subroutine la_srb_xs( this,gridWareHouse,ProfileWareHouse,Airvcol,Airscol,o2xs )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Compute equivalent optical depths for O2 absorption, and O2 effective    =*
!=  absorption cross sections, parameterized in the Lyman-alpha and SR bands =*
!-----------------------------------------------------------------------------* 
!=  parameterS:                                                              =*
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
!=            grid                                                           =*
!=  Z       - REAL, specified altitude working grid (km)                  (I)=*
!=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
!=            wavelength grid                                                =*
!=  WL      - REAL, vector of lxower limits of wavelength intervals in    (I)=*
!=            working wavelength grid                                        =*
!=  CZ      - REAL, number of air molecules per cm^2 at each specified    (I)=*
!=            altitude layer                                                 =*
!=  ZEN     - REAL, solar zenith angle                                    (I)=*
!=                                                                           =*
!=  O2XS1   - REAL, O2 cross section from rdo2xs                          (I)=*
!=                                                                           =*
!=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
!=            vertical layer at each specified wavelength                    =*
!=  O2XS    - REAL, molecular absorption cross section in SR bands at     (O)=*
!=            each specified altitude and wavelength.  Includes Herzberg     =*
!=            continuum.                                                     =*
!-----------------------------------------------------------------------------*

    use micm_Profile_warehouse,   only : Profile_warehouse_t
    use micm_Profile,             only : base_profile_t
    use micm_grid_warehouse,      only : grid_warehouse_t
    use micm_1d_grid,             only : base_grid_t
    use musica_string,            only : string_t

    !> Arguments
    class(la_srb_t), intent(inout)  :: this
    !> Vertical Profile warehouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse
    !> grid warehouse
    type(grid_warehouse_t), intent(inout)    :: gridWareHouse

    real(dk), intent(in)    :: Airvcol(:), Airscol(:)
    real(dk), intent(inout) :: o2xs(:,:)

    !> Local variables
    real(dk), parameter :: o2Vmr = .2095_dk

    integer(ik) :: nz, nzm1, nw, i, iz, iw
    real(dk)    :: secchi(size(Airscol))
    real(dk)    :: o2scol(size(Airscol))
    class(base_grid_t), pointer :: zGrid
    class(base_grid_t), pointer :: lambdaGrid
    class(base_profile_t), pointer :: temperature
    type(string_t) :: Handle

!----------------------------------------------------------------------
! Lyman-alpha variables
! O2 optical depth and equivalent cross section in the Lyman-alpha region
! Wavelengths for Lyman alpha and SRB parameterizations:
!----------------------------------------------------------------------
    real(dk)    :: o2xsla(size(Airscol),nla)
!----------------------------------------------------------------------
! grid on which Koppers' parameterization is defined
! O2 optical depth and equivalent cross section on Koppers' grid
!----------------------------------------------------------------------
    real(dk)    :: o2xsk(size(Airscol),nsrb)

has_la_srb: &
    if( this%has_la_srb ) then
!-----------------------------------------------------------------------------
!> get specific grids and vertical profiles
!-----------------------------------------------------------------------------
      Handle = 'Vertical Z'
      zGrid => gridWareHouse%get_grid( Handle )
      Handle = 'Photolysis, wavelength'
      lambdaGrid => gridWareHouse%get_grid( Handle )

      Handle = 'Temperature'
      temperature => ProfileWareHouse%get_Profile( Handle )

      nw   = lambdaGrid%ncells_ + iONE
      nzm1 = zGrid%ncells_
      nz   = nzm1 +  iONE
      !>----------------------------------------------------------------------
      !> Slant O2 column
      !>----------------------------------------------------------------------
      o2scol(:) = o2Vmr * Airscol(:)
      !>----------------------------------------------------------------------
      !> Effective secant of solar zenith angle.  
      !> Use 2.0 if no direct sun (value for isotropic radiation)
      !> For nz, use value at nz-1
      !>----------------------------------------------------------------------
      where( Airscol(1:nzm1) <= .1_dk*largest )
        secchi(1:nzm1) = Airscol(1:nzm1)/Airvcol(1:nzm1)
      elsewhere
        secchi(1:nzm1) = rTWO
      endwhere
      secchi(nz) = secchi(nzm1)
!---------------------------------------------------------------------
! Lyman-Alpha parameterization, output values of O2 optical depth
! and O2 effective (equivalent) cross section
!----------------------------------------------------------------------
      if( this%has_la ) then
        call this%lymana_xs( o2scol,secchi,o2xsla )
        do iw = this%ila, this%ila + nla - iONE
          o2xs(:,iw) = o2xsla(:,iw-this%ila+iONE)
        enddo
      endif
!------------------------------------------------------------------------------
! Koppers' parameterization of the SR bands, output values of O2
! optical depth and O2 equivalent cross section 
!------------------------------------------------------------------------------
      if( this%has_srb ) then
        call this%schum_xs( o2scol,temperature%edge_val_,secchi,o2xsk )
        do iw = this%isrb, this%isrb + nsrb - iONE
          o2xs(:,iw) = o2xsk(:,iw-this%isrb+iONE)
        enddo
      endif
    endif has_la_srb

    end subroutine la_srb_xs

    subroutine lymana_OD( this,o2col,secchi,dto2la )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate the effective absorption cross section of O2 in the Lyman-Alpha=*
!=  bands and an effective O2 optical depth at all altitudes.  Parameterized =*
!=  after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the  =*
!=  absorption of the solar Lyman-Alpha line, Geophysical Research Letters,  =*
!=  Vol.24, No.21, pp 2659-2662, 1997.                                       =*
!-----------------------------------------------------------------------------*
!=  parameterS:                                                              =*
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
!=            grid                                                           =*
!=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
!=            altitude                                                       =*
!=  DTO2LA  - REAL, optical depth due to O2 absorption at each specified  (O)=*
!=            vertical layer                                                 =*
!=  O2XSLA  - REAL, molecular absorption cross section in LA bands        (O)=*
!-----------------------------------------------------------------------------*

      !> Arguments
      class(la_srb_t), intent(inout)  :: this
      real(dk), intent(in) :: o2col(:)
      real(dk), intent(in) :: secchi(:)

      real(dk), intent(out) :: dto2la(:,:)

      !> Local variables
      real(dk), parameter :: xsmin     = 1.e-20_dk
      real(dk), parameter :: tiny_val  = 1.e-100_dk
      real(dk), parameter :: exp_lim   = 100.e8_dk
      real(dk), parameter :: large_od  = 1000._dk
      real(dk), parameter :: b(3) = (/ 6.8431e-01_dk, 2.29841e-01_dk,  8.65412e-02_dk /)
      real(dk), parameter :: c(3) = (/8.22114e-21_dk, 1.77556e-20_dk,  8.22112e-21_dk /)
      real(dk), parameter :: d(3) = (/ 6.0073e-21_dk, 4.28569e-21_dk,  1.28059e-20_dk /)
      real(dk), parameter :: e(3) = (/ 8.21666e-21_dk, 1.63296e-20_dk,  4.85121e-17_dk /)

      integer(ik) :: nz
      integer(ik) :: iz, i
      real(dk) :: coldens
      real(dk) :: sigma(3), tau(3)
      real(dk) :: rm(size(o2col))

! calculate reduction factors at every altitude

      rm(:)  = rZERO
      nz = size(o2col)
      do iz = iONE, nz
        coldens = o2col(iz)
        sigma   = c * coldens
        where( sigma < exp_lim )
          tau = EXP( -sigma )
        elsewhere
          tau = rZERO
        endwhere
        rm(iz) = dot_product( tau,b )
      enddo

! calculate effective O2 optical depth
      do iz = iONE, nz-iONE
        if( rm(iz) > tiny_val ) then
          if( rm(iz+iONE) > rZERO) then
            dto2la(iz,iONE) = LOG(rm(iz+1)) / secchi(iz+1) - LOG(rm(iz)) / secchi(iz)
          else
            dto2la(iz,iONE) = large_od
          endif
        else
          dto2la(iz,iONE) = large_od
        endif
      enddo

    end subroutine lymana_OD

    subroutine lymana_xs( this,o2col,secchi,o2xsla )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate the effective absorption cross section of O2 in the Lyman-Alpha=*
!=  bands and an effective O2 optical depth at all altitudes.  Parameterized =*
!=  after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the  =*
!=  absorption of the solar Lyman-Alpha line, Geophysical Research Letters,  =*
!=  Vol.24, No.21, pp 2659-2662, 1997.                                       =*
!-----------------------------------------------------------------------------*
!=  parameterS:                                                              =*
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
!=            grid                                                           =*
!=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
!=            altitude                                                       =*
!=  DTO2LA  - REAL, optical depth due to O2 absorption at each specified  (O)=*
!=            vertical layer                                                 =*
!=  O2XSLA  - REAL, molecular absorption cross section in LA bands        (O)=*
!-----------------------------------------------------------------------------*

      !> Arguments
      class(la_srb_t), intent(inout)  :: this
      real(dk), intent(in) :: o2col(:)
      real(dk), intent(in) :: secchi(:)

      real(dk), intent(out) :: o2xsla(:,:)

      !> Local variables
      real(dk), parameter :: xsmin     = 1.e-20_dk
      real(dk), parameter :: tiny_val  = 1.e-100_dk
      real(dk), parameter :: exp_lim   = 100.e8_dk
      real(dk), parameter :: large_od  = 1000._dk
      real(dk), parameter :: b(3) = (/ 6.8431e-01_dk, 2.29841e-01_dk,  8.65412e-02_dk /)
      real(dk), parameter :: c(3) = (/8.22114e-21_dk, 1.77556e-20_dk,  8.22112e-21_dk /)
      real(dk), parameter :: d(3) = (/ 6.0073e-21_dk, 4.28569e-21_dk,  1.28059e-20_dk /)
      real(dk), parameter :: e(3) = (/ 8.21666e-21_dk, 1.63296e-20_dk,  4.85121e-17_dk /)

      integer(ik) :: nz
      integer(ik) :: iz, i
      real(dk) :: coldens
      real(dk) :: sigma(3), tau(3)
      real(dk) :: rm(size(o2col)), ro2(size(o2col))

! calculate reduction factors at every altitude

      rm(:)  = rZERO
      ro2(:) = rZERO
      nz = size(o2col)
      do iz = iONE, nz
        coldens = o2col(iz)
        sigma   = c * coldens
        where( sigma < exp_lim )
          tau = EXP( -sigma )
        elsewhere
          tau = rZERO
        endwhere
        rm(iz) = dot_product( tau,b )

        sigma   = e * coldens
        where( sigma < exp_lim )
          tau = EXP( -sigma ) 
        elsewhere
          tau = rZERO
        endwhere
        ro2(iz) = dot_product( tau,d )
      enddo

! calculate effective O2 cross sections
      do iz = iONE, nz-iONE
        if( rm(iz) > tiny_val ) then
          if( ro2(iz) > tiny_val ) then
            o2xsla(iz,iONE) = ro2(iz)/rm(iz)
          else
            o2xsla(iz,iONE) = xsmin
          endif
        else
          o2xsla(iz,iONE) = xsmin
        endif
      enddo

! do top layer separately
      if(rm(nz) > tiny_val ) then
         o2xsla(nz,1) = ro2(nz)/rm(nz)
      else
         o2xsla(nz,1) = xsmin
      endif

      end subroutine lymana_xs

      subroutine schum_OD( this, o2col, tlev, secchi, dto2k )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate the equivalent absorption cross section of O2 in the SR bands. =*
!=  The algorithm is based on parameterization of G.A. Koppers, and          =*
!=  D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]                         =*
!=  Final values do include effects from the Herzberg continuum.             =*
!-----------------------------------------------------------------------------*
!=  parameterS:                                                              =*
!=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
!=            altitude                                                       =*
!=  TLEV    - tmeperature at each level                                   (I)=*
!=  SECCHI  - ratio of slant to vertical o2 columns                       (I)=*
!=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
!=            vertical layer at each specified wavelength                    =*
!=  O2XSK  - REAL, molecular absorption cross section in SR bands at     (O)=*
!=            each specified wavelength.  Includes Herzberg continuum        =*
!-----------------------------------------------------------------------------*
    
      !> Arguments
      class(la_srb_t), intent(inout)  :: this
      real(dk), intent(in) :: o2col(:)
      real(dk), intent(in) :: tlev(:), secchi(:)

      real(dk), intent(out) :: dto2k(:,:)

      !> Local variables
      real(dk), parameter :: colmin = exp( 38._dk )
      real(dk), parameter :: xslod(nsrb) = &
                (/ 6.2180730E-21_dk, 5.8473627E-22_dk, 5.6996334E-22_dk, &
                   4.5627094E-22_dk, 1.7668250E-22_dk, 1.1178808E-22_dk, &
                   1.2040544E-22_dk, 4.0994668E-23_dk, 1.8450616E-23_dk, &
                   1.5639540E-23_dk, 8.7961075E-24_dk, 7.6475608E-24_dk, &
                   7.6260556E-24_dk, 7.5565696E-24_dk, 7.6334338E-24_dk, &
                   7.4371992E-24_dk, 7.3642966E-24_dk /)

      integer(ik) :: nz
      integer(ik) :: k, kp1, ktop, ktop1, kbot, lambdaNdx

      real(dk) :: X, NORM
      real(dk) :: o2col1(size(o2col))
      real(dk), allocatable :: o2xsk(:,:)

      allocate( o2xsk(size(o2col),nsrb) )
!------------------------------------------
!sm	 Initialize cross sections to values
!sm	 at large optical depth
!------------------------------------------
      do lambdaNdx = iONE, nsrb
        o2xsk(:,lambdaNdx) = xslod(lambdaNdx)
      enddo

!------------------------------------------
!     Calculate cross sections
!sm:  Set smallest O2col = exp(38.) molec cm-2
!sm     to stay in range of parameterization
!sm     given by Koppers et al. at top of atm.
!------------------------------------------
      nz = size(o2col)
      ktop = nz
      kbot = 0_ik
      NORM = rONE/REAL(nz-iONE,dk)
      do k = iONE,nz
        o2col1(k) = MAX(o2col(k),colmin)
        x  = LOG(o2col1(k))
        if (x < 38.0_dk) then
          ktop1 = k - iONE
          ktop  = MIN(ktop1,ktop)
        else if (x > 56.0_dk) then
          kbot = k
        else
          o2xsk(k,:nsrb) = this%effxs( x, tlev(k) )
        endif
      enddo

!------------------------------------------
!  fill in cross section where X is out of range 
!  by repeating edge table values
!------------------------------------------
!sm do not allow kbot = nz to avoid division by zero in no light case.
      if(kbot == nz) then
        kbot = nz - iONE
      endif

      do lambdaNdx = iONE,nsrb
        o2xsk(iONE:kbot,lambdaNdx)    = o2xsk(kbot+iONE,lambdaNdx)
        o2xsk(ktop+iONE:nz,lambdaNdx) = o2xsk(ktop,lambdaNdx)
      enddo
      
!------------------------------------------
!  Calculate incremental optical depth
!------------------------------------------
      do lambdaNdx = iONE,nsrb
        do k = iONE,nz-iONE
          kp1 = k + iONE
!... calculate an optical depth weighted by density
!sm:  put in mean value estimate, if in shade
          if( abs(rONE - o2col1(kp1)/o2col1(k)) <= rTWO*precis ) then
            dto2k(k,lambdaNdx) = o2xsk(kp1,lambdaNdx)*o2col1(kp1)*NORM
          else
            dto2k(k,lambdaNdx) = &
                abs( (o2xsk(kp1,lambdaNdx)*o2col1(kp1) - o2xsk(k,lambdaNdx)*o2col1(k)) &
                / (rONE + LOG(o2xsk(kp1,lambdaNdx)/o2xsk(k,lambdaNdx))/LOG(o2col1(kp1)/o2col1(k))) )
!... change to vertical optical depth
            dto2k(k,lambdaNdx) = rTWO * dto2k(k,lambdaNdx) / (secchi(k) + secchi(kp1))
          endif
        enddo
      enddo 

      end subroutine schum_OD

      subroutine schum_xs( this, o2col, tlev, secchi, o2xsk )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate the equivalent absorption cross section of O2 in the SR bands. =*
!=  The algorithm is based on parameterization of G.A. Koppers, and          =*
!=  D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]                         =*
!=  Final values do include effects from the Herzberg continuum.             =*
!-----------------------------------------------------------------------------*
!=  parameterS:                                                              =*
!=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
!=            altitude                                                       =*
!=  TLEV    - tmeperature at each level                                   (I)=*
!=  SECCHI  - ratio of slant to vertical o2 columns                       (I)=*
!=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
!=            vertical layer at each specified wavelength                    =*
!=  O2XSK  - REAL, molecular absorption cross section in SR bands at     (O)=*
!=            each specified wavelength.  Includes Herzberg continuum        =*
!-----------------------------------------------------------------------------*
    
      !> Arguments
      class(la_srb_t), intent(inout)  :: this
      real(dk), intent(in) :: o2col(:)
      real(dk), intent(in) :: tlev(:), secchi(:)

      real(dk), intent(out) :: o2xsk(:,:)

      !> Local variables
      real(dk), parameter :: colmin = exp( 38._dk )
      real(dk), parameter :: xslod(nsrb) = &
                (/ 6.2180730E-21_dk, 5.8473627E-22_dk, 5.6996334E-22_dk, &
                   4.5627094E-22_dk, 1.7668250E-22_dk, 1.1178808E-22_dk, &
                   1.2040544E-22_dk, 4.0994668E-23_dk, 1.8450616E-23_dk, &
                   1.5639540E-23_dk, 8.7961075E-24_dk, 7.6475608E-24_dk, &
                   7.6260556E-24_dk, 7.5565696E-24_dk, 7.6334338E-24_dk, &
                   7.4371992E-24_dk, 7.3642966E-24_dk /)

      integer(ik) :: nz
      integer(ik) :: i, k, kp1, ktop, ktop1, kbot, lambdaNdx

      real(dk) :: X, NORM
      real(dk) :: o2col1(size(o2col))

!------------------------------------------
!sm	 Initialize cross sections to values
!sm	 at large optical depth
!------------------------------------------
      do lambdaNdx = iONE, nsrb
        o2xsk(:,lambdaNdx) = xslod(lambdaNdx)
      enddo

!------------------------------------------
!     Calculate cross sections
!sm:  Set smallest O2col = exp(38.) molec cm-2
!sm     to stay in range of parameterization
!sm     given by Koppers et al. at top of atm.
!------------------------------------------
      nz = size(o2col)
      ktop = nz
      kbot = 0_ik
      NORM = rONE/REAL(nz-iONE,dk)
      do k = iONE,nz
        o2col1(k) = MAX(o2col(k),colmin)
        x  = LOG(o2col1(k))
        if (x < 38.0_dk) then
          ktop1 = k - 1_ik
          ktop  = MIN(ktop1,ktop)
        else if (x > 56.0_dk) then
          kbot = k
        else
          o2xsk(k,:nsrb) = this%effxs( x, tlev(k) )
        endif
      enddo

!------------------------------------------
!  fill in cross section where X is out of range 
!  by repeating edge table values
!------------------------------------------
!sm do not allow kbot = nz to avoid division by zero in
!   no light case.
       
      if(kbot == nz) then
        kbot = nz - iONE
      endif

      do lambdaNdx = iONE,nsrb
        o2xsk(iONE:kbot,lambdaNdx)    = o2xsk(kbot+iONE,lambdaNdx)
        o2xsk(ktop+iONE:nz,lambdaNdx) = o2xsk(ktop,lambdaNdx)
      enddo

      end subroutine schum_xs

      subroutine init_srb_xs( this )
!-------------------------------------------------------------
!	polynomial coeffs necessary to calculate O2 effective
!       cross-sections
!-------------------------------------------------------------

      use musica_assert, only : die_msg

      class(la_srb_t), intent(inout)  :: this
      integer(ik) :: in_lun             ! file unit number
      integer(ik) :: i, j, istat

      character(len=*), parameter :: filespec = 'odat/DATAE1/O2/effxstex.txt'

      in_lun = 11_ik

      OPEN( UNIT=in_lun, FILE=filespec, FORM='FORMATTED', iostat=istat )
      if( istat /= 0 ) then
        call die_msg( 34056,'failed to open ' // trim(filespec) )
      endif

      READ(in_lun,901)
      do I = iONE,nPoly
        READ(in_lun,903,iostat=istat) (this%AC(I,J), J=1,nsrb)
        if( istat /= 0 ) then
          call die_msg( 34057,'failed to read ' // trim(filespec) )
        endif
      enddo
      READ(in_lun,901)
      do I = iONE,nPoly
        READ(in_lun,903,iostat=istat) (this%BC(I,J), J=1,nsrb)
        if( istat /= 0 ) then
          call die_msg( 34058,'failed to read ' // trim(filespec) )
        endif
      enddo

 901  FORMAT( / )
 903  FORMAT( 17(E23.14,1x))

      CLOSE( in_lun )

      end subroutine init_srb_xs

      function effxs( this, X, T ) result( xs )
!-------------------------------------------------------------
!     Subroutine for evaluating the effective cross section
!     of O2 in the Schumann-Runge bands using parameterization
!     of G.A. Koppers, and D.P. Murtagh [ref. Ann.Geophys., 14
!     68-79, 1996]
!      
!     method:
!     ln(xs) = A(X)[T-220]+B(X)
!     X = log of slant column of O2
!     A,B calculated from Chebyshev polynomial coeffs
!     AC and BC using NR routine chebev.  Assume interval
!     is 38<ln(NO2)<56.
!
!     Revision History:
!
!     drm 2/97  initial coding
!-------------------------------------------------------------

      !> Arguments
      class(la_srb_t), intent(inout)  :: this
      real(dk), intent(in)  :: T, X
      real(dk), allocatable :: XS(:)

      !> Local variables
      real(dk), parameter :: T0 = 220._dk
      REAL(dk)     :: A(nsrb), B(nsrb) 

      call this%calc_params( X, A, B )

      XS = exp( A*(T - T0) + B )

      end function effxs

      subroutine calc_params( this, X, A, B )
!-------------------------------------------------------------
!       calculates coefficients (A,B), used in calculating the
!	effective cross section, for 17 wavelength intervals
!       as a function of log O2 column density (X)
!       Wavelength intervals are defined in WMO1985
!-------------------------------------------------------------

      class(la_srb_t), intent(inout)  :: this
      REAL(dk), intent(in)  :: X
      REAL(dk), intent(out) :: A(:), B(:)

      INTEGER(ik) :: I

!-------------------------------------------------------------
!       call Chebyshev Evaluation routine to calc A and B from
!	set of 20 coeficients for each wavelength
!-------------------------------------------------------------
      do I = 1,size(A)
        A(I) = this%chebev(38.0_dk , 56.0_dk, this%AC(1,I), nPoly, X)
        B(I) = this%chebev(38.0_dk , 56.0_dk, this%BC(1,I), nPoly, X)
      enddo

      end subroutine calc_params

      function chebev(this,a,b,c,m,x) result(poly)
!-------------------------------------------------------------
!     Chebyshev evaluation algorithm
!     See Numerical recipes p193
!-------------------------------------------------------------

      class(la_srb_t), intent(inout)  :: this
      integer(ik), intent(in) :: M
      real(dk), intent(in)    :: A, B, X
      real(dk), intent(in)    :: C(M)

      real(dk) :: poly

      integer(ik) :: J
      real(dk)    :: D,DD,SV,Y,Y2


      if( (X-A)*(X-B) > rZERO ) then
        write(*,*) 'X NOT IN RANGE IN CHEBEV', X
        poly = rZERO
      else
         D  = rZERO
         DD = rZERO
         Y  = (rTWO*X - (A + B))/(B-A)
         Y2 = rTWO*Y
         do J = M,2_ik,-1_ik
           SV = D
           D  = Y2*D-DD + C(J)
           DD = SV
         enddo
         poly = Y*D - DD + 0.5_dk*C(1)
      endif
      
      end function chebev

      end module la_srb_type
