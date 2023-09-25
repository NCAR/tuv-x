   module debug

   use musica_constants, only : dk => musica_dk

   implicit none

   interface diagout
     module procedure :: diagnostic_1d
     module procedure :: diagnostic_1d_dk
     module procedure :: diagnostic_2d
     module procedure :: diagnostic_2d_dk
   end interface diagout

   contains

   subroutine diagnostic_1d( filename, variable )
   
   character(len=*), intent(in) :: filename
   real(4), intent(in)          :: variable(:)

   integer :: ios

   write(*,*) 'diagnostic_1d: entering'

   open(unit=44,file='odat/OUTPUTS/'//filename,form='unformatted',iostat=ios)
   if( ios /= 0 ) then
     write(*,*) 'diagnostic_1d: failed to open ',filename,'; error = ',ios
     stop 'OpnErr'
   endif
   write(unit=44,iostat=ios) variable
   if( ios /= 0 ) then
     write(*,*) 'diagnostic_1d: failed to write ',filename,'; error = ',ios
     stop 'OpnErr'
   endif

   write(*,*) 'diagnostic_1d: exiting'

   end subroutine diagnostic_1d

   subroutine diagnostic_1d_dk( filename, variable )
   
   character(len=*), intent(in) :: filename
   real(dk), intent(in)         :: variable(:)

   integer :: ios

   write(*,*) 'diagnostic_1d_dk: entering'

   open(unit=44,file='odat/OUTPUTS/'//filename,form='unformatted',iostat=ios)
   if( ios /= 0 ) then
     write(*,*) 'diagnostic_1d: failed to open ',filename,'; error = ',ios
     stop 'OpnErr'
   endif
   write(unit=44,iostat=ios) variable
   if( ios /= 0 ) then
     write(*,*) 'diagnostic_1d: failed to write ',filename,'; error = ',ios
     stop 'OpnErr'
   endif

   write(*,*) 'diagnostic_1d_dk: exiting'

   end subroutine diagnostic_1d_dk

   subroutine diagnostic_2d( filename, variable )
   
   character(len=*), intent(in) :: filename
   real(4), intent(in)          :: variable(:,:)

   integer :: ios

   open(unit=44,file='odat/OUTPUTS/'//filename,form='unformatted',iostat=ios)
   if( ios /= 0 ) then
     write(*,*) 'diagnostic_2d: failed to open ',filename,'; error = ',ios
     stop 'OpnErr'
   endif
   write(unit=44,iostat=ios) variable
   if( ios /= 0 ) then
     write(*,*) 'diagnostic_2d: failed to write ',filename,'; error = ',ios
     stop 'OpnErr'
   endif

   end subroutine diagnostic_2d

   subroutine diagnostic_2d_dk( filename, variable )
   
   character(len=*), intent(in) :: filename
   real(dk), intent(in)         :: variable(:,:)

   integer :: ios

   open(unit=44,file='odat/OUTPUTS/'//filename,form='unformatted',iostat=ios)
   if( ios /= 0 ) then
     write(*,*) 'diagnostic_2d: failed to open ',filename,'; error = ',ios
     stop 'OpnErr'
   endif
   write(unit=44,iostat=ios) variable
   if( ios /= 0 ) then
     write(*,*) 'diagnostic_2d: failed to write ',filename,'; error = ',ios
     stop 'OpnErr'
   endif

   end subroutine diagnostic_2d_dk

   end module debug
