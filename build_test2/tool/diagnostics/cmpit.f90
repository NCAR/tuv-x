
   program cmpit

   implicit none

   real, parameter    :: rZERO = 0.0
   integer, parameter :: unit = 20
   real, parameter    :: thresh = .001
!  real, parameter    :: thresh = 1.0
!  real, parameter    :: minvalue = 1.e2
   real, parameter    :: minvalue = rZERO

   integer :: ios, m, filndx, nOpts, nArrays
   integer :: delimNdx, sLower, arrNdx
   integer :: arrayDims(2)
   integer :: maxdims(2), mindims(2)
   real :: maxValue
   real :: accum, rms_error
   real, allocatable :: array(:,:,:,:)
   real, allocatable :: reldiff(:,:)
   character(len=1), allocatable :: visual(:,:)
   character(len=1), allocatable :: marker(:)
   character(len=1), allocatable :: srb(:)
   character(len=50),allocatable :: header(:,:)
   character(len=64), allocatable :: command_tokens(:)
   character(len=64), allocatable :: command_keys(:)
   character(len=64), allocatable :: command_vals(:)
   character(len=256) :: filespec(2)

   logical :: doMask
   logical, allocatable :: mask(:,:)
   logical, allocatable :: masklasrb(:,:)
   logical, allocatable :: combomask(:,:)

   nOpts = command_argument_count()
   if( nOpts < 2 ) then
     write(*,*) 'cmpit: must have two files to compare'
     stop 3
   endif

   allocate( command_tokens(nOpts), command_keys(nOpts), command_vals(nOpts) )

option_loop: &
   do m = 1,nOpts
     call get_command_argument( m,command_tokens(m) )
     delimNdx = index( trim(command_tokens(m)),'=' )
     if( delimNdx == 0 ) then
       write(*,*) 'cmpit: command option must be key=val form'
       stop 3
     endif
     command_keys(m) = command_tokens(m)(:delimNdx-1)
     command_vals(m) = command_tokens(m)(delimNdx+1:)
     select case( command_keys(m) )
       case( 'file1' )
         filespec(1) = trim(command_vals(m))
       case( 'file2' )
         filespec(2) = trim(command_vals(m))
       case( 'dims' )
         nArrays  = 1
         delimNdx = index( trim(command_vals(m)),',' )
         read(command_vals(m)(:delimNdx-1),*) arrayDims(1)
         sLower = delimNdx+1
         delimNdx = index( trim(command_vals(m)),'x' )
         if( delimNdx == 0 ) then
           read(command_vals(m)(sLower:),*) arrayDims(2)
         else
           read(command_vals(m)(sLower:delimNdx-1),*) arrayDims(2)
           read(command_vals(m)(delimNdx+1:),*) nArrays
         endif
         if( any(arrayDims < 1) .or. any(arrayDims > 9999) ) then
           write(*,*) 'cmpit: array dimensions must be > 0 and < 10000'
           stop 3
         endif
         if( nArrays < 1 .or. nArrays > 999 ) then
           write(*,*) 'cmpit: 0 < number arrys < 999'
           stop 3
         endif
       case default
         write(*,*) 'cmpit: command option keys must be file1,file2,dims'
         stop 3
     end select
   enddo option_loop

   write(*,*) ' '
   write(*,*) 'cmpit: filespec #1 = ',trim(filespec(1))
   write(*,*) 'cmpit: filespec #2 = ',trim(filespec(2))
   write(*,*) 'cmpit: array dimensions = ',arrayDims
   write(*,*) 'cmpit: number arrays    = ',nArrays

   allocate( header(nArrays,2) )
   allocate( array(arrayDims(1),arrayDims(2),2,nArrays) )
   allocate( reldiff(arrayDims(1),arrayDims(2)) )
   allocate( visual(arrayDims(1),arrayDims(2)) )
   allocate( marker(arrayDims(2)) )
   allocate( srb(arrayDims(2)) )
   allocate( mask(arrayDims(1),arrayDims(2)) )
   allocate( masklasrb(arrayDims(1),arrayDims(2)) )
   allocate( combomask(arrayDims(1),arrayDims(2)) )

   file_loop: do filndx = 1,2

     open(unit=33,file=trim(filespec(filndx)),form='unformatted',iostat=ios)
     if( ios /= 0 ) then
       write(*,*) 'cmpit: can not locate ',trim(filespec(filndx))
       stop 3
     endif
     do arrNdx = 1,nArrays
       read(unit=33,iostat=ios) header(arrNdx,filndx)
       if( ios /= 0 ) then
         write(*,*) 'cmpit: can not read ',trim(filespec(filndx))
         stop 3
       endif
       read(unit=33,iostat=ios) array(:,:,filndx,arrNdx)
       if( ios /= 0 ) then
         write(*,*) 'cmpit: can not read ',trim(filespec(filndx))
         stop 3
       endif
     enddo
     close(unit=33)

   enddo file_loop

   write(*,*) 'cmpit: input array #1'
   write(*,'(1p10g15.7)') array(:,:,1,1)
   write(*,*) 'cmpit: input array #2'
   write(*,'(1p10g15.7)') array(:,:,2,1)

   array_loop: do arrNdx = 1,nArrays
     write(*,*) '--------------------------------------------------'
     write(*,'('' Analysis for array '',a)') header(arrNdx,1)
     write(*,*) '--------------------------------------------------'
     if( maxval(array(:,:,:,arrNdx)) == rZERO ) then
       write(*,*) ' '
       write(*,*) 'cmpit: both arrays are zero'
     else
       mask = .true.
       reldiff = array(:,:,2,arrNdx) - array(:,:,1,arrNdx)
       where( array(:,:,1,arrNdx) /= rZERO )
         where( array(:,:,1,arrNdx) > minvalue )
           reldiff = 100.*reldiff/array(:,:,1,arrNdx)
         elsewhere
           reldiff = rZERO
         endwhere
       elsewhere
         where( array(:,:,2,arrNdx) == rZERO )
           reldiff = rZERO
         elsewhere
           reldiff = -100.
           mask = .false.
         endwhere
       endwhere

       maxvalue = maxval( array(:,:,1,arrNdx) )
       masklasrb = .true.
!  masklasrb(:,2)     = .false.
!  masklasrb(:,21:36) = .false.

       combomask = mask .and. masklasrb

       maxdims = maxloc( abs(reldiff),mask=combomask )
       mindims = minloc( abs(reldiff),mask=combomask.and.reldiff/=0. )

       write(*,*) ' '
       write(*,*) 'there are ',count(.not.mask),' with differences from zero'
       write(*,*) 'there are ',count(mask),' with valid differences'
       write(*,*) ' '
       write(*,*) 'minimum % diff = ',minval(abs(reldiff),mask=combomask .and. reldiff/=0.)
       write(*,*) 'minimum % diff is @ ',mindims
       write(*,*) 'actual values at min: ',array(mindims(1),mindims(2),1,arrNdx), &
                                           array(mindims(1),mindims(2),2,arrNdx)
       write(*,*) ' '
       write(*,*) 'maximum % diff = ',maxval(abs(reldiff),mask=combomask)
       write(*,*) '% maximum value  = ',100.*array(maxdims(1),maxdims(2),1,arrNdx)/maxValue
       write(*,*) 'maximum % diff is @ ',maxdims
       write(*,*) 'actual values at max: ',array(maxdims(1),maxdims(2),1,arrNdx), &
                                           array(maxdims(1),maxdims(2),2,arrNdx)
       accum = sum(reldiff*reldiff,mask=combomask)
       rms_error = sqrt( accum/real(count(combomask)) )

       write(*,*) ' '
       write(*,*) 'rms % diff = ',rms_error

       write(*,*) ' '
       write(*,*) '% within ',thresh,'% relative error = ',100.*(count(abs(reldiff) <= thresh))/real(count(combomask))

       write(*,*) 'cmpit: % difference'
       write(*,'(1p10g15.7)') reldiff
     endif
   enddo array_loop

   doMask = .false.
   if( doMask ) then
     write(*,*) ' '
     write(*,*) '.01% mask'

     marker = 'c'
     srb = ' '
     srb(21:36) = 'S'
     visual = ' '
     where( combomask .and. abs(reldiff) > thresh )
       visual = 'X'
     endwhere

     write(*,'(3x,200a)') srb
     write(*,'(3x,200a)') marker
     do m = 1,arrayDims(1)
       write(*,'(i3.3,200a)') m,visual(m,:)
     enddo
     write(*,'(3x,200a)') marker
     write(*,'(3x,200a)') srb
   endif

   end program cmpit
