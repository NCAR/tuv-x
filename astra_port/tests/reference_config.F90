module reference_config
  use iso_c_binding
  implicit none
contains
  subroutine ref_configuration(length,text,n,number,integer_value,flag,values,string_value) bind(c)
    use musica_config
    use musica_string
    integer(c_int), value :: length,n
    character(c_char), intent(in) :: text(length)
    real(c_double), intent(out) :: number,values(n)
    integer(c_int), intent(out) :: integer_value,flag
    character(c_char), intent(out) :: string_value(128)
    type(config_t) :: config
    type(string_t) :: string
    character(len=length) :: input
    real(c_double), allocatable :: array(:)
    logical :: boolean
    integer :: i
    do i=1,length
      input(i:i)=text(i)
    enddo
    config=input
    call config%get('number',number,'reference')
    call config%get('integer',integer_value,'reference')
    call config%get('flag',boolean,'reference')
    flag=merge(1,0,boolean)
    call config%get('values',array,'reference')
    if(size(array)/=n) error stop 'wrong array length'
    values=array
    call config%get('string',string,'reference')
    string_value=c_null_char
    do i=1,min(127,len(string%val_))
      string_value(i)=string%val_(i:i)
    enddo
  end subroutine
end module
