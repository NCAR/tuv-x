! Original LINPACK factorization, condition estimate, and both solve directions.
module reference_linalg
  use iso_c_binding
  implicit none
contains
  subroutine ref_lu(layout,n,ld,ml,mu,a,pivots,rcond,z,rhs,trans) bind(c)
    use tuvx_linear_algebra_linpack
    integer(c_int), value :: layout,n,ld,ml,mu,trans
    real(c_double), intent(inout) :: a(ld,n),rhs(n)
    real(c_double), intent(out) :: rcond,z(n)
    integer(c_size_t), intent(out) :: pivots(n)
    integer :: ipvt(n)
    type(linear_algebra_linpack_t) :: algebra
    if(layout==0)then
      call algebra%SGECO(a,n,ipvt,rcond,z)
      call algebra%SGESL(a,n,ipvt,rhs,trans)
    else
      call algebra%SGBCO(a,n,ml,mu,ipvt,rcond,z)
      call algebra%SGBSL(a,n,ml,mu,ipvt,rhs,trans)
    endif
    pivots=int(ipvt-1,c_size_t)
  end subroutine
end module
