module nquad_module
implicit none
integer, parameter :: dp = kind(0.d0)

integer :: max_depth

real(db), dimension(:), allocatable :: x

function nquad(func,ranges,args,opts)
  real(db), external :: func
  real(db), dimension(:), external :: ranges
  real(db), dimension(:), optional :: args
  integer, dimension(:) :: opts
  

  ! Initialize the module for integration
  max_depth = size(ranges)
  allocate( x(max_depth) )

  deallocate( x )

end function nquad

function f6(x6)

end function f6

function f5(x5)

end function f5

function f4(x4)

end function f4

function f3(x3)

end function f3

function f2(x2)

end function f2

function f1(x1)

end function f1

function f0(x0)
  return func
end function f0

end module nquad_module
