subroutine f1(x, y, z)
    real, intent(in) :: x,y
    real, intent(out) :: z
    z = sin(x+y)
end subroutine f1
