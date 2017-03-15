      subroutine cube(n, x)

      integer n
      double precision x
      dimension x(n)
      integer i

      do 100 i = 1, n
        x(i) = i ** 3
        !write(*,*) x(i)
 100  continue
      end

