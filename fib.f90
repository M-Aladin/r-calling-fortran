INTEGER FUNCTION fib(n)
  integer n
  integer, parameter :: fib0 = 0, fib1 = 1
  integer back1, back2, i
  select case (n)
    case (:0);      fib = fib0
    case (1);       fib = fib1
    case default
      fib = fib1
      back1 = fib0
      do i = 2, n
        back2 = back1
        back1 = fib
        fib   = back1 + back2
      end do
  end select
END FUNCTION

SUBROUTINE fib_R_wrapper(n, answer)
  INTEGER n, answer, fib
  EXTERNAL fib
  answer = fib(n)
END SUBROUTINE