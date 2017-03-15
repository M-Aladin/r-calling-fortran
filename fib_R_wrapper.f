SUBROUTINE fib_R_wrapper(n, answer)
  INTEGER n, answer, fib
  EXTERNAL fib
  answer = fib(n)
END SUBROUTINE