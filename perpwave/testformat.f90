       program testformat
       implicit none
       doubleprecision xx(3)
       xx(1) = 29304
       xx(2) = 12
       xx(3) = 4.8

       write(6,101) xx(1:3)

101    format(3p3e15.6)

       end