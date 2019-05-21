using Dierckx
x = [1,3,5,7]
y = [0,2,4, 6]

spl = Spline1D(x, y; k=2, bc="nearest")


