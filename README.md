# Numerical-Wave-Equation-Solver

Uses the finite difference method to solve the wave equation.

# Executing

```
pip install -r requirements.txt
python wave.py
```

As of now, the code has fixed boundary conditions such that the string is fixed on both ends

# Notes

You can change the paramenters, dt,T, dx (indirectly), L and c.
```
c: the wave speed
L: the string length
T: the total time to simulate
dt: the difference between lines in the time grid
dx: difference between points in the space grid
```

Note: You must make sure the courant number (`c*dt/dx`) is less than or equal to 1. Otherwise the finite difference method
blows up. This book explains it in much detail.  [analysis of difference equations for wave equation](https://hplgit.github.io/fdm-book/doc/pub/book/sphinx/._book007.html#wave-pde1-num-dispersion)
