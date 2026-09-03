 * Solving problem of Hestavan's book -> pg-27 == Ex-2.4
 * ----------------------------------------------------
   U_t + a U_x = 0, x in [0, 2*PI]
   a = -2 * PI

   with periodic boundary conditions and initial condition
   u(x, 0) = sin(l x), l = (2 * PI / lambda)
   lambda == wavelength
   nodes chosen ==> GLL

   Numerical scheme:
   Upwind flux-RK4(explicit)

   keep delta_t so small to ensure that timestep errors can be neglected

   K = number of elements [2, 4, 8, 16, 32, 64]
   h = (2 * PI) / K
   N = local order of approximation [1, 2, 4, 8]

Objectives:
    1.) Plot global "L2 error vs K elements" at final time T = PI
    2.) Plot "L2 error norm vs N (local order of approximation)"
    3.) Create a table of Global L2-errors when solving the wave equation using
K elements with each N(local order of approximation). Also check p+1 convergence
rate for a fixed N? 4.) Create a table of scaled execution times when solving
wave equation using K elements each with a different N(local order of
approximation) take time {N = 1, K = 2} as reference.

    5.) Table for Global L2 errors as a function of final time T when solving
wave eq. using K elements each with a local order of approximation, N.

    Final time (T) = PI, 10 PI, 100 PI, 1000 PI, 2000 PI;
    (N,K) = (2,4)
    (N,K) = (4,2)
    (N,K) = (4,4)

General advice -> Keep it raw.


Some key points:
1.) GLL points are chosen for nodal values of solution.
2.) Since same GLL points can be thus used for quadrature rule.
3.) Let number of quadrature points inside the finite element be Nq.
4.) GLL quadrature yeilds exact integration of a polynomial upto degree (2*Nq - 3), If we take Nq points at GLL nodes.
5.) Within an element for a polynomial of degree N, we require N+1 nodal values. (irrespective of where it is situated)
6.) If these N+1 nodes where solution is defined, are also taken as quadrature nodal locations, then GLL quadrature would yeild quadrature that is exact for polynomial degree upto 2*(N+1)-3 = 2N-1.
7.) this is a problem in evaluating Mass matrix, which has an integrand of l_j * l_i.
max_degree(l_j * l_i) is 2*N.

Cases:
1.) This leads to Mass matrix being not evaluated exactly, if quadrature points and solution points are the same. If quadrature was taken to be inexact using N+1 points, mass matrix would be diagonal.
2.) If GL points were taken instead of GLL, then we need to find nodal values using interpolation. Hence, requiring interpolation matrix. Exact integration of Mass matrix, but full and illconditioned
