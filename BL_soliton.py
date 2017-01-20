
# Solving Benney-Luke equations for an incident solitary wave interacting with an oblique wall.
# Floriane Gidel, 2016.

from firedrake import *
import numpy as np
import time
op2.init()
parameters["coffee"]["O2"] = False

# Timing
t00 = time.time()

""" ________________ Parameters ________________ """
Lw = 500.0          # Length of the wall
Lc = 5.0            # Length of the incident channel
psi = pi/6          # Angle of incidence
Ld = Lw*sin(psi)    # Width of the incident channel
ep = 0.19           # Small amplitude parameter
mu = 0.02           # Small dispersion parameter
t=0.0               # Initial time
t_save = t          # Saving time
dt = 0.0028         # Time step
Tend = 150.0        # Final time
psi_inc = 0.0


""" ___________________ Mesh ___________________ """
mesh = Mesh("horizontal.msh")   # Load the mesh file
coords = mesh.coordinates       # access to coordinates

""" ______________ Function Space ______________ """
V = FunctionSpace(mesh, "CG", 2)    # Vector space

""" ___________ Define the functions ___________ """
eta_n0 = Function(V)    # eta(n)
phi_n0 = Function(V)    # phi(n)
q_n0 = Function(V)      # q(n)
eta_n1 = Function(V)    # eta(n+1)
phi_n1 = Function(V)    # phi(n+1)
q_n1 = Function(V)      # q(n+1)
eta_half = Function(V)  # eta(n+1/2)
phi_half = Function(V)  # phi(n+1/2)
q_half = Function(V)    # q(n+1/2)
inf_y=Function(V).interpolate(Expression("0.5*(1+copysign(1.0,x[1]-96.0))")) # used to measure the incident amplitude

# Unknown
eta = TrialFunction(V)
phi = TrialFunction(V)
q = TrialFunction(V)

# Test function
v = TestFunction(V)



""" _____________ Initial solution _____________ """
# Soliton's parameters
A = 1.0                             # Amplitude
C=0.5*(A + tan(psi)*tan(psi)/ep)    # Constant C (cf def of soliton in the paper)
dist = 0.5                          # distance (%) from the incident channel boundary
alpha = pi/2-psi
x0 = 0.5*Lc                         # Point (x0,y0) where to place the soliton
y0 = 0.0

# Expression of eta and phi
expr_eta = Expression("A*pow(cosh(sqrt(3*ep*A/(4*mu))*((x[0]-x0)  + (x[1]-y0)*tan(psi) -(t-t0)*(1+C*ep))),-2)",  A=A, x0=x0, y0=y0, psi=psi_inc,ep=ep, mu=mu, C=C, t=t, t0=0.0 )
expr_phi = Expression("A*sqrt(4*mu/(3*ep*A))*(tanh(sqrt(3*ep*A/(4*mu))*((x[0]-x0) + (x[1]-y0)*tan(psi)-(t-t0)*(1+C*ep)))+1) ", A=A, x0=x0, y0=y0, psi=psi_inc, ep=ep, mu=mu, C=C, t=t, t0=0.0)

## Initialization
phi_n0.interpolate(expr_phi)
eta_n0.interpolate(expr_eta)



""" ________________ Saving data ________________ """
phi_file = File('data/phi.pvd')                # potential phi numerical solution
eta_file = File('data/eta.pvd')                # surface deviation eta numerical solution
Ampl_file = open('data/amplitudes.txt', 'w')   # Incident and stem waves' amplitudes


""" _____________ Weak formulations _____________ """
# Get phi(n+1/2)
F_phi_half = (v*(phi_half - phi_n0)/(dt/2.0) + 0.5*mu*inner(grad(v),grad((phi_half-phi_n0)/(dt/2.0))) + v*eta_n0 + 0.5*ep*v*inner(grad(phi_half),grad(phi_half)))*dx
phi_problem_half = NonlinearVariationalProblem(F_phi_half,phi_half)
phi_solver_half = NonlinearVariationalSolver(phi_problem_half)

# Get q(n+1/2)
a_q_half = v*q*dx
L_q_half = 2.0/3.0*inner(grad(v),grad(phi_half))*dx
q_problem_half = LinearVariationalProblem(a_q_half,L_q_half,q_half)
q_solver_half  = LinearVariationalSolver(q_problem_half)

# Get eta(n+1)
a_eta = (v*eta/dt + 0.5*mu*inner(grad(v),grad(eta/dt)) - 0.5*inner(grad(v),grad(phi_half))*ep*eta)*dx
L_eta = (v*eta_n0/dt + 0.5*mu*inner(grad(v),grad(eta_n0)/dt) +  mu*inner(grad(v),grad(q_half)) + 0.5*inner(grad(v),grad(phi_half))*(2+ep*eta_n0))*dx
eta_problem = LinearVariationalProblem(a_eta,L_eta, eta_n1)
eta_solver  = LinearVariationalSolver(eta_problem)

# Get phi(n+1)
a_phi_n1 = (v*phi/(dt/2) + 0.5*mu*inner(grad(v),grad(phi/(dt/2))))*dx
L_phi_n1 = (v*phi_half/(dt/2) + 0.5*mu*inner(grad(v),grad(phi_half)/(dt/2)) - v*eta_n1 - 0.5*ep*v*inner(grad(phi_half),grad(phi_half)))*dx
phi_problem_n1 = LinearVariationalProblem(a_phi_n1,L_phi_n1, phi_n1)
phi_solver_n1  = LinearVariationalSolver(phi_problem_n1)


""" _____________ Time loop _____________ """
while(t < Tend):
    
    # --- Saving data --- #
    if t >= t_save:
        eta_max = max(eta_n0.dat.data)                  # Get the stem wave's amplitude
        init_max = max(inf_y.dat.data*eta_n0.dat.data)  # Get the incident wave's amplitude
        Ampl_file.write('%-10s %-10s %-10s %-10s %-10s\n' % (t,'', eta_max,'' ,init_max))
        eta_file << eta_n0                              # Save the surface deviation solution
        phi_file << phi_n0                              # Save the potential flow solution
        t_save = t_save + 10*dt                         # Update the saving time
        print t/Tend

    # --- Update time --- #
    t += dt

    # --- Solve the weak formulations --- #
    phi_solver_half.solve()
    q_solver_half.solve()
    eta_solver.solve()
    phi_solver_n1.solve()
    
    # --- Update the solutions --- #
    phi_n0.assign(phi_n1)
    eta_n0.assign(eta_n1)

Ampl_file.close()           # Close data file
print time.time() - t00     # Print computational time (s)
