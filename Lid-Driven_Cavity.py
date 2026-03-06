import numpy as np
import matplotlib.pyplot as plt

def Ghia_data(Re):
    # Data set from Ghia
    if Re == 100:
        U = np.array([
        0,-0.03717,-0.04192,-0.04775,-0.06434,-0.10150,-0.15662,
        -0.21090,-0.20581,-0.13641,0.00332,0.23151,0.68717,
        0.73722,0.78871,0.84123,1.0000])
        y1 = np.array([0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,
                0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1])
        V = np.array([
        0,0.09233,0.10091,0.10890,0.12317,0.16077,0.17507,
        0.17527,0.05454,-0.24533,-0.22445,-0.16914,
        -0.10313,-0.08864,-0.07391,-0.05906,0])
        x1 = np.array([0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5,
                0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1])

    if Re == 400:
        U = np.array([
        0,-0.08186,-0.09266,-0.10338,-0.14612,-0.24299,-0.32726,
        -0.17119,-0.11777,0.02135,0.16256,0.29093,0.55892,
        0.61756,0.68439,0.75837,1.0000])
        y1 = np.array([0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,
                0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1])
        V = np.array([
        0,0.18360,0.19713,0.20920,0.22695,0.28124,0.30203,
        0.30174,0.05186,-0.38598,-0.44993,-0.38227,
        -0.22847,-0.12954,-0.15663,-0.12146,0])
        x1 = np.array([0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5,
                0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1])
        
    if Re == 1000:
        U = np.array([
        0,-0.18109,-0.20196,-0.22220,-0.29730,-0.38289,-0.27805,
        -0.10648,-0.06080,0.05702,0.18719,0.33304,0.46604,
        0.51117,0.57492,0.65928,1.0000])
        y1 = np.array([0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,
                0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1])
        V = np.array([
        0,0.27485,0.29012,0.30353,0.32627,0.37095,0.33075,
        0.32235,0.02526,-0.31966,-0.42665,-0.51550,
        -0.39188,-0.33714,-0.27669,-0.21388,0])
        x1 = np.array([0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5,
                0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1]) 
    if Re == 3200:
        U = np.array([
        0,-0.32407,-0.35344,-0.37827,-0.41933,-0.34323,-0.24427,
        -0.086636,-0.04272,0.07156,0.19791,0.34682,0.46101,
        0.46547,0.48296,0.53236,1.0000])
        y1 = np.array([0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,
                0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1])
        V = np.array([
        0,0.39560,0.40917,0.41906,0.42768,0.37119,0.29030,
        0.28188,0.00999,-0.31184,-0.37401,-0.44307,
        -0.54053,-0.52357,-0.47425,-0.39017,0])
        x1 = np.array([0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5,
                0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1]) 
             
    return x1, y1, U, V

def poisson_solver_for_psi(omega, psi):
    # Solving for psi (Poisson equation) using Jacobi Iteration
    psi_new = np.copy(psi)
    iteration_psi, convergence = 0, 0
    while True:
        iteration_psi += 1
        psi_temp = np.copy(psi_new)
        psi_new[1:-1, 1:-1] =(0.25 * (omega[1:-1, 1:-1] * dx**2 + psi_new[:-2, 1:-1] + 
                                      psi_new[2:, 1:-1] + psi_new[1:-1, :-2] + psi_new[1:-1, 2:]))
        error = np.abs(np.max(psi_new - psi_temp))
        if iteration_psi > 200:
            convergence = 1
            break
        if error < 1e-3:
            break     
            
    # Updating Boundary condition for psi
    psi_new[:,0] = psi_new[:,1]
    psi_new[:,Nx+1] = psi_new[:,Nx]
    psi_new[Nx+1,:] = psi_new[Nx,:]
    psi_new[0,:] = dy + psi_new[1,:]
    return psi_new, convergence

def explicit_solver_for_omega(u, v, omega, psi_new, dx, dy, dt, Re):
    omega_new = np.copy(omega)
    # Solving for omega using Explicit Method (Upwind scheme for convective terms and central difference for diffusive terms)
    u_pos = np.maximum(u[1:-1,1:-1], 0)
    u_neg = np.minimum(u[1:-1,1:-1], 0)
    convective_x = (u_pos * (omega[1:-1,1:-1] - omega[1:-1, :-2]) + 
                u_neg * (omega[1:-1, 2:]   - omega[1:-1,1:-1])) / dx
    v_pos = np.maximum(v[1:-1,1:-1], 0)
    v_neg = np.minimum(v[1:-1,1:-1], 0)
    convective_y = (v_pos * (omega[1:-1,1:-1] - omega[:-2, 1:-1]) + 
                v_neg * (omega[2:, 1:-1]   - omega[1:-1,1:-1])) / dy

    diffusive_x = (omega[1:-1, 2:] - 2 * omega[1:-1, 1:-1] + omega[1:-1, :-2]) / dx**2
    diffusive_y = (omega[:-2, 1:-1] - 2 * omega[1:-1, 1:-1] + omega[2:, 1:-1]) / dy**2
    omega_new[1:-1, 1:-1] = omega[1:-1, 1:-1] - dt * (convective_x + convective_y) +  \
                                                    (dt / Re) * (diffusive_x + diffusive_y)
                
    # Update Boundary Condition for omega
    omega_new[:,0] = - 16 * psi_new[:,0] / dx**2 - omega_new[:,1]
    omega_new[:,Nx+1] = - 16 * psi_new[:,Nx] / dx**2 - omega_new[:,Nx]
    omega_new[Nx+1,:] = - 16 * psi_new[Nx+1,:] / dy**2 - omega_new[Nx,:]
    omega_new[0,:] = - 8 / dy - 16 * psi_new[1,:] / dy**2 - omega_new[1,:]  
    return omega_new

def plotter(Re, t, iteration, convergence, psi, u, v, x1, y1, U, V):  
    print(f"Re={Re}, Time={t:.2f}, Iterations={iteration}, Convergence error={convergence: .2e}")
    x = np.linspace(0, Lx, Nx+2)
    y = np.linspace(0, Ly, Ny+2)
    X, Y = np.meshgrid(x, y,)

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
            
    # For Streamlines
    axs[0].contour(X, Y, psi, levels=30, colors='black', linestyles='solid', linewidths=0.5)
    axs[0].set_title(f'streamlines_Re_{Re}_t_{t:.2f}.png')
    axs[0].set_xlabel('X')
    axs[0].set_ylabel('Y')
    axs[0].invert_yaxis()

    # For Comparing U-velocity with Ghia
    axs[1].plot(y1, U, 'o', color='red', label='Ghia results', linewidth=2)
    axs[1].plot(y, u[::-1, (Ny + 1) // 2], '-g', label='Present code results', linewidth=2)
    axs[1].set_xlabel('Y')
    axs[1].set_ylabel('U velocity')
    axs[1].set_title(f'U_velocity_comparison_Re_{Re}_t_{t:.2f}')
    axs[1].legend()
    axs[1].grid(True)

    # For Comparing V-velocity with Ghia
    axs[2].plot(x1, V, 'o', color='red', label='Ghia results', linewidth=2)
    axs[2].plot(x, -v[(Nx + 1) // 2, :], '-g', label='Current code results', linewidth=2)
    axs[2].set_xlabel('X')
    axs[2].set_ylabel('V velocity')
    axs[2].set_title(f'V_velocity_comparison_Re_{Re}_t_{t:.2f}')
    axs[2].legend()
    axs[2].grid(True)
    plt.savefig(f'Lid_Driven_Cavity_Re_{Re}_t_{t:.2f}.png')

# Parameters
Re = [100, 400, 1000, 3200]       # Reynolds number
Nx = Ny = 129       # Number of grid points in x and y direction
Lx = Ly = 1      # Length of the cavity in x and y direction
dx = dy = Lx / (Nx-1)

for Re in Re:
    x1, y1, U, V = Ghia_data(Re)

    # Varible definition
    t, iteration = 0, 0
    omega = np.zeros((Ny+2, Nx+2))
    psi = np.zeros_like(omega)
    u = np.zeros_like(omega)
    v = np.zeros_like(omega)

    # Boundary Conditions at top wall
    u[0,:] = 1   
    omega[0, :] = - 8 / dy - 16 * psi[1,:] / dy**2 - omega[1,:]
    psi[0, :] = dy + psi[1,:]

    # Solver
    while True:
        # Adaptive time stepping based on CFL condition for stability
        convective_limit = 1.0 / ((np.max(np.abs(u))/dx) + (np.max(np.abs(v))/dy) + 1e-12)
        diffusive_limit  = Re * dx**2 / 4
        dt = 0.2 * min(convective_limit, diffusive_limit) 
        t += dt
        iteration += 1

        omega_new = explicit_solver_for_omega(u, v, omega, psi, dx, dy, dt, Re)
        convergence = np.abs(np.max(omega_new - omega))
        omega = np.copy(omega_new)

        psi_new, convergence_psi = poisson_solver_for_psi(omega_new, psi)
        if convergence_psi == 1:
            print("Jacobi didn't converging after", iteration, "iterations")
            break
        psi = np.copy(psi_new)

        # Update velocities
        u[1:-1, 1:-1] = (psi[:-2, 1:-1] - psi[2:, 1:-1]) / (2 * dy)
        v[1:-1, 1:-1] = -(psi[1:-1, :-2] - psi[1:-1, 2:]) / (2 * dx)

        print(f"Re={Re}, Time={t:.2f}, Iterations={iteration}, Convergence error={convergence: .2e}")

        if t > 500 or (convergence < 1e-6 and t > 1):
            plotter(Re, t, iteration, convergence, psi, u, v, x1, y1, U, V)
            break
