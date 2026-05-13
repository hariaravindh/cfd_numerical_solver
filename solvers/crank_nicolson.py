def crank_nicolson(dx, dt, total_time, alpha):
    N = int(1 / dx)
    x = np.linspace(0, 1, N + 1)
    T = initial_T(x)
    num_steps = int(total_time / dt)
    r = alpha * dt / dx**2
    # Tridiagonal matrix A
    A = np.diag((1 + r) * np.ones(N - 1)) + \
        np.diag(-r / 2 * np.ones(N - 2), k=-1) + \
        np.diag(-r / 2 * np.ones(N - 2), k=1)
    for _ in range(num_steps):
        b = (r / 2) * T[:-2] + (1 - r) * T[1:-1] + (r / 2) * T[2:]
        T_new = np.linalg.solve(A, b)
        T[1:N] = T_new
    return x, T
