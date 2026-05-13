def implicit_laasonen(dx, dt, total_time, alpha):
    N = int(1 / dx)
    x = np.linspace(0, 1, N + 1)
    T = initial_T(x)
    num_steps = int(total_time / dt)
    r = alpha * dt / dx**2
    # Tridiagonal matrix A
    A = np.diag((1 + 2 * r) * np.ones(N - 1)) + \
        np.diag(-r * np.ones(N - 2), k=-1) + \
        np.diag(-r * np.ones(N - 2), k=1)
    for _ in range(num_steps):
        b = T[1:N]
        T_new = np.linalg.solve(A, b)
        T[1:N] = T_new  # T[0] and T[N] remain 0
    return x, T
