def explicit_ftcs(dx, dt, total_time, alpha):
    N = int(1 / dx)
    x = np.linspace(0, 1, N + 1)
    T = initial_T(x)
    num_steps = int(total_time / dt)
    r = alpha * dt / dx**2
    if r > 0.5:
        print(f"Warning: r = {r} > 0.5, may be unstable")
    for _ in range(num_steps):
        T_new = T.copy()
        for i in range(1, N):
            T_new[i] = T[i] + r * (T[i + 1] - 2 * T[i] + T[i - 1])
        T = T_new

  def explicit_ftcs_unstable(dx, dt, total_time, alpha):
    N = int(1 / dx)
    x = np.linspace(0, 1, N + 1)
    T = initial_T(x)
    num_steps = int(total_time / dt)
    r = alpha * dt / dx**2
    T_center = []
    for _ in range(num_steps):
        T_new = T.copy()
        for i in range(1, N):
            T_new[i] = T[i] + r * (T[i + 1] - 2 * T[i] + T[i - 1])
        T = T_new
        T_center.append(T[N // 2])  # x = 0.5 m
    return T_center
