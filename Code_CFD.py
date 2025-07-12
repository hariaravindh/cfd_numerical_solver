import numpy as np
import matplotlib.pyplot as plt

alpha = 0.02
def initial_T(x):
    return np.where(x <= 0.5, 200 * x, 200 - 200 * x)
def exact_T(x, t, alpha, m_max=100):
    sum_series = 0
    for m in range(1, m_max + 1):
        term = (1 / m**2) * np.sin(m * np.pi / 2) * np.sin(m * np.pi * x) * np.exp(-m**2 * np.pi**2 * alpha * t)
        sum_series += term
    return (800 / np.pi**2) * sum_series
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
    return x, T
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

dx = 0.1 #1
dt = 0.1
total_time = 10
x, T_exp = explicit_ftcs(dx, dt, total_time, alpha)
_, T_laa = implicit_laasonen(dx, dt, total_time, alpha)
_, T_cn = crank_nicolson(dx, dt, total_time, alpha)
T_exact = exact_T(x, total_time, alpha)

plt.figure()
plt.plot(x, T_exact, 'k-', label='Exact')
plt.plot(x, T_exp, 'r--', label='Explicit')
plt.plot(x, T_laa, 'g-.', label='Laasonen')
plt.plot(x, T_cn, 'b:', label='Crank-Nicolson')
plt.xlabel('x (m)')
plt.ylabel('T (K)')
plt.title('t = 10 h, dx = 0.1 m, dt = 0.1 h, r = 0.2')
plt.legend()
plt.show()

dt = 0.25 #2
x, T_exp = explicit_ftcs(dx, dt, total_time, alpha)
_, T_laa = implicit_laasonen(dx, dt, total_time, alpha)
_, T_cn = crank_nicolson(dx, dt, total_time, alpha)
T_exact = exact_T(x, total_time, alpha)

plt.figure()
plt.plot(x, T_exact, 'k-', label='Exact')
plt.plot(x, T_exp, 'r--', label='Explicit')
plt.plot(x, T_laa, 'g-.', label='Laasonen')
plt.plot(x, T_cn, 'b:', label='Crank-Nicolson')
plt.xlabel('x (m)')
plt.ylabel('T (K)')
plt.title('t = 10 h, dx = 0.1 m, dt = 0.25 h, r = 0.5')
plt.legend()
plt.show()

dx = 0.05 #3
dt = 0.0625
x, T_exp = explicit_ftcs(dx, dt, total_time, alpha)
_, T_laa = implicit_laasonen(dx, dt, total_time, alpha)
_, T_cn = crank_nicolson(dx, dt, total_time, alpha)
T_exact = exact_T(x, total_time, alpha)

plt.figure()
plt.plot(x, T_exact, 'k-', label='Exact')
plt.plot(x, T_exp, 'r--', label='Explicit')
plt.plot(x, T_laa, 'g-.', label='Laasonen')
plt.plot(x, T_cn, 'b:', label='Crank-Nicolson')
plt.xlabel('x (m)')
plt.ylabel('T (K)')
plt.title('t = 10 h, dx = 0.05 m, dt = 0.0625 h, r = 0.5')
plt.legend()
plt.show()

dx = 0.1 #4
dt = 0.3
total_time = 20
T_center = explicit_ftcs_unstable(dx, dt, total_time, alpha)
num_steps = int(total_time / dt)
time = dt * np.arange(1, num_steps + 1)

plt.figure()
plt.plot(time, T_center, 'r-')
plt.xlabel('Time (h)')
plt.ylabel('T at x = 0.5 m (K)')
plt.title('Centerline T, Explicit, dx = 0.1 m, dt = 0.3 h, r = 0.6')
plt.show()

dx = 0.1 #5&6
total_time = 10
for r in [1, 2, 3, 4]:  # Test until agreement worsens
    dt = 0.5 * r
    x, T_laa = implicit_laasonen(dx, dt, total_time, alpha)
    _, T_cn = crank_nicolson(dx, dt, total_time, alpha)
    T_exact = exact_T(x, total_time, alpha)
    
    plt.figure()
    plt.plot(x, T_exact, 'k-', label='Exact')
    plt.plot(x, T_laa, 'g-.', label='Laasonen')
    plt.plot(x, T_cn, 'b:', label='Crank-Nicolson')
    plt.xlabel('x (m)')
    plt.ylabel('T (K)')
    plt.title(f't = 10 h, dx = 0.1 m, dt = {dt} h, r = {r}')
    plt.legend()
    plt.show()
