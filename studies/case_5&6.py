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
