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
