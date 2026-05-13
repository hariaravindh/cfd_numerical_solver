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
