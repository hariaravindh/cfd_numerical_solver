def initial_T(x):
    return np.where(x <= 0.5, 200 * x, 200 - 200 * x)
