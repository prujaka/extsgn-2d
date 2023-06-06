import matplotlib.pylab as plt
import numpy as np


n_x = 800
n_y = 800

with open('out/res.dat') as f:
    lines = [line.strip() for line in f.readlines()]
    x = np.array([float(line.split()[0]) for line in lines])
    y = np.array([float(line.split()[1]) for line in lines])
    h = np.array([float(line.split()[2]) for line in lines])
    u = np.array([float(line.split()[3]) for line in lines])
    v = np.array([float(line.split()[4]) for line in lines])
    eta = np.array([float(line.split()[5]) for line in lines])
    w = np.array([float(line.split()[6]) for line in lines])

x = x.reshape((n_x, n_y))
y = y.reshape((n_x, n_y))
h = h.reshape((n_x, n_y))
u = u.reshape((n_x, n_y))
v = v.reshape((n_x, n_y))
eta = eta.reshape((n_x, n_y))
w = w.reshape((n_x, n_y))

k = n_y // 2 - 1

fig, axs = plt.subplots(3, 2, figsize=(10, 10))
axs[0, 0].plot(x[:, k], h[:, k], label='Depth h')
axs[0, 0].set_ylabel('m')

axs[0, 1].plot(x[:, k], u[:, k], label='Horizontal velocity u')
axs[0, 1].set_ylabel('m/s')

axs[1, 0].plot(x[:, k], v[:, k], label='Vertical velocity u')
axs[1, 0].set_ylabel('m/s')

axs[1, 1].plot(x[:, k], eta[:, k], label='eta')
axs[1, 1].set_ylabel('m')

axs[2, 0].plot(x[:, k], w[:, k], label='w')
axs[2, 0].set_ylabel('m/s')

for i in range(3):
    for j in range(2):
        if i == 2 and j != 0:
            break
        axs[i, j].set_xlabel('x (m)')
        axs[i, j].legend()
        axs[i, j].grid()

plt.savefig('huvetaw-1d.png', dpi=600)
plt.close()
