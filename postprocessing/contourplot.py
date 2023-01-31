import random
import matplotlib.pylab as plt
import numpy as np


def create_testdata(file, nx=10, ny=5):
    with open(file, 'w') as f:
        for i in range(nx):
            for j in range(ny):
                rand = round(random.uniform(-0.2, 0.2), 2)
                s = ' '.join(map(str, [i, j, 1 + rand]))
                f.write(s + '\n')


# create_testdata('data.dat')

n_x = 500
n_y = 500
with open('res.dat') as f:
    lines = [line.strip() for line in f.readlines()]
    x = np.array([float(line.split()[0]) for line in lines])
    y = np.array([float(line.split()[1]) for line in lines])
    z = np.array([float(line.split()[2]) for line in lines])

x = x.reshape((n_x, n_y))
y = y.reshape((n_x, n_y))
z = z.reshape((n_x, n_y))
grad_z = np.gradient(z)
abs_grad_z = np.sqrt(grad_z[0]**2 + grad_z[1]**2)

fig, ax = plt.subplots()
fig.set_size_inches(5, 5)

ax.contourf(x, y, np.log(1 + np.log(1 + 25*abs_grad_z)), levels=100,
            cmap=plt.get_cmap('Blues'))

plt.savefig('schlieren-2d.png', dpi=1200)
plt.show()
