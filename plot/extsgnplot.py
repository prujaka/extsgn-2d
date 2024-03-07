import numpy as np
import matplotlib.pyplot as plt
import random

n_x = 100
n_y = 100


def create_testdata(file, nx=10, ny=5):
    with open(file, 'w') as f:
        for i in range(nx):
            for j in range(ny):
                rand = round(random.uniform(-0.2, 0.2), 2)
                s = ' '.join(map(str, [i, j, 1 + rand]))
                f.write(s + '\n')


class Solution:
    def __init__(self, file):
        with open(file) as f:
            lines = [line.strip() for line in f.readlines()]
            x = np.array([float(line.split()[0]) for line in lines])
            y = np.array([float(line.split()[1]) for line in lines])
            h = np.array([float(line.split()[2]) for line in lines])
            u = np.array([float(line.split()[3]) for line in lines])
            v = np.array([float(line.split()[4]) for line in lines])
            eta = np.array([float(line.split()[5]) for line in lines])
            w = np.array([float(line.split()[6]) for line in lines])

            self.x = x.reshape((n_x, n_y))
            self.y = y.reshape((n_x, n_y))
            self.h = h.reshape((n_x, n_y))
            self.u = u.reshape((n_x, n_y))
            self.v = v.reshape((n_x, n_y))
            self.eta = eta.reshape((n_x, n_y))
            self.w = w.reshape((n_x, n_y))
            self.t = float(lines[0].split()[7])

            grad_h = np.gradient(self.h)
            grad_h_abs = np.sqrt(grad_h[0]**2 + grad_h[1]**2)
            schlieren = np.log(1 + np.log(1 + 25 * grad_h_abs))
            self.schlieren = schlieren

    def plot_schlieren(self, file):
        fig, ax = plt.subplots()
        fig.set_size_inches(5, 5)

        ax.contourf(self.x, self.y, self.schlieren, levels=100,
                    cmap=plt.get_cmap('Blues'))

        plt.savefig(file, dpi=1200)
        plt.close()


if __name__ == '__main__':
    file = '../out/res.dat'
    png = '../img/schlieren-2d.png'
    solution = Solution(file)
    solution.plot_schlieren(png)
