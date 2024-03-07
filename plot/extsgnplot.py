import numpy as np
import matplotlib.pyplot as plt

n_x = 100
n_y = 100


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

            grad_h = np.gradient(h)
            grad_h_abs = np.sqrt(grad_h[0]**2 + grad_h[1]**2)
            schlieren = np.log(1 + np.log(1 + 25 * grad_h_abs))

            self.x = x.reshape((n_x, n_y))
            self.y = y.reshape((n_x, n_y))
            self.h = h.reshape((n_x, n_y))
            self.u = u.reshape((n_x, n_y))
            self.v = v.reshape((n_x, n_y))
            self.eta = eta.reshape((n_x, n_y))
            self.w = w.reshape((n_x, n_y))
            self.schlieren = schlieren
            self.t = float(lines[0].split()[7])


if __name__ == '__main__':
    file = '../out/res.dat'
    solution = Solution(file)
    print(solution.t)
