import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.colors import LinearSegmentedColormap

n_x = 500
n_y = 500


def truncate_colormap(cmap_name, cutoff_percentage=0.8):
    original_cmap = plt.get_cmap(cmap_name)
    num_colors = original_cmap.N
    cutoff_index = int(num_colors * cutoff_percentage)
    new_colors = original_cmap(np.linspace(0, cutoff_index / num_colors,
                                           cutoff_index))
    truncated_cmap_name = f'truncated_{cmap_name}'
    truncated_cmap = LinearSegmentedColormap.from_list(truncated_cmap_name,
                                                       new_colors)

    return truncated_cmap


def create_testdata(file, nx=10, ny=5):
    with open(file, 'w') as f:
        for i in range(nx):
            for j in range(ny):
                rand = round(random.uniform(-0.2, 0.2), 2)
                s = ' '.join(map(str, [i, j, 1 + rand]))
                f.write(s + '\n')


class Solution:
    def __init__(self, file):
        x, y, h, u, v, eta, w, t = np.loadtxt(file, unpack=True)

        self.x = x.reshape((n_x, n_y))
        self.y = y.reshape((n_x, n_y))
        self.h = h.reshape((n_x, n_y))
        self.u = u.reshape((n_x, n_y))
        self.v = v.reshape((n_x, n_y))
        self.eta = eta.reshape((n_x, n_y))
        self.w = w.reshape((n_x, n_y))
        self.t = t[0]

        mid_index = n_y // 2 - 1

        self.section_x = self.x[:, mid_index]
        self.section_h = self.h[:, mid_index]
        self.section_u = self.u[:, mid_index]
        self.section_v = self.v[:, mid_index]
        self.section_eta = self.eta[:, mid_index]
        self.section_w = self.w[:, mid_index]

        grad_h = np.gradient(self.h)
        grad_h_abs = np.sqrt(grad_h[0] ** 2 + grad_h[1] ** 2)
        schlieren = np.log(1 + np.log(1 + 25 * grad_h_abs))
        self.grad_h_abs = grad_h_abs
        self.schlieren = schlieren

    def plot_schlieren(self, file):
        fig, ax = plt.subplots()
        fig.set_size_inches(5, 5)

        ax.contourf(self.x, self.y, self.schlieren, levels=100,
                    cmap=plt.get_cmap('Blues'))

        plt.savefig(file, dpi=1200)
        plt.close()

    def plot_artsy(self, file, cmap='bone_r'):
        fig, ax = plt.subplots()
        fig.set_size_inches(5, 5)

        ax.contourf(self.x, -self.y, self.h.T, levels=100,
                    cmap=plt.get_cmap(cmap))

        ax.axis('off')
        # fig.patch.set_alpha(0)
        plt.savefig(file, dpi=1200)
        plt.close()

    def plot_sections(self, file):
        fig, axs = plt.subplots(3, 2, figsize=(10, 10))

        x = self.section_x
        h = self.section_h
        u = self.section_u
        v = self.section_v
        eta = self.section_eta
        w = self.section_w

        axs[0, 0].plot(x, h, label='Depth h')
        axs[0, 0].set_ylabel('m')

        axs[0, 1].plot(x, u, label='Horizontal velocity u')
        axs[0, 1].set_ylabel('m/s')

        axs[1, 0].plot(x, v, label='Vertical velocity u')
        axs[1, 0].set_ylabel('m/s')

        axs[1, 1].plot(x, eta, label='eta')
        axs[1, 1].set_ylabel('m')

        axs[2, 0].plot(x, w, label='w')
        axs[2, 0].set_ylabel('m/s')

        for i in range(3):
            for j in range(2):
                if i == 2 and j != 0:
                    break
                axs[i, j].set_xlabel('x (m)')
                axs[i, j].legend()
                axs[i, j].grid()

        plt.savefig(file, dpi=600)
        plt.close()


if __name__ == '__main__':
    file = 'out/res.dat'
    png = 'img/schlieren-2d.png'
    png_artsy = 'img/artsy.png'
    png_1d = 'img/huvetaw-1d.png'
    solution = Solution(file)

    cmap = truncate_colormap('ocean_r', cutoff_percentage=0.7)
    solution.plot_artsy(png_artsy, cmap=cmap)
    # solution.plot_schlieren(png)
    # solution.plot_sections(png_1d)
