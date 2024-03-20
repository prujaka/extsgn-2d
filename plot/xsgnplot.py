import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import plotly.graph_objects as go
import random


def get_meshsize(parameters_path='src/m_parameters.f90'):
    with open(parameters_path, 'r') as f:
        lines = f.readlines()
    nx_line, ny_line = [line for line in lines
                        if 'NX' in line or 'NY' in line][:2]
    nx = int(nx_line.split()[-1])
    ny = int(ny_line.split()[-1])
    return nx, ny


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

        nx, ny = get_meshsize()

        self.x = x.reshape((nx, ny))
        self.y = y.reshape((nx, ny))
        self.h = h.reshape((nx, ny))
        self.u = u.reshape((nx, ny))
        self.v = v.reshape((nx, ny))
        self.eta = eta.reshape((nx, ny))
        self.w = w.reshape((nx, ny))
        self.t = t[0]

        mid_index = ny // 2 - 1

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

    def plot_schlieren(self, file, cmap='Blues'):
        fig, ax = plt.subplots()
        fig.set_size_inches(5, 5)

        ax.contourf(self.x, self.y, self.schlieren, levels=100,
                    cmap=plt.get_cmap(cmap))

        plt.savefig(file, dpi=300)
        plt.close()

    def plot_artsy(self, file, cmap='bone_r', norm=None):
        fig, ax = plt.subplots()
        fig.set_size_inches(5, 5)

        ax.contourf(self.x, -self.y, self.h.T, levels=100,
                    cmap=plt.get_cmap(cmap), norm=norm)

        ax.axis('off')
        # fig.patch.set_alpha(0)
        plt.savefig(file, dpi=300)
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

        plt.savefig(file, dpi=300)
        plt.close()

    def plot_3d_surface(self, file):
        fig = go.Figure(data=[go.Surface(x=-self.x, y=-self.y, z=self.h,
                                         colorscale='blues')])

        camera = dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=0, y=0, z=1.5),
        )

        fig.update_layout(
            scene=dict(zaxis=dict(nticks=1, range=[0, 4], title='')),
            width=700, height=700,
            scene_camera=camera,
            margin=dict(r=10, l=10, b=10, t=10)
        )

        fig.show()
        fig.write_image(file)


if __name__ == '__main__':
    file = 'out/res.dat'
    png_1d_sections = 'img/plot_1d_sections.png'
    png_2d_schlieren = 'img/plot_2d_schlieren.png'
    png_2d_artsy = 'img/plot_2d_artsy.png'
    png_3d_surface = 'img/plot_3d_surface.png'
    solution = Solution(file)

    cmap = truncate_colormap('ocean_r', cutoff_percentage=0.8)
    solution.plot_sections(png_1d_sections)
    solution.plot_artsy(png_2d_artsy, cmap=cmap)
    solution.plot_schlieren(png_2d_schlieren)
    solution.plot_3d_surface(png_3d_surface)
