from xsgnplot import Solution, truncate_colormap
from os import listdir
from os.path import isfile, join
import numpy as np


if __name__ == '__main__':
    data_dir = 'out'
    data_files = [join(data_dir, f) for f in listdir(data_dir)
                  if isfile(join(data_dir, f)) and 'res_t=' in f]
    data_files = sorted(data_files)

    cmap = truncate_colormap('ocean_r', cutoff_percentage=0.7)

    matrix_file = 'out/init_matrix.dat'
    matrix = np.loadtxt(matrix_file)

    init_data_max = matrix.max()
    init_data_min = matrix.min()
    print(init_data_max, init_data_min)

    for (i, file) in enumerate(data_files):
        res = Solution(file)
        img_file = file.replace('out', 'img').replace('.dat', '.png')
        res.plot_artsy(img_file, cmap=cmap)
        print(f'Generating {img_file}, {i + 1} out of {len(data_files)}')
