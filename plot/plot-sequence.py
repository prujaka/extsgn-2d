from xsgnplot import Solution, truncate_colormap
from os import listdir
from os.path import isfile, join
from matplotlib.colors import Normalize


if __name__ == '__main__':
    data_dir = 'out'
    data_files = [join(data_dir, f) for f in listdir(data_dir)
                  if isfile(join(data_dir, f)) and 'res_t=' in f]
    data_files = sorted(data_files)

    maxs = []
    mins = []
    for file in data_files:
        res = Solution(file)
        maxs.append(res.h.max())
        mins.append(res.h.min())

    h_max = max(maxs)
    h_min = min(mins)

    print(h_min, h_max)

    color_norm = Normalize(vmin=h_min, vmax=h_max)

    cmap = truncate_colormap('ocean_r', cutoff_percentage=0.8)

    for (i, file) in enumerate(data_files):
        res = Solution(file)
        img_file = file.replace('out', 'img').replace('.dat', '.png')
        res.plot_artsy(img_file, cmap=cmap, norm=color_norm)
        print(f'Generating {img_file}, {i + 1} out of {len(data_files)}')
