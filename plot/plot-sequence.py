from xsgnplot import Solution
from os import listdir
from os.path import isfile, join

if __name__ == '__main__':
    data_dir = 'out'
    data_files = [join(data_dir, f) for f in listdir(data_dir)
                  if isfile(join(data_dir, f)) and 'res_t=' in f]
    data_files = sorted(data_files)
    for (i, file) in enumerate(data_files):
        res = Solution(file)
        img_file = file.replace('out', 'img').replace('.dat', '.png')
        res.plot_schlieren(img_file)
        print(f'Generating {file}, {i + 1} out of {len(data_files)}')
