import numpy as np
from PIL import Image
from matplotlib import pyplot as plt
from xsgnplot import get_meshsize


nx, ny = get_meshsize()
image = Image.open('img/input.jpg')
image = image.resize((nx, ny))
image = image.convert('L')
data = np.asarray(image)

# threshold = 50
# matrix = (data < threshold).astype(dtype=np.uint8)
# image_test = Image.fromarray(matrix*255)
# image_test.save('image-test.jpg')

matrix = ((np.amax(data) - data) / np.amax(data)) + 1
plt.imshow(matrix, interpolation='nearest')
plt.savefig('img/image-test.png')
# image_test = Image.fromarray(matrix*255)
# image_test.save('img/image-test.png')

with open('out/init_matrix.dat', 'w') as f:
    for i in range(nx):
        for j in range(ny):
            f.write(str(matrix[i, j]))
            f.write('\n')

print(matrix)
print(np.amax(data))
