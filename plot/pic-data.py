import numpy as np
from PIL import Image

Nx = Ny = 100
image = Image.open('img/input.png')
image = image.resize((Nx, Ny))
image = image.convert('L')
data = np.asarray(image)

# threshold = 50
# matrix = (data < threshold).astype(dtype=np.uint8)
# image_test = Image.fromarray(matrix*255)
# image_test.save('image-test.jpg')

matrix = ((np.amax(data) - data) / np.amax(data))
# image_test = Image.fromarray(matrix*255)
# image_test.save('image-test.jpg')

with open('out/init_matrix.dat', 'w') as f:
    for i in range(Nx):
        for j in range(Ny):
            f.write(str(matrix[i, j]))
            f.write('\n')

print(matrix)
print(np.amax(data))
