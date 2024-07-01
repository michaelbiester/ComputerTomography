# testRadonUtils.py

import numpy as np
import matplotlib.pyplot as plt
from RadonUtils import intersections, digitalProjectionLine

fname = fname1 = './images/radon_tf_img1.npy'
img1 = np.load(fname)
Ny, Nx = img1.shape
print(f"nr. rows (Ny): {Ny}; nr. columns (Nx): {Nx}")

# define physical dimensions
x_l = -2
x_u = 3
y_l = -1
y_u = 3

# define parameters of projection line
t_val = 1.0
theta_deg = 30.0

# compute intersections
intersectionCount, intersectionPoints = intersections(t_val, theta_deg, x_l, x_u, y_l, y_u)
point_x1y1, point_x2y2 = intersectionPoints
print(f"point (x1,y1): {point_x1y1}; point (x2,y2): {point_x2y2}")

# load image
img1 = np.load(fname1)
Ny, Nx = img1.shape
xIndices, yIndices, points = digitalProjectionLine(intersectionPoints, Nx, Ny, x_l, x_u, y_l, y_u)
print(f"points: {points}")
n_x1, n_y1, n_x2, n_y2 = points

fig1 = plt.figure(1, figsize=[8, 8])
ax_f1 = fig1.add_subplot(1, 1, 1)
# plot of image
ax_f1.imshow(img1, cmap='Greys' )
ax_f1.plot([n_x1, n_x2], [Ny - 1 - n_y1, Ny - 1 - n_y2], color='b')
ax_f1.plot(xIndices, Ny - 1 - yIndices, color='g')
ax_f1.grid(True)
ax_f1.set_ylabel('y')
ax_f1.set_title('test image1')

plt.show()