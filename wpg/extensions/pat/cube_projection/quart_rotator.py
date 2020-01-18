import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


from pyquaternion import Quaternion

nx = 250

cube_nx = 30

sample = np.zeros((nx, nx, nx))

cube_s = int((nx - cube_nx) / 2)
cube_f = int((nx + cube_nx) / 2)

sample[cube_s:cube_f, cube_s:cube_f, cube_s:cube_f] = 1

s0i = np.argwhere(sample == 1)

sample = np.zeros((nx, nx, nx))

q1 = Quaternion(axis=[0, 1, 0], angle=-3.14159265 / 4)

for i in range(cube_nx ** 3):
    s0i[i, :] = np.round(q1.rotate(s0i[i, :]))

min_x = np.min(s0i[:, 0])
min_y = np.min(s0i[:, 1])
min_z = np.min(s0i[:, 2])

s0i[:, 0] -= min_x - 1
s0i[:, 1] -= min_y - 1
s0i[:, 2] -= min_z - 1


sample[s0i[:, 0], s0i[:, 1], s0i[:, 2]] = 1

generation_ims = []
fig = plt.figure()
for i in range(nx):
    generation_im = plt.imshow(sample[:, :, i], animated=True)

    generation_ims.append([generation_im])

ani = animation.ArtistAnimation(fig, generation_ims, interval=100, blit=False, repeat_delay=0)

plt.show()
