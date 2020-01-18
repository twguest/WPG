import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import scipy.ndimage.morphology as scipym

from pyquaternion import Quaternion

cube_nx = 50

cube = np.ones((cube_nx, cube_nx, cube_nx))

cube_points = np.argwhere(cube == 1)


# Problem angles
# q1 = Quaternion(axis=[1, 1, 1], angle=-3.14159265 / 4)
# q2 = Quaternion(axis=[0, 0, 1], angle=-3.14159265 / 4)
# q3 = q1*q2

# q1 = Quaternion(axis=[1, 0, 1], angle=-3.14159265 / 4)
# q2 = Quaternion(axis=[0, 0, 1], angle=-3.14159265 / 4)
# q3 = q1*q2


# q1 = Quaternion(axis=[1, 2, 1], angle=-3.14159265 / 3)
# q2 = Quaternion(axis=[0, 0, 1], angle=-3.14159265 / 6)
# q3 = q1*q2


# Good angles
q1 = Quaternion(axis=[2, 1, 0], angle=-3.14159265 / 7)
q2 = Quaternion(axis=[0, 1, 5], angle=-3.14159265 / 6)
q3 = q1*q2

q1 = Quaternion(axis=[3, 1, 0], angle=-3.14159265 / 7)
q2 = Quaternion(axis=[0, 1, 0], angle=-3.14159265 / 6)
q3 = q1*q2

q1 = Quaternion(axis=[0, 1, 1], angle=-3.14159265 / 7)
q2 = Quaternion(axis=[1, 1, 1], angle=-3.14159265 / 6)
q3 = q1*q2

q1 = Quaternion(axis=[1, 0, 0], angle=0)
q2 = Quaternion(axis=[1, 1, 1], angle=-2*3.14159265 / 3)
q3 = q1*q2

q1 = Quaternion(axis=[1, 0, 0], angle=0)
q2 = Quaternion(axis=[1, 1, 1], angle=-3.14159265 / 3)
q3 = q1*q2

rotated_cube = np.dot(cube_points, q3.rotation_matrix)

rotated_cube = np.round(rotated_cube - np.min(rotated_cube)).astype(int)

sample_nx = np.max(rotated_cube)+2

sample = np.zeros((sample_nx, sample_nx, sample_nx))

for i in range(cube_nx**3):
    sample[rotated_cube[i,0], rotated_cube[i,1],rotated_cube[i,2]] = 1




generation_ims = []
fig = plt.figure()
for i in range(np.max(rotated_cube)+2):

    sample_im = sample[:,:,i]

    #scipym.binary_fill_holes(sample_im,structure=np.ones((2,2)),output=sample_im)
    scipym.binary_dilation(sample_im, structure=np.ones((2, 2)), output=sample_im)

    generation_im = plt.imshow(sample_im, animated=True)

    generation_ims.append([generation_im])

ani = animation.ArtistAnimation(fig, generation_ims, interval=100, blit=False, repeat_delay=0)

plt.show()
