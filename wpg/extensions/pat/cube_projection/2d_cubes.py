import numpy as np
import matplotlib.pyplot as plt
import random as rn
import imageio
import scipy.misc


nx = 1000  # number of pixels
dx = 1  # size of the detector side (um)
px = dx / nx  # size of pixel


cube_dx = 0.1 # size of a cube (um)

cube_nx = int(cube_dx * nx/dx) # number of pixels for a cube

def make_rot_cube(cube2_length=20):


    # side length of a matrix that contains the cube
    cube2_corner_length = int(np.round(np.cos(np.pi / 4) * cube2_length))

    # one corner of the matrix that contains the cube
    cube2_corner = np.zeros((cube2_corner_length, cube2_corner_length))

    min_ij = cube2_corner.size
    #loop through the matrix and only put values in the bottom right corner, increasing closer to the corner
    for i in range(cube2_corner_length):
        for j in range(cube2_corner_length):
            if i > cube2_corner_length - j:
                cube2_corner[i, j] = i+j

                if min_ij > i+j:
                    min_ij = i+j
                    print(min_ij)


    cube2_corner =  np.where(cube2_corner>min_ij, cube2_corner - min_ij, 0)


    #make a matrix that will contain the whole cube
    cube2 = np.zeros((2 * cube2_corner_length, 2 * cube2_corner_length))

    #top left corner of the cube matrix will be the coner we first calculated
    cube2[0:cube2_corner_length, 0:cube2_corner_length] = cube2_corner

    #for every other corner, rotate the matrix we made
    cube2[cube2_corner_length:2 * cube2_corner_length, 0:cube2_corner_length] = np.rot90(cube2_corner)
    cube2[cube2_corner_length:2 * cube2_corner_length, cube2_corner_length:2 * cube2_corner_length] = np.rot90(cube2_corner,                                                                                                    2)
    cube2[0:cube2_corner_length, cube2_corner_length:2 * cube2_corner_length] = np.rot90(cube2_corner, 3)

    #Scale the values in the cube so the longest distance (corner to corner) matches that of a cube of length 1 (root(3))
    cube2 = (cube2 / np.max(cube2)) * np.sqrt(3)

    return cube2, cube2_corner

def make_cube(cube_length=20):
    cube = np.ones((cube_length,cube_length))
    return cube

def add_cube_to_sample(cube, sample):

    new_sample = np.copy(sample)

    pos = (rn.randint(0, sample.shape[0]), rn.randint(0,sample.shape[1]))
    count = 0
    while pos[0] + cube.shape[0] > sample.shape[0] or pos[1] + cube.shape[1] > sample.shape[1]:
        print('Finding new position')
        count +=1
        pos = (rn.randint(0, sample.shape[0]), rn.randint(0, sample.shape[1]))
        if count>20:
            print("Can't find new position, try a different cube size")
            exit(1)





    new_sample[pos[0]:pos[0]+cube.shape[0],pos[1]:pos[1]+cube.shape[1]] = sample[pos[0]:pos[0]+cube.shape[0],pos[1]:pos[1]+cube.shape[1]]+ cube

    sample = np.where(new_sample==0, sample, new_sample)

    return sample










sample_plane = np.zeros((nx, nx))

cube2, corner = make_rot_cube(cube_nx)
cube = make_cube(cube_nx)



sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)

sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)

sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)

sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)

sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)

sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)
sample_plane = add_cube_to_sample(cube, sample_plane)
sample_plane = add_cube_to_sample(cube2, sample_plane)



scipy.misc.imsave('cube_samples.tiff', sample_plane)


plt.imshow(sample_plane,cmap='plasma')

plt.show()



