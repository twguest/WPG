import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import scipy.ndimage.morphology as scipym
import random as rn
import imageio

from pyquaternion import Quaternion


def make_cube_stamp(cube_nx, axis=[1, 1, 1], angle=np.pi / 6, save=False):
    cube = np.ones((cube_nx, cube_nx, cube_nx))  # init a matrix that defines the cube

    qs = Quaternion(
        axis=axis, angle=angle
    )  # make an object that defines the rotation of the points

    # Get the coordinates in the cube
    # cube_points = [ [0,0,0], [0,0,1], ..., [nx, nx, nx] ]
    cube_points = np.argwhere(cube == 1)

    rotated_cube = np.dot(
        cube_points, qs.rotation_matrix
    )  # rotate the points to find the rotated cube

    rotated_cube = np.round(rotated_cube - np.min(rotated_cube)).astype(
        int
    )  # translate the points so everything is +

    for i in range(3):
        rotated_cube[:, i] -= np.min(rotated_cube[:, i]) - 0

    # turn all the points to integers
    rotated_cube = rotated_cube.astype(int)

    # taking the max value of the rotated cube will give the size of matrix required to contain the rotated cube
    # +2 just for good measure
    sample_nx = np.max(rotated_cube)
    sample_nx = int(cube_nx * 2)

    # make a matrix that will contain the rotated cube, sample
    sample = np.zeros((sample_nx, sample_nx, sample_nx))

    # loop through each of the rotated points
    # make those points =1 in the sample matrix
    for i in range(rotated_cube.shape[0]):
        sample[rotated_cube[i, 0], rotated_cube[i, 1], rotated_cube[i, 2]] = 1

    # if we took a slice of the sample matrix, the rotated cube projection looks like swiss cheese
    # use binary dilation of fill the holes in the matrix
    for i in range(sample_nx):
        sample_im = sample[:, :, i]

        scipym.binary_dilation(sample_im, structure=np.ones((2, 2)), output=sample_im)

        sample[:, :, i] = sample_im

        if save:
            imageio.imwrite("sample_slices\\slice_" + str(i) + ".tiff", sample_im)

    return sample


def make_n_slices(sample, n=10, save=False, fname="slice_chunk"):
    # number of arrays in z we have
    sample_z_size = sample.shape[2]

    # number of arrays in z for every slice
    slice_size = int(sample_z_size / n)

    # Extra work to pad the last array (not required yet)
    # remainder = sample_z_size%n
    # pad = np.zeros( (sample.shape[0], sample.shape[1], slice_size-remainder))

    # empty array where we will put our slices
    slices = np.zeros((sample.shape[0], sample.shape[1], n))

    # for every slice
    for i in range(n):
        # take the relevant z arrays from sample
        chunk = sample[:, :, slice_size * i : slice_size * i + slice_size]
        # sum them to make 1 2d array
        chunk_2d = np.sum(chunk, axis=2)
        # put that in our empty matrix
        slices[:, :, i] = chunk_2d

        # save the output to tiff
        if save:
            imageio.imwrite(
                "sample_slices\\" + fname + str(i) + ".tiff", chunk_2d.astype(int)
            )

    return slices


def movie_slices(movie):
    # generates a movie plot through the slices
    generation_ims = []
    fig = plt.figure()
    for i in range(movie.shape[2]):
        movie_frame = movie[:, :, i]

        generation_im = plt.imshow(movie_frame, animated=True)

        generation_ims.append([generation_im])

    ani = animation.ArtistAnimation(fig, generation_ims, interval=1, blit=False)

    plt.show()
    return generation_ims


def stamp_cube_into_sample(cube, sample):
    # init an empty sample array
    new_sample = np.zeros(sample.shape)

    # x,y,z location of where we will stamp the cube into sample
    stamp_loc_xs = rn.randint(2, sample.shape[0] - cube.shape[0] - 2)
    stamp_loc_ys = rn.randint(2, sample.shape[1] - cube.shape[1] - 2)
    stamp_loc_zs = rn.randint(2, sample.shape[2] - cube.shape[2] - 2)

    # stamp the cube into the new sample array
    new_sample[
        stamp_loc_xs : stamp_loc_xs + cube.shape[0],
        stamp_loc_ys : stamp_loc_ys + cube.shape[1],
        stamp_loc_zs : stamp_loc_zs + cube.shape[2],
    ] = cube

    # for every place where new sample=0, make it equal to what the old sample was
    new_sample = np.where(new_sample == 0, sample, new_sample)

    return new_sample


sample = np.zeros((1000, 1000, 500))

num_cubes = 10
for i in range(num_cubes):
    angle = 2 * np.pi * rn.random()
    axis = [rn.randint(1, 3) for i in range(3)]

    print("Cube:", str(i + 1), "/", str(num_cubes))
    print("\tAngle: ", str(angle * 180 / np.pi), sep="\t")
    print("\tRotation Axis: ", str(axis), sep="\t")

    cube_stamp = make_cube_stamp(20, axis=axis, angle=angle)
    sample = stamp_cube_into_sample(cube_stamp, sample)


sliced = make_n_slices(sample, sample.shape[2])  # int(sample.shape[2]/10))

generation_ims = []
fig = plt.figure()
for i in range(sliced.shape[2]):
    movie_frame = sliced[:, :, i]

    generation_im = plt.imshow(movie_frame, animated=True)

    generation_ims.append([generation_im])

ani = animation.ArtistAnimation(fig, generation_ims, interval=100, blit=False)

plt.show()
