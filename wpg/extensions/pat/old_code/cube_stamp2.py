import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import scipy.ndimage.morphology as scipym
import random as rn
import imageio

from pyquaternion import Quaternion


class Cube:

    # A compound cube class object that creates an x-ray scattering event

    def __init__(self, length, sample):

        self.length = length  # side length of the cube [px]
        self.sample = sample  # the sample area that is attached to the cube
        self.qs = self.make_quart()  # the quarternion that defines the random rotation for the cube
        self.stamp = self.make_stamp()  # the 3d stamp of the cube area thats placed in the sample object
        self.loc = self.make_loc(self.sample)  # loctation of the cube in the sample space [px coordinates]

    def make_quart(self):

        angle = 2 * np.pi * rn.random()  # make a random angle
        axis = [rn.randint(1, 3) for i in range(3)]  # random axis of rotation [111] to [333]
        qs = Quaternion(axis=axis, angle=angle)  # make the quarternion
        return qs

    def grow(self, amnt):

        self.length += amnt  # increase the cube length by amnt
        self.stamp = self.make_stamp()  # remake the stamp

    def make_stamp(self):

        cube_stamp = np.ones((self.length, self.length, self.length))  # make the initial un rotated cube
        cube_points = np.argwhere(cube_stamp == 1)  # list the points in the cube
        rotated_cube_points = np.dot(cube_points,
                                     self.qs.rotation_matrix)  # rotate the list of points about the quarternion

        # translate all the points to positive and become integers
        for i in range(3):
            rotated_cube_points[:, i] -= np.min(rotated_cube_points[:, i])

        rotated_cube_points = rotated_cube_points.astype(int)

        # taking the max value of the rotated cube will give the size of matrix required to contain the rotated cube
        # +10 and +5 allow for a 5 pixel border surrounding the cube
        rotated_cube_nx = np.max(rotated_cube_points) + 10
        rotated_cube_points += 5

        # make a matrix that will contain the rotated cube, rot_cube
        rot_cube = np.zeros((rotated_cube_nx, rotated_cube_nx, rotated_cube_nx))

        # loop through each of the rotated points
        # make those points =1 in the rot_cube matrix
        for i in range(rotated_cube_points.shape[0]):
            rot_cube[rotated_cube_points[i, 0], rotated_cube_points[i, 1], rotated_cube_points[i, 2]] = 1

        # if we took a slice of the sample matrix, the rotated cube projection looks like swiss cheese
        # use binary dilation of fill the holes in the matrix
        for i in range(rotated_cube_nx):
            sample_im = rot_cube[:, :, i]

            scipym.binary_dilation(sample_im, structure=np.ones((2, 2)), output=sample_im)

            rot_cube[:, :, i] = sample_im

        self.stamp = rot_cube
        return self.stamp

    def make_loc(self, sample):

        sample.cubes.append(self)  # connect the cube and the sample object
        # stamp location index of the top left corner facing the user
        stamp_loc_x = rn.randint(sample.growth_border, self.sample.nx - self.length - sample.growth_border)
        stamp_loc_y = rn.randint(sample.growth_border, self.sample.ny - self.length - sample.growth_border)
        stamp_loc_z = rn.randint(sample.growth_border, self.sample.nz - self.length - sample.growth_border)

        # index coordinates of the stamp in the sample
        loc = (stamp_loc_x, stamp_loc_y, stamp_loc_z)

        self.loc = loc
        return loc


class Sample:
    def __init__(self, nx, ny, nz, growth_border=50):
        self.nx = nx  # number of pixels in x,y,z
        self.ny = ny
        self.nz = nz
        self.cubes = []  # list of cube objects attached to the sample space
        self.volume = np.zeros((nx, ny, nz))  # initalize a volume for the sample space
        # define a border for which no cubes can start in (allows space for growth)
        self.growth_border = growth_border

    def make_array(self):
        print('Rendering cubes into sample')
        # initailze the stamp volume
        self.volume = np.zeros((self.nx, self.ny, self.nz))
        
        count=0
        # render each cube in the sample volume
        for cube in self.cubes:
            
            print('Rendering Cube:', str(count),'/',str(len(self.cubes)))
            count=count+1
            # create a space to only include the current cube
            new_volume = np.zeros((self.nx, self.ny, self.nz))

            # stamp diminsions locations to make uniform expansion about center of mass
            stamp_x_start = np.round(cube.loc[0] - np.shape(cube.stamp)[0] / 2).astype(int)
            stamp_y_start = np.round(cube.loc[1] - np.shape(cube.stamp)[1] / 2).astype(int)
            stamp_z_start = np.round(cube.loc[2] - np.shape(cube.stamp)[2] / 2).astype(int)

            stamp_x_fin = stamp_x_start + np.shape(cube.stamp)[0]
            stamp_y_fin = stamp_y_start + np.shape(cube.stamp)[1]
            stamp_z_fin = stamp_z_start + np.shape(cube.stamp)[2]

            # stamp the cube into the blank volume space
            try:
                new_volume[stamp_x_start:stamp_x_fin,
                stamp_y_start:stamp_y_fin,
                stamp_z_start:stamp_z_fin,
                ] = cube.stamp

            except ValueError:
                # Catch an error if cube tries to be stamped outside of the blank volume
                print('ERROR: Trying to expand cube outside of sample volume')
                print('Try increaseing growth border or use less growth steps')
                exit(1)

            # Where the blank volume is 0, there is no new cube. we can safely put the value that was in old sample
            # volume into these indices
            new_volume = np.where(new_volume == 0, self.volume, new_volume)
            # save the volume
            self.volume = new_volume
        # if you want to output the array, do x=sam.make_array()
        # if you just want to update, do sam.make_array()
        return new_volume

    def grow(self, amnt):
        #grow every cube connected to the sample volume
        for cube in self.cubes:
            cube.grow(amnt)

    def save_to_tiff(self, n, fname=''):


        print('saving array')
        # number of arrays in z for every slice
        slice_size = int(self.nz / n)

        # empty array where we will put our slices
        slices = np.zeros((self.nx, self.ny, n))

        # for every slice
        for i in range(n):
            # take the relevant z arrays from sample
            chunk = self.volume[:, :, slice_size * i:slice_size * i + slice_size]
            # sum them to make 1 2d array
            chunk_2d = np.sum(chunk, axis=2)
            # put that in our empty matrix
            slices[:, :, i] = chunk_2d

            # save the output to tiff
            imageio.imwrite(fname + str(i)+'.TIFF',  chunk_2d.astype(np.float32))

        return slices


#   'fully grown' iron crystal is about 200nm
#   for a fully grown crystal to take up 20 pixel, thats 10nm per pixel.
resolution = 10e-9

## inputs
seed_int = 30
sample_nx = 256
sample_ny = 256
sample_nz = 256
growth_border = 50
num_cubes =25
init_cube_size = 5
grow_amnt = 45
fname = None
num_tiffs = 5

rn.seed(seed_int)
# make the sample
print('Making the sample')



print('Making the cubes')
# make the cubes
for i in range(num_cubes):
    print('cube:', str(i), '/', str(num_cubes))
    c1 = Cube(init_cube_size, sam)


#make the sample array to image it
for i in range(0,grow_amnt,1):
    print('grown\t', str(i))
    sam.make_array()
    file_name = str(i) +'samples_256_pix_25cubes' 
    sam.save_to_tiff(1,'sample_slices/'+file_name)
    
    sam.grow(1)
##sum to get 2d projection
#sam_ave = np.sum(sam.volume, 2)
##plot the projection
#plt.figure()
#plt.imshow(sam_ave)


##make the sample array to image it

##sum to get 2d projection
##sum to get 2d projection
#sam_ave = np.sum(sam.volume, 2)
##plot the projection
#plt.figure()
#plt.imshow(sam_ave)


#####Animation stuff
#
#
# # emtpy array of frames for the move (list of 2d image stack in z direction)
# movie_frames = np.zeros((sample_nx, sample_ny, 25))
#
# for i in range(np.shape(movie_frames)[2]):  # for every frame in the movie
#
#     print('growth step', str(i))  # print the current frame
#
#     sam.make_array()  # render the cubes in the sample volume
#
#     sam_ave = np.sum(sam.volume, 2)  # sum through the z axis to get the projection
#     movie_frames[:, :, i] = sam_ave  # insert the projection as a movie frame
#     sam.grow(1)  # grow the cubes in the sample plane
#
# # code to write an animation to the figure
# generation_ims = []
# fig = plt.figure()
#
# for i in range(np.shape(movie_frames)[2]):
#     movie_frame = movie_frames[:, :, i]
#
#     generation_im = plt.imshow(movie_frame, animated=True)
#
#     generation_ims.append([generation_im])
#
# ani = animation.ArtistAnimation(fig, generation_ims, interval=100, blit=False)
#
# plt.show()
