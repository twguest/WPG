

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import scipy.ndimage.morphology as scipym
import random as rn
import imageio
from wpg.srwl_uti_smp import srwl_opt_setup_transm_from_file
import pickle 



from pyquaternion import Quaternion

rn.seed(1000)

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
        
        while True: #added to avoid  division by zero errors when axis is [0,0,0]
            axis = [rn.randint(0, 3) for i in range(3)]  # random axis of rotation [111] to [333]
            mag_sq = np.dot(axis, axis)
            if mag_sq > 0:
                break
      
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
        # initailze the stamp volume
        self.volume = np.zeros((self.nx, self.ny, self.nz))

        # render each cube in the sample volume
        i=0
        for cube in self.cubes:
            
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
            
            else:
                print('Cube {} placed in array'.format(i))
            
            i +=1

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

    def save_to_tiff(self, n, fname=None):


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
            if type(fname) == str:
                imageio.imwrite(fname + str(i) + '.tiff', chunk_2d.astype(int))

        return slices





## inputs
pxnm = 1 #number of pixels per nm

num_cubes =10
init_cube_size = 200*pxnm
growth_border = init_cube_size
sample_nx = 3000*pxnm #was 1000
sample_ny = 1500*pxnm  #was 200
sample_nz =int(2*growth_border+((init_cube_size+num_cubes)*2**0.5))*pxnm 



sam = Sample(sample_nx, sample_ny, sample_nz, growth_border)

for i in range(num_cubes):
    c1 = Cube(init_cube_size, sam)
    print ('Cube {} of {} done.'.format(i, num_cubes))



x = sam.make_array()

x= np.sum(x,axis=2)

maxThickness = np.max(x)
print('Max. thickness = {}'.format(maxThickness))

# create optical element transmission function
hv = 7374.0
delta = 2.51144647e-5 # CXRO
attenLength = 3.39642e-6  # CXRO
thickness =  303.e-9


tr = srwl_opt_setup_transm_from_file('/opt/wpg/wpg/extensions/pat/io/sample/sample.tiff',
                    resolution=pxnm*1e-9,
                    thickness=maxThickness,  #this is max thickness
                    delta=delta,
                    atten_len=attenLength,
                    xc=0.0, yc=0.0,
                    area=None,
                    rotate_angle=0,
                    rotate_reshape=False,
                    cutoff_background_noise=None,
                    background_color=0,
                    tile=None,
                    shift_x=0,
                    shift_y=0,
                    invert=False
                    )


pickle.dump( tr, open( "/opt/wpg/wpg/extensions/pat/io/sample/tr.p", "wb" ) )

im =Image.fromarray(x)
im.save('/opt/wpg/wpg/extensions/pat/io/sample/sample.tiff')


plt.figure()
plt.imshow(x)#, extent=[3,-3,15,-15])
plt.show()
