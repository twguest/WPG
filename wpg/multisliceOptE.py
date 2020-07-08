#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 21:05:32 2019
Upated 18/9/2019.  

@author: gvanriessen
"""

import math
import types

# from wpg.wpg_uti_oe import show_transmission

from wpg.beamline import Beamline
from wpg.wavefront import Wavefront
from wpg.useful_code.wfrutils import print_beamline

from wpg.srwlib import SRWLOptD, SRWLOptC, SRWLOptA, SRWLOptT

from wpg.gvrutils import propagationParameters, plotWavefront, plot_wfront2
from wpg.gvrutils import (
    writeIntensity,
    pixelScale,
    get_transmissionFunctionCplx,
    get_intensity_on_axis,
)
from wpg.gvrutils import show_transmissionFunctionCplx

from wpg.optical_elements import Use_PP


import numexpr as ne
from PIL import Image
from beautifultable import BeautifulTable
import timeit
import pickle

use_gpu = False

if use_gpu:
    import afnumpy as np
    import afnumpy.fft as fft

    use = "afnumpy/GPU"
else:
    import numpy as np
    import numpy.fft as fft
use = "numpy/CPU"

import os

try:
    from wpg import srwlpy as srwl
except ImportError:
    import srwlpy as srwl  #  Hack for read the docs

from wpg.uti_io import *


def blFromOE(elements, parameters, description=None):
    # construct beamline
    bl = Beamline()  # (description=description)

    for e, p in zip(elements, parameters):
        bl.append(e, p)

    return bl


# ********************** The class for Samples:
class SRWLUtiSmp:
    """MOdified class for Samples from image file (e.g., .tif), NumPy array (.npy), etc.
    !! still restricted to 8-bit representation in input source

    :param inputData: full path to the image or the saved NumPy array,
                      OR a numpy array (to avoid file IO)

    :param area: the coordinates of the rectangle area listed in the following order: x_start, x_end, y_start, y_end.
    :param rotate_angle: the angle [deg] to rotate the read image counterclockwise. See scipy.ndimage.interpolation.rotate() for details.
    :param rotate_reshape: if reshape is true, the output shape is adapted so that the input array is contained completely in the output.
    :param cutoff_background_noise: the ratio for cutoff the background noise (between 0 and 1).
    :param background_color: the background color code to use instead of the background noise (0=black, 255=white).
    :param tile: the list/tuple (rows, columns) to tile the cut area of the image. See numpy.tile() for details.
    :param shift_x: shift the whole image horizontally. Positive value shifts the image to the right, negative - to the left. See numpy.pad() for details.
    :param shift_y: shift the whole image vertically. Positive value shifts the image to the top, negative - to the bottom. See numpy.pad() for details.
    :param invert: invert the image. See numpy.invert() for details.
    :param is_show_images: a flag to show the initial and processed images.
    :param is_save_images: a flag to save the initial and processed images.
    :param raw_image_name: the name of the raw file in case if it's saved.
    :param processed_image_name: the name of the processed file in case if it's saved.
    :param prefix: the prefix to add to the names of the saved image files.
    :param output_image_format: the format of the output file. If not specified, the input format is used.
    """

    def __init__(
        self,
        inputData,
        area=None,
        rotate_angle=0,
        rotate_reshape=True,
        cutoff_background_noise=0.25,
        background_color=0,
        tile=None,
        shift_x=None,
        shift_y=None,
        invert=None,
        invert_bitwise=None,
        is_show_images=True,
        is_save_images=True,
        raw_image_name="raw",
        processed_image_name="processed",
        prefix="",
        output_image_format=None,
        crop=None,
        pad=None,
    ):
        # Input parameters:
        self.inputData = inputData
        self.area = area
        self.crop = crop
        self.rotate_angle = rotate_angle if rotate_angle is not None else 0
        self.rotate_reshape = True if rotate_reshape else False
        self.cutoff_background_noise = (
            cutoff_background_noise if cutoff_background_noise is not None else 0
        )
        self.background_color = background_color if background_color is not None else 0
        self.invert = invert
        self.invert_bitwise = invert_bitwise
        self.tile = tile
        self.shift_x = shift_x
        self.shift_y = shift_y
        self.is_show_images = is_show_images
        self.is_save_images = is_save_images
        output_image_format = (
            os.path.splitext(inputData)[1].replace(".", "")
            if not output_image_format
            else output_image_format
        )
        self.raw_image_name = self._add_prefix(
            prefix, raw_image_name, output_image_format
        )
        self.processed_image_name = self._add_prefix(
            prefix, processed_image_name, output_image_format
        )
        self.padding = pad

        # Output parameters:
        self.data = None
        self.raw_image = None
        self.processed_image = None
        self.nx = None
        self.ny = None
        self.limit_value = None

        # Check input type automatically:
        self.input_type = self._check_input_type()

        # Set the dir where to save the files:
        self.save_dir = os.path.abspath(os.path.dirname(inputData))

        # Check the input file(s):
        if self.input_type != np.ndarray:
            self._check_files()

            # Process input:
        self.read_sample()

        # print ('DATA SHAPE')
        # print (self.data.shape)

        self.preprocess_input()

        # Show the resulted images:
        if self.is_show_images:
            self.show_images()

        # Save the resulted images:
        if self.is_save_images:
            self.save_images()

    def preprocess_input(self):
        try:
            from PIL import Image
        except ImportError:
            raise ImportError(import_err_msg.format("PIL", "pillow"))

        # Remove background noise:
        assert (
            0 <= self.background_color <= 255
        ), "Background color ({}) should be between 0 and 255".format(
            self.background_color
        )
        assert (
            0 <= self.cutoff_background_noise <= 1
        ), "Cutoff background noise ({}) should be between 0 and 1".format(
            self.cutoff_background_noise
        )
        # self.data[np.where(self.data < self.limit_value * self.cutoff_background_noise)] = np.uint16(
        #   self.background_color)

        if self.area:
            self.crop()

        if self.padding is not None:
            self.pad()

        if self.tile:
            self.tile()

        if self.shift_x:
            assert (
                type(self.shift_x) is int
            ), "Type of shift_x ({}) should be int".format(type(self.shift_x).__name__)
            if self.shift_x > 0:
                self.data = np.pad(self.data, ((0, 0), (self.shift_x, 0)), mode="edge")[
                    :, : -self.shift_x
                ]
                # self.data = np.pad(self.data, ((0, 0), (self.shift_x, 0)), mode='constant',
                #                   constant_values=(self.background_color))[:, :-self.shift_x]
                # self.data[:,-self.shift_x:] = self.background_color # unnecessary?, but included for testing
            else:
                self.data = np.pad(
                    self.data, ((0, 0), (0, -self.shift_x)), mode="edge"
                )[:, -self.shift_x :]
                # self.data = np.pad(self.data, ((0, 0), (0, -self.shift_x)), mode='constant',
                #                   constant_values=(self.background_color))[:, -self.shift_x:]
                # self.data[:,:-self.shift_x]  = self.background_color  # unnecessary?, but included for testing

        if self.shift_y:
            assert (
                type(self.shift_y) is int
            ), "Type of shift_y ({}) should be int".format(type(self.shift_y).__name__)
            if self.shift_y < 0:
                self.data = np.pad(
                    self.data,
                    ((-self.shift_y, 0), (0, 0)),
                    mode="constant",
                    constant_values=(self.background_color),
                )[: self.shift_y, :]
            else:
                self.data = np.pad(
                    self.data,
                    ((0, self.shift_y), (0, 0)),
                    mode="constant",
                    constant_values=(self.background_color),
                )[self.shift_y :, :]

        if self.invert_bitwise:
            # Originally 'invert'. Implemewnt bit-wise NOT of the underlying binary representation of the integers in the input arrays.
            self.data = np.invert(self.data)

        if self.invert:
            # invert a greyscale image (swapping black and white)
            self.data = self.limit_value - self.data

        self.processed_image = Image.fromarray(self.data)

    def crop(self):

        assert (
            type(self.area) in [tuple, list] and len(self.area) == 4
        ), "The area should be a list/tuple and contain 4 elements (x_start, x_end, y_start, y_end)"
        assert (
            self.area[0] <= self.data.shape[1]
            and self.area[1] <= self.data.shape[1]
            and self.area[2] <= self.data.shape[0]
            and self.area[3] <= self.data.shape[0]
        ), "x_start ({}) and x_end ({}) should be less than {} y_start ({}) and y_end ({}) should be less than {}".format(
            self.area[0],
            self.area[1],
            self.data.shape[1],
            self.area[2],
            self.area[3],
            self.data.shape[0],
        )
        assert (
            self.area[0] < self.area[1] and self.area[2] < self.area[3]
        ), "x_start ({}) should be less than x_end ({}) y_start ({}) should be less than y_end ({})".format(
            self.area[0], self.area[1], self.area[2], self.area[3],
        )
        self.data = self.data[self.area[2] : self.area[3], self.area[0] : self.area[1]]

    def pad(self):
        assert (
            type(self.padding) in [tuple, list] and len(self.padding) == 4
        ), "The pad parameters should be a list/tuple and contain 4 elements (top,bottom,left, right,)"
        assert (
            self.padding[0] >= 0
            and self.padding[1] >= 0
            and self.padding[2] >= 0
            and self.padding[3] >= 0
        ), "Pad values must be >= 0"
        self.data = np.pad(
            self.data,
            [(self.padding[0], self.padding[1]), (self.padding[2], self.padding[3])],
            "edge",
        )
        #'constant', constant_values=self.background_color)

    def tile(self):
        assert type(self.tile) in [
            list,
            tuple,
        ], "The type of tile ({}) should be list/tuple".format(type(self.tile).__name__)
        assert (
            len(self.tile) == 2
        ), "The size ({}) of the list/tuple should be 2".format(len(self.tile))
        self.data = np.tile(self.data, self.tile)

        if self.rotate_angle:
            assert (
                -360 < self.rotate_angle < 360
            ), "The angle should be from -360 to 360 degrees."
            self.data = rotate(
                self.data,
                self.rotate_angle,
                axes=(1, 0),
                reshape=self.rotate_reshape,
                output=None,
                order=0,
                mode="constant",
                cval=self.background_color,
                prefilter=False,
            )

    def get_data_from_image(self):
        """
        import_err_msg = '"{}" library cannot be imported. Please install it first with "pip install {}".'

        try:
            from PIL import Image
        except ImportError:
            raise ImportError(import_err_msg.format('PIL', 'pillow'))
        try:
            from scipy.ndimage.interpolation import rotate
        except ImportError:
            raise ImportError(import_err_msg.format('SciPy', 'scipy'))
        """

        d = read_image(self.inputData)
        self.data = d["data"]
        self.raw_image = d["raw_image"]
        self.limit_value = d["limit_value"]

        print("Limit Value set to {}".format(self.limit_value))

    def get_image_from_data(self):
        data = np.load(self.inputData)
        self.limit_value = 255
        self.data = (data - data.min()) / (data.max() - data.min()) * self.limit_value
        self.raw_image = Image.fromarray(self.data.astype(np.uint8))

    def get_image_from_np_array(self):

        self.limit_value = 255
        self.data = (data - data.min()) / (data.max() - data.min()) * self.limit_value
        self.raw_image = Image.fromarray(self.data.astype(np.uint8))

    def read_sample(self):
        if self.input_type == "image":
            self.get_data_from_image()
        elif self.input_type == "npy":
            self.get_image_from_data()
        elif self.input_type == np.ndarray:
            self.get_image_from_np_array()
        else:
            raise NotImplementedError(
                'Processing of the "{}" input type is not implemented yet.'.format(
                    self.input_type
                )
            )
        self.nx = self.data.shape[1]
        self.ny = self.data.shape[0]

    def save_images(self):
        # self.raw_image.save(os.path.join(self.save_dir, self.raw_image_name))
        self.processed_image.save(
            os.path.join(self.save_dir, self.processed_image_name)
        )
        print("Wrote: " + os.path.join(self.save_dir, self.processed_image_name))

    def show_images(self):
        self.raw_image.show()
        self.processed_image.show()

    def _add_prefix(self, prefix, name, image_format):
        output_name = "{}_{}".format(prefix, name) if prefix else name
        return "{}.{}".format(output_name, image_format)

    def _check_files(self):
        if not os.path.isfile(self.inputData):
            raise ValueError(
                'Provided file "{}" does not exist.'.format(self.inputData)
            )

    def _check_input_type(self):
        if type(self.inputData) is np.ndarray:
            return np.ndarray
        else:
            return self._check_input_file_type()

    def _check_input_file_type(self):
        self.possible_extensions = {
            "image": ["tif", "tiff", "png", "bmp", "gif", "jpg", "jpeg"],
            "npy": ["npy"],
        }
        extension = os.path.splitext(self.inputData)[1][1:].lower()
        for k in self.possible_extensions.keys():
            for e in self.possible_extensions[k]:
                if extension == e:
                    return k
        all_extensions = [
            x for x_list in self.possible_extensions.values() for x in x_list
        ]
        all_extensions += [x.upper() for x in all_extensions]
        raise ValueError(
            "Incorrect extension: {}. Possible values: {}.".format(
                extension, ", ".join(all_extensions)
            )
        )


# ********************** Create transmission element from the data from an image file:
def srwl_opt_setup_transm_from_file(
    inputData,
    resolution,
    thickness,
    delta,
    atten_len,
    arTr=None,
    extTr=0,
    fx=1e23,
    fy=1e23,
    xc=0,
    yc=0,
    ne=1,
    e_start=0,
    e_fin=0,
    area=None,
    rotate_angle=None,
    rotate_reshape=None,
    cutoff_background_noise=None,
    background_color=None,
    tile=None,
    shift_x=None,
    shift_y=None,
    invert=None,
    is_save_images=False,
    prefix="",
    output_image_format=None,
    pad=None,
):
    """Setup Sample element.

    :param inputData: path to the input file (image or .npy).
    :param resolution: resolution of the image [m/pixel].
    :param thickness: thickness of the sample [m].
    :param delta: refractive index decrement.
    :param atten_len: attenuation length [m].
    :param arTr: complex C-aligned data array (of 2*ne*nx*ny length) storing amplitude transmission and optical path difference as function of transverse coordinates.
    :param extTr: transmission outside the grid/mesh is zero (0), or it is same as on boundary (1).
    :param fx: estimated focal length in the horizontal plane [m].
    :param fy: estimated focal length in the vertical plane [m].
    :param xc: horizontal coordinate of center [m].
    :param yc: vertical coordinate of center [m].
    :param ne: number of transmission data points vs photon energy.
    :param e_start: initial photon energy [eV].
    :param e_fin: final photon energy [eV].
    :param area: the coordinates of the rectangle area listed in the following order: x_start, x_end, y_start, y_end.
    :param rotate_angle: the angle [deg] to rotate the read image counterclockwise. See scipy.ndimage.interpolation.rotate() for details.
    :param rotate_reshape: if reshape is true, the output shape is adapted so that the input array is contained completely in the output.
    :param cutoff_background_noise: the ratio for cutoff the background noise (between 0 and 1).
    :param background_color: the background color code to use instead of the background noise (0=black, 255=white).
    :param tile: the list/tuple (rows, columns) to tile the cut area of the image. See numpy.tile() for details.
    :param shift_x: shift the whole image horizontally. Positive value shifts the image to the right, negative - to the left. See numpy.pad() for details.
    :param shift_y: shift the whole image vertically. Positive value shifts the image to the top, negative - to the bottom. See numpy.pad() for details.
    :param invert: invert the image. See numpy.invert() for details.
    :param is_save_images: a flag to save the initial and processed images.
    :param prefix: the prefix to add to the names of the saved image files.
    :param output_image_format: the format of the output file. If not specified, the input format is used.
    :return: transmission (SRWLOptT) type optical element which simulates the Sample.
    """

    input_parms = {
        "type": "sample",
        "resolution": resolution,
        "thickness": thickness,
        "refractiveIndex": delta,
        "attenuationLength": atten_len,
        "horizontalCenterCoordinate": xc,
        "verticalCenterCoordinate": yc,
        "initialPhotonEnergy": e_start,
        "finalPhotonPnergy": e_fin,
        "area": area,
        "rotateAngle": rotate_angle,
        "rotateReshape": rotate_reshape,
        "cutoffBackgroundNoise": cutoff_background_noise,
        "backgroundColor": background_color,
        "tile": tile,
        "shiftX": shift_x,
        "shiftY": shift_y,
        "invert": invert,
        "outputImageFormat": output_image_format,
        "padding": pad,
    }

    s = SRWLUtiSmp(
        inputData=inputData,
        area=area,
        rotate_angle=rotate_angle,
        rotate_reshape=rotate_reshape,
        cutoff_background_noise=cutoff_background_noise,
        background_color=background_color,
        tile=tile,
        shift_x=shift_x,
        shift_y=shift_y,
        invert=invert,
        is_show_images=False,
        is_save_images=is_save_images,
        prefix=prefix,
        output_image_format=output_image_format,
        pad=pad,
    )

    # Input parameters to SRWLOptT:
    nx = s.nx
    ny = s.ny

    print("image dimensions: %s x %s" % (str(nx), str(ny)))
    rx = nx * resolution
    ry = ny * resolution

    opT = SRWLOptT(
        _nx=nx,
        _ny=ny,
        _rx=rx,
        _ry=ry,
        _arTr=arTr,
        _extTr=extTr,
        _Fx=fx,
        _Fy=fy,
        _x=xc,
        _y=yc,
        _ne=ne,
        _eStart=e_start,
        _eFin=e_fin,
    )
    opT.pathInBody = []
    data = s.data

    # Same data alignment as for wavefront: outmost loop vs y, inmost loop vs e
    offset = 0
    for iy in range(ny):
        for ix in range(nx):
            for ie in range(ne):
                # In images Y=0 corresponds from upper-left corner, in SRW it's lower-left corner:
                pathInBody = thickness * data[ny - iy - 1, ix] / s.limit_value
                opT.pathInBody.append(pathInBody)
                opT.arTr[offset] = math.exp(
                    -0.5 * pathInBody / atten_len
                )  # amplitude transmission
                opT.arTr[offset + 1] = -delta * pathInBody  # optical path difference
                offset += 2

    opT.input_parms = input_parms

    return opT


class greyscaleToSlices:

    """MOdified class for Samples from image file (e.g., .tif), NumPy array (.npy), etc.
    !! still restricted to 8-bit representation in input source

    :param inputData: full path to the image or the saved NumPy array,
                      OR a numpy array (to avoid file IO)
    : param slices - the number of slices
    :param area: the coordinates of the rectangle area listed in the following order: x_start, x_end, y_start, y_end.
    :param rotate_angle: the angle [deg] to rotate the read image counterclockwise. See scipy.ndimage.interpolation.rotate() for details.
    :param rotate_reshape: if reshape is true, the output shape is adapted so that the input array is contained completely in the output.
    :param cutoff_background_noise: the ratio for cutoff the background noise (between 0 and 1).
    :param background_color: the background color code to use instead of the background noise (0=black, 255=white).
    :param tile: the list/tuple (rows, columns) to tile the cut area of the image. See numpy.tile() for details.
    :param shift_x: shift the whole image horizontally. Positive value shifts the image to the right, negative - to the left. See numpy.pad() for details.
    :param shift_y: shift the whole image vertically. Positive value shifts the image to the top, negative - to the bottom. See numpy.pad() for details.
    :param invert: invert the image. See numpy.invert() for details.
    :param is_show_images: a flag to show the initial and processed images.
    :param is_save_images: a flag to save the initial and processed images.
    :param raw_image_name: the name of the raw file in case if it's saved.
    :param processed_image_name: the name of the processed file in case if it's saved.
    :param prefix: the prefix to add to the names of the saved image files.
    :param output_image_format: the format of the output file. If not specified, the input format is used.
    
    
    Changes:
        included propagationparametes as optional argument.  Note that the typically short propagation distances
        between slices allow agressive zoom, but it may be more sensible to apply to wavefront before driftvol propagation.
    """

    def __init__(
        self,
        inputData,
        propagationParameters=None,
        slices=1,
        area=None,
        rotate_angle=0,
        rotate_reshape=True,
        cutoff_background_noise=0,
        background_color=0,
        tile=None,
        shift_x=None,
        shift_y=None,
        invert=None,
        invert_bitwise=None,
        is_show_images=False,
        is_save_images=False,
        raw_image_name="raw",
        processed_image_name="processed",
        prefix="",
        output_image_format=None,
        crop=None,
        pad=None,
        resolution=[10.0e-9, 10.0e-9, None],
        thickness=657.85e-9,  # total thickness of object, i.e. sum of thickness of each slice (m)
        delta=0.0031406,  # delta value of material in slices
        attenLength=0.277984e-6,  # attenuation length in material in slices
    ):

        # Input parameters:

        self.inputData = inputData
        self.greyLevels = slices
        self.area = area
        self.crop = crop
        self.rotate_angle = rotate_angle if rotate_angle is not None else 0
        self.rotate_reshape = True if rotate_reshape else False
        self.cutoff_background_noise = (
            cutoff_background_noise if cutoff_background_noise is not None else 0
        )
        self.background_color = background_color if background_color is not None else 0
        self.invert = invert
        self.invert_bitwise = invert_bitwise
        self.tile = tile
        self.shift_x = shift_x
        self.shift_y = shift_y
        self.is_show_images = is_show_images
        self.is_save_images = is_save_images
        output_image_format = (
            os.path.splitext(inputData)[1].replace(".", "")
            if not output_image_format
            else output_image_format
        )
        self.raw_image_name = self._add_prefix(
            prefix, raw_image_name, output_image_format
        )
        self.processed_image_name = self._add_prefix(
            prefix, processed_image_name, output_image_format
        )
        self.padding = pad

        # for now we assume same material in every slice.  Improve in next version
        self.resolution = resolution
        self.thickness = thickness  # total thickness of object, i.e. sum of thickness of each slice (m)
        self.delta = delta  # delta value of material in slices
        self.attenLength = attenLength  # attenuation length in material in slice

        if propagationParameters is None:
            self.pp = use_PP()
        else:
            self.pp = propagationParameters

        # Output parameters:
        self.data = None
        self.raw_image = None
        self.processed_image = None
        self.nx = None
        self.ny = None
        self.limit_value = None

        # Check input type automatically:
        self.input_type = self._check_input_type()

        # Set the dir where to save the files:
        self.save_dir = os.path.abspath(os.path.dirname(self.inputData))

        # Check the input file(s):
        if self.input_type != np.ndarray:
            self._check_files()

            # Process input:
        self.read_sample()

        self.preprocess_input()

        self.slicer()

        self.multiSliceTr()

        # Show the resulted images:
        if self.is_show_images:
            self.show_images()

        # Save the resulted images:
        if self.is_save_images:
            self.save_images()

    def multiSliceTr(self):

        material = [
            (self.thickness / self.greyLevels / self.greyLevels),
            self.delta,
            self.attenLength,
        ]

        self.trSlices = []
        for s, i in zip(self.sliceFiles, range(0, self.greyLevels)):

            trans = [0, 0, i * self.thickness / self.greyLevels, 0]

            self.trSlices.append(
                [s, trans, self.resolution, material, self.padding, self.pp]
            )

        self.tr = multisliceOptE(objectSlices=self.trSlices)

        self.tr.printSliceParams()

    def propagate(self, wf=None):

        self.tr.propagateSlices(wf)

    def slicer(self):
        try:
            from PIL import Image
        except ImportError:
            raise ImportError(import_err_msg.format("PIL", "pillow"))

        import matplotlib.pyplot as plt
        import copy
        #from pymp.shared import array as psa
        import pymp as pymp

        gs = copy.deepcopy(self.data)
        s = copy.deepcopy(self.data)

        greyLevelUpper = np.round(
            np.linspace(self.limit_value, 0, self.greyLevels + 1)[0 : self.greyLevels]
        )
        greyLevelLower = np.round(
            np.linspace(self.limit_value, 0, self.greyLevels + 1)[
                1 : self.greyLevels + 1
            ]
        )
        # index=range(0,self.greyLevels)
        # greyLevelRanges = zip(greyLevelLower,greyLevelUpper,index)
        greyLevelRanges = []
        for i in range(self.greyLevels):
            greyLevelRanges.append([greyLevelLower[i], greyLevelUpper[i], i])

        parproc = False
        if parproc:
            NumProcesses = 10
            self.sliceFiles = psa(len(list(greyLevelRanges)))  # , dtype=np.string_)
            with pymp.Parallel(NumProcesses) as p:
                for index in p.range(0, len(list(greyLevelRanges)) - 1):
                    p.print("Slicing, slice number %d" % index)
                    gs[
                        (gs > rl) & (gs <= rh)
                    ] = rl  # object minus slice, i.e. what is left of the cake
                    s = s - gs  # the slice
                    img = Image.fromarray(s.astype(np.uint8))
                    outPath = os.path.join(
                        self.save_dir, self.processed_image_name + str(index) + ".tif"
                    )
                    img.save(outPath)
                    self.sliceFiles[index] = outPath

        else:  # no parallel processing, and display of progress

            # self.slice = []
            self.sliceFiles = []
            for rl, rh, i in greyLevelRanges:

                print(
                    "Slicing, slice number %d - greyvalues in interval[%d,%d]"
                    % (i, rl, rh)
                )
                s = np.where(gs > rl, gs - rl, 0)  # the slice
                print("rl is {}".format(rl))
                # s = np.subtract(s, rl)
                if i == (len(greyLevelRanges) - 1):  # last slice
                    s = np.subtract(
                        s, np.min(s[s != 0] - 2)
                    )  # correction, but why +2??, i.e. why do negative values appear in last slice
                    s[s == 255] = 0

                gs = np.where((gs > rl), rl, gs)  # the residual

                print(
                    "Greyscale values in slice in interval [%d, %d]"
                    % (np.min(s), np.max(s))
                )
                print(
                    "Greyscale values in residual object in interval [%d, %d]"
                    % (np.min(gs), np.max(gs))
                )

                # self.slice.append(s)          # append to list of slices

                plt.ion()
                fig = plt.figure()
                axCake = fig.add_subplot(121)
                axSlice = fig.add_subplot(122)
                # fig.colorbar()
                imCake = axCake.imshow(gs, cmap="gray")  # , vmin=0, vmax=rh-rl)
                imSlice = axSlice.imshow(
                    s, cmap="gray"
                )  # , vmin=0, vmax=self.limit_value)
                fig.canvas.draw()

                img = Image.fromarray(s.astype(np.uint8))
                outPath = os.path.join(
                    self.save_dir, self.processed_image_name + str(i) + ".tif"
                )
                img.save(outPath)
                print("Saved {}".format(outPath))
                self.sliceFiles.append(outPath)

    def preprocess_input(self):

        try:
            from PIL import Image
        except ImportError:
            raise ImportError(import_err_msg.format("PIL", "pillow"))

        # Remove background noise:
        assert (
            0 <= self.background_color <= 255
        ), "Background color ({}) should be between 0 and 255".format(
            self.background_color
        )
        assert (
            0 <= self.cutoff_background_noise <= 1
        ), "Cutoff background noise ({}) should be between 0 and 1".format(
            self.cutoff_background_noise
        )
        self.data[
            np.where(self.data < self.limit_value * self.cutoff_background_noise)
        ] = np.uint16(self.background_color)

        if self.area:
            self.crop()

        if self.padding:
            self.pad()

        if self.tile:
            self.tile()

        if self.shift_x:
            assert (
                type(self.shift_x) is int
            ), "Type of shift_x ({}) should be int".format(type(self.shift_x).__name__)
            if self.shift_x > 0:
                self.data = np.pad(
                    self.data,
                    ((0, 0), (self.shift_x, 0)),
                    mode="constant",
                    constant_values=(self.background_color),
                )[:, : -self.shift_x]
            else:
                self.data = np.pad(
                    self.data,
                    ((0, 0), (0, -self.shift_x)),
                    mode="constant",
                    constant_values=(self.background_color),
                )[:, -self.shift_x :]

        if self.shift_y:
            assert (
                type(self.shift_y) is int
            ), "Type of shift_y ({}) should be int".format(type(self.shift_y).__name__)
            if self.shift_y < 0:
                self.data = np.pad(
                    self.data,
                    ((-self.shift_y, 0), (0, 0)),
                    mode="constant",
                    constant_values=(self.background_color),
                )[: self.shift_y, :]
            else:
                self.data = np.pad(
                    self.data,
                    ((0, self.shift_y), (0, 0)),
                    mode="constant",
                    constant_values=(self.background_color),
                )[self.shift_y :, :]

        if self.invert_bitwise:
            # Originally 'invert'. Implemewnt bit-wise NOT of the underlying binary representation of the integers in the input arrays.
            self.data = np.invert(self.data)

        if self.invert:
            # invert a greyscale image (swapping black and white)
            self.data = self.limit_value - self.data

        self.processed_image = Image.fromarray(self.data)

    def area(self):

        assert (
            type(self.area) in [tuple, list] and len(self.area) == 4
        ), "The area should be a list/tuple and contain 4 elements (x_start, x_end, y_start, y_end)"
        assert (
            self.area[0] <= self.data.shape[1]
            and self.area[1] <= self.data.shape[1]
            and self.area[2] <= self.data.shape[0]
            and self.area[3] <= self.data.shape[0]
        ), "x_start ({}) and x_end ({}) should be less than {} y_start ({}) and y_end ({}) should be less than {}".format(
            self.area[0],
            self.area[1],
            self.data.shape[1],
            self.area[2],
            self.area[3],
            self.data.shape[0],
        )
        assert (
            self.area[0] < self.area[1] and self.area[2] < self.area[3]
        ), "x_start ({}) should be less than x_end ({}) y_start ({}) should be less than y_end ({})".format(
            self.area[0], self.area[1], self.area[2], self.area[3],
        )
        self.data = self.data[self.area[2] : self.area[3], self.area[0] : self.area[1]]

    def pad(self):
        assert (
            type(self.padding) in [tuple, list] and len(self.padding) == 4
        ), "The pad parameters should be a list/tuple and contain 4 elements (left, right, top bottom)"
        assert (
            self.padding[0] >= 0
            and self.padding[1] >= 0
            and self.padding[2] >= 0
            and self.padding[3] >= 0
        ), "Pad values must be >= 0"

        self.data = np.padding(
            self.data,
            ((self.padding[0], self.padding[1]), (self.padding[2], self.padding[3])),
            "constant",
            constant_values=self.background_color,
        )

    def tile(self):
        assert type(self.tile) in [
            list,
            tuple,
        ], "The type of tile ({}) should be list/tuple".format(type(self.tile).__name__)
        assert (
            len(self.tile) == 2
        ), "The size ({}) of the list/tuple should be 2".format(len(self.tile))
        self.data = np.tile(self.data, self.tile)

        if self.rotate_angle:
            assert (
                -360 < self.rotate_angle < 360
            ), "The angle should be from -360 to 360 degrees."
            self.data = rotate(
                self.data,
                self.rotate_angle,
                axes=(1, 0),
                reshape=self.rotate_reshape,
                output=None,
                order=0,
                mode="constant",
                cval=self.background_color,
                prefilter=False,
            )

    def get_data_from_image(self):
        """import_err_msg = '"{}" library cannot be imported. Please install it first with "pip install {}".'

        try:
            from PIL import Image
        except ImportError:
            raise ImportError(import_err_msg.format('PIL', 'pillow'))
        try:
            from scipy.ndimage.interpolation import rotate
        except ImportError:
            raise ImportError(import_err_msg.format('SciPy', 'scipy'))

        """

        d = read_image(self.inputData)
        self.data = d["data"]
        self.raw_image = d["raw_image"]
        self.limit_value = d["limit_value"]
        print("Limit value is %s" % str(self.limit_value))

    def get_image_from_data(self):
        data = np.load(self.inputData)
        self.limit_value = 255
        self.data = (data - data.min()) / (data.max() - data.min()) * self.limit_value
        self.raw_image = Image.fromarray(self.data.astype(np.uint8))

    def get_image_from_np_array(self):

        self.limit_value = 255
        self.data = (data - data.min()) / (data.max() - data.min()) * self.limit_value
        self.raw_image = Image.fromarray(self.data.astype(np.uint8))

    def read_sample(self):
        if self.input_type == "image":
            self.get_data_from_image()
        elif self.input_type == "npy":
            self.get_image_from_data()
        elif self.input_type == np.ndarray:
            self.get_image_from_np_array()
        else:
            raise NotImplementedError(
                'Processing of the "{}" input type is not implemented yet.'.format(
                    self.input_type
                )
            )
        self.nx = self.data.shape[1]
        self.ny = self.data.shape[0]

    def save_images(self):
        # self.raw_image.save(os.path.join(self.save_dir, self.raw_image_name))
        self.processed_image.save(
            os.path.join(self.save_dir, self.processed_image_name)
        )

    def show_images(self):
        self.raw_image.show()
        self.processed_image.show()

    def _add_prefix(self, prefix, name, image_format):
        output_name = "{}_{}".format(prefix, name) if prefix else name
        return "{}.{}".format(output_name, image_format)

    def _check_files(self):
        if not os.path.isfile(self.inputData):
            raise ValueError(
                'Provided file "{}" does not exist.'.format(self.inputData)
            )

    def _check_input_type(self):
        if type(self.inputData) is np.ndarray:
            return np.ndarray
        else:
            return self._check_input_file_type()

    def _check_input_file_type(self):
        self.possible_extensions = {
            "image": ["tif", "tiff", "png", "bmp", "gif", "jpg", "jpeg"],
            "npy": ["npy"],
        }
        extension = os.path.splitext(self.inputData)[1][1:].lower()
        for k in self.possible_extensions.keys():
            for e in self.possible_extensions[k]:
                if extension == e:
                    return k
        all_extensions = [
            x for x_list in self.possible_extensions.values() for x in x_list
        ]
        all_extensions += [x.upper() for x in all_extensions]
        raise ValueError(
            "Incorrect extension: {}. Possible values: {}.".format(
                extension, ", ".join(all_extensions)
            )
        )


def testblObjectSliced(
    N,
    trFilePath,
    orientation="x",
    resolution=[10.0e-9, 10.0e-9, 10.0e-9],
    thickness=131.0e-9,  # thickness of ech slice (m) 657.85e-9 , 492.6e-9, 985.2e-9, 448.93e-9
    delta=0.0031406,  # delta value of material in slice
    attenLength=0.277984e-6,  # attenuation length in material in slice
    trPreFilePath=None,  #'extensions/in/cs.tif'
    trPreProperties=None,  # expect [thickness delta, attenLength]
):
    """
     Function testing/demonstrating the generation of a list of N slices
     in structure expected for multisliceOptE

     orientation: direction 'x' or 'y' along which slices are arranged. imgF
                  is rotated 90 degrees if orientation='y'

     In this example the slices are offset from eachother to represent rotation
of the object the x or y axis.

     !! Note that delta and attenLength must match the energy used for propagation!  (this should be
                                                                                      managed better)
     
     trFilePathL: half left side of the image
     trFilePathR: half right side of the image
    """

    rx, ry, rz = resolution
    material = [thickness, delta, attenLength]

    # For this demonstraion we will use many slices of the same object with an offset
    # appropriate for a LTZP at angle ang.
    zoneRowWidth = 400e-9  # 420.0e-9 # 210.0e-9 #392.e-9
    ang = 27.0  # 27 # degrees
    # zStep = zoneRowWidth*math.cos(math.pi*ang/180)  #349.27 nm at 27D   #372.81 nm at 18D
    # lateralStep = zoneRowWidth*math.cos(math.pi*(90-ang)/180)  #177.96 nm at 27 D    #121.13 nm at 18 D
    zStep = 356.4e-9  # 492.6e-9# zoneRowWidth*math.tan(math.pi*(90-ang)/180)
    lateralStep = 181.6e-9  # zoneRowWidth *math.cos(math.pi*ang/180)
    # print (lateralStep)

    pp = propagationParameters(SemiAnalyt=1.0)

    slices = []

    # add a central stop
    if trPreFilePath is not None:
        trP = trPreProperties  # [1e-6, 0.00451796176,5.568403E-8]  # 1 um W: delta = 0.00451796176, att. len = 5.568403E-8 for 600 eV
        slices.append([trPreFilePath, [0.0, 0.0, 0.0, 0.0], resolution, trP, None, pp])

    for i in range(N):
        zi = i * zStep  #

        if orientation == "x":
            yi = 0.0
            xi = i * lateralStep - (lateralStep * N / 2.0)
            r = 0.0
            pad = None  # (int(((N+1)*(lateralStep)/rx))+500,int(((N+1)*(lateralStep)/rx))+500,500,500)

        if orientation == "y":
            xi = 0.0
            yi = i * lateralStep - (lateralStep * N / 2.0)
            r = 90.0
            pad = None  # (500, 500, int(((N+1)*(lateralStep)/ry))+500, int(((N+1)*(lateralStep)/ry))+500)

        trans = [xi, yi, zi, r]
        slices.append([trFilePath, trans, resolution, material, pad, pp])
    #
    #        if i% 2 == 0:
    #            slices.append( [trFilePathL,trans,resolution,material,pad,pp] )
    #        else:
    #            slices.append( [trFilePathR,trans,resolution,material,pad,pp] )
    #
    return slices


class multisliceOptE:
    """
    Multislice object definition and propagation.
    This started as an optical element but evolved into something like a variant
    of the beeamline class.  Might be better to have extended beamline class so that its functionality could be inherited.
    Slices are generated from images immediately before propagation through each slice
    to reduce memory overheads and/or file i/o associated with pre-creating
    many transmission function objects.



    objectSlices: list of  images and relative positions x, y and z (units !?) and rotation r (degrees):
                  [image File Path, trans=[x,y,z,r],resolution = [rx,ry,rz] (m/pixel),material[thickness, delta, attenuationLength]]

    """

    def __init__(self, objectSlices=None, description=None):
        # cutoff_background_noise=0.0,
        # background_color=0,
        # show=True)

        self.description = description

        NoneType = type(None)
        if isinstance(objectSlices, types.FunctionType):
            self.slices = objectSlices()
        elif isinstance(objectSlices, NoneType):
            self.slices = self.generateSlicesFromImage()
            self.numberSlices = len(self.slices)
        else:  # let's hope that it contains slices in expected structure!
            self.slices = objectSlices
            self.numberSlices = len(objectSlices)

        self.sliceIndex = 0
        self.sumTrans = None
        self.sumThickness = None

    def getAmplitudeTransmission(self):
        # return amplitude transmission characteristic of accumulative transmission function
        val = self.sumTrans.get_data(_typ=1)

        mesh = self.sumTrans.mesh
        nx, ny = mesh.nx, mesh.ny
        return np.array(val.reshape((ny, nx)))

    def getAmplitudeTransmission(self):
        # return intensity transmission characteristic of accumulative transmission function
        val = np.asarray(self.sumTrans.get_data(_typ=2))

        mesh = self.sumTrans.mesh
        nx, ny = mesh.nx, mesh.ny
        return np.array(val.reshape((ny, nx)))

    def getOpticalPathLength(self):
        # return Optical Path Difference transmission characteristic of accumulative transmission function
        val = np.asarray(self.sumTrans.get_data(_typ=3))

        mesh = self.sumTrans.mesh
        nx, ny = mesh.nx, mesh.ny
        return np.array(val.reshape((ny, nx)))

    def getSRWTransmissionFunction(self):
        # return accumulative transmission function as repersented in SRWLOptT, see srwlib.h
        # return self.sumTrans.arTr
        return self.sumTrans

    #    def getThickness(self,attenuationLength=None, delta = None):
    #        #return thickness derived from accumulative transmission function, under assumption
    #        # of uniform material characterised by attenuationLength
    #        # INPUT:  specify attenuationLength OR delta
    #
    #        thickness = None
    #        if attenuationLength is not None:
    #            thickness =  2*self.getAttenuationLength()*math.log(self.getAmplitudeTransmission())
    #
    #        if delta is not None:
    #            thickness = -1*self.getOpticalPathDifference()/delta
    #
    #        return thickness

    def getThickness(self):
        mesh = self.sumTrans.mesh
        nx, ny = mesh.nx, mesh.ny
        thickness = np.array(self.sumThickness.reshape((ny, nx)))
        return thickness

    def getDelta(self):
        # return equivalent delta value for accumulative transmission function

        val = -1 * np.divide(self.getOpticalPathLength(), self.getThickness())

        mesh = self.sumTrans.mesh
        nx, ny = mesh.nx, mesh.ny
        return np.array(val.reshape((ny, nx)))

    def getAttenuationLength(self):
        """
        return equivalent AttenuationLength for accumulative transmission function at the wavelength
        for which the transmission function was constructed.
        """
        val = -0.5 * np.divide(
            self.getThickness(), math.log(getAmplitudeTransmission())
        )

        mesh = self.sumTrans.mesh
        nx, ny = mesh.nx, mesh.ny
        return np.array(val.reshape((ny, nx)))

    def getPhaseShift(self, wavelength):
        """
        return phase shift through accumulative transmission function 
        phase shift is (2*pi*delta / wavelength) * pathInBody
        """
        val = (
            2 * math.pi * np.multiply(self.getDelta(), self.getThickness()) / wavelength
        )

        # if val > 2*math.pi:       #Alaleh Added
        #     val = val - 2*math.pi

        mesh = self.sumTrans.mesh
        nx, ny = mesh.nx, mesh.ny
        return np.array(val.reshape((ny, nx)))

    def writeTransmissionFunction(self, filename, comments=""):

        print("Saving transmission function to " + filename)
        # mkdir_p(os.path.dirname(incremental_filename))
        with open(filename, "wb") as f:
            pickle.dump(self.sumTrans, f)

    def addToSumTrans(self, tr):

        # this is old and wrong!!  don't us it

        import numpy as np

        if self.sumTrans is None:

            mesh = tr.mesh
            nx = mesh.nx
            ny = mesh.ny

            print("creating container with dimensions %d,%d" % (nx, ny))

            self.sumTrans = SRWLOptT(_nx=nx, _ny=ny)

        if len(tr.arTr) == len(self.sumTrans.arTr):
            self.sumTrans.arTr = np.add(self.sumTrans.arTr, tr.arTr)
        else:
            print(
                "Warning - refusing to add transmission funnctions of different dimensions"
            )

        self.addToSumThickness(tr)

    def addToSumThickness(self, tr):
        if self.sumThickness is None:
            self.sumThickness = np.asarray(tr.pathInBody)
        else:
            if len(tr.pathInBody) == len(self.sumThickness):
                self.sumThickness = np.add(self.sumThickness, np.asarray(tr.pathInBody))
            else:
                print(
                    "Warning - cowardly refusing to add thickness of slice with different number of elements"
                )

        # for testing only...
        # mesh = tr.mesh
        # nx = mesh.nx
        # ny = mesh.ny

        # import matplotlib.pyplot as plt
        # fig,x  = plt.subplots(nrows=1)
        # th = plt.imshow(np.array(self.sumThickness.reshape((ny, nx))))
        # fig.colorbar(th)
        # plt.show()
        print("Thickness (max, accumulative): %f um" % np.max(self.sumThickness * 1e6))

    def showTransmissionFunction(self):
        show_transmissionFunctionCplx(self.getTransmissionFunction())

    def print_beamline(self):
        print(self.description)  # Sufficient for compatibility, can do better

    def generateSlicesFromImage(
        self,
        imgF,  # one image or a list of N images
        N,
        translation,  # one or list of N x,y shifts for images
        resolution,  # one or a list of N resolution values
        material,
        pad,
    ):
        # =============================================================================
        #
        #     Uses image, or list of images specified in imgF to build a multislice object
        #
        #     imgF: image file name  or list of file names
        #       N:  number of slices
        #        resolution:  =[rx,ry,rz)], resolution in x, y and z direction
        #                     where rz is the slice separation
        #       material = [thickness= 657.85e-9, # not necessarily same a rz
        #                   delta= 0.0031406,
        #                     attenLength= 0.277984e-6]

        #
        # =============================================================================

        print("%d slices of %s coming up..." % (N, imgF))

        if not isinstance(imgF, list):
            #    assert (len(imgF) == N and len(translation)==N and len(resolution)==N,
            #            'The number of values in  imgF, translation, resolution  and materialmust be the same')
            imgF = [imgF for i in range(N)]

        if not isinstance(translation, list):
            translation = [translation for i in range(N)]

        if not isinstance(resolution, list):
            translation = [resolution for i in range(N)]

        if not isinstance(material[0], list):
            material = [material for i in range(N)]

        if not isinstance(pad[0], list):
            pad = [pad for i in range(N)]

        slices = []
        for i in range(N):

            slices.append(
                [
                    imgF[i],
                    translation[i],
                    resolution[i],
                    material[i],
                    pad[i],
                    propagationParameters(),
                ]
            )

        self.slices = slices

        return slices

    #    def getNext(self):
    #
    #        self.sliceIndex +=1
    #        return self.getSlice(self.sliceIndex)

    def getSliceParam(self, n):
        # return the nth slice parameters
        return self.slices[n]

    def printSliceParams(self):
        # pretty print  slice parameters
        table = BeautifulTable()
        table.column_headers = [
            "Image Path",
            "Translation",
            "Resolution",
            "Material",
            "Pad",
            "PropagationParameters",
        ]
        for s in self.slices:
            table.append_row(s)
        print(table)

    def loadImageFromFile(self):
        import_err_msg = '"{}" library cannot be imported. Please install it first with "pip install {}".'
        try:
            from PIL import Image
        except ImportError:
            raise ImportError(import_err_msg.format("PIL", "pillow"))

        d = Image.open(self.inputData)
        self.trdata = d["data"]

    def getSlice(
        self,
        n,
        show=False,
        outFilePath=None,
        saveProcessedImage=False,
        driftSpaces=True,
    ):

        # return the nth slice as a beamline after creating the transmission function

        #          # we assume that values of z in trans are absolute, so we need the difference
        #          if n>0:
        #              imgPath0,trans0,resolution0,material0,pad0,propagationParams0 = self.getSliceParam(n-1)
        #          else:
        #              z0 = 0
        #              imgPath0=''
        #
        (
            imgPath,
            trans,
            resolution,
            material,
            pad,
            propagationParams,
        ) = self.getSliceParam(n)

        x, y, z, r = trans
        zstep = z
        rx, ry, rz = resolution
        thickness, delta, attenLength = material

        #   propagationDistance = z-z0

        el, pp = [], []

        # For some situations it might be appropriate to add a drift space between slices.
        # This requres careful consideration...
        if driftSpaces == True:
            if z != 0:
                # insert drift space to element
                # el.append(SRWLOptD(propagationDistanceF))
                el.append(SRWLOptD(zstep))  # for LTZP
                pp.append(propagationParameters(SemiAnalyt=1.0))

        # if imgPath != imgPath0:
        #    self.loadImageFromFile(imgPath)

        # create optical element transmission function
        tr = srwl_opt_setup_transm_from_file(
            inputData=imgPath,  # self.trdata,
            resolution=rx,  # None that ry is ignored here!
            thickness=thickness,  # this is max thickness!
            delta=delta,
            atten_len=attenLength,
            xc=0.0,
            yc=0.0,
            area=None,
            rotate_angle=r,
            rotate_reshape=False,
            cutoff_background_noise=None,
            background_color=255,
            tile=None,
            shift_x=int(pixelScale(x, rx)),
            shift_y=int(pixelScale(y, ry)),
            invert=True,
            pad=pad,
        )
        # is_save_images=saveProcessedImage,
        # prefix=imgPath[:-4]+'_processed_' + str(n),
        # output_image_format='tif',
        # )

        # keep accumulative sum of slice transmission functions,
        self.addToSumTrans(tr)

        el.append(tr)
        pp.append(
            propagationParameters()
        )  # use simplest assumption for drift space propagation paramters.  To be improved!

        if show:
            # show_transmission(tr)
            show_transmissionFunctionCplx(get_transmissionFunctionCplx(tr))

        if outFilePath:
            incremental_filename = os.path.join(
                outFilePath, "%04d_%s.tr" % (n, "slice")
            )
            print("Saving transmission function to " + incremental_filename)
            # mkdir_p(os.path.dirname(incremental_filename))
            with open(incremental_filename, "wb") as f:
                pickle.dump(tr, f)

        return blFromOE(el, pp, description="Slice %d" % self.sliceIndex)

    # def propagate(self,wfr):  # wrapper for compatibility
    #    self.propagateSlices(wfr=wfr)

    def propagateSlices(
        self, wfr, outFilePath=None, showSliceSummary=False, plot=False
    ):
        """
        Propagate wavefront through multislice object (beamline).

        :param wf: Input wavefront (rewritten after propagation)
        :type wf: wpg.wavefront.Wavefront
        """
        self.sumThickness = None  # re-initialise the thickness map
        for i in range(self.numberSlices):

            bl = self.getSlice(i)

            # print_beamline(bl)

            # Switch to frequency domain.
            # srwl.SetRepresElecField(wfr._srwl_wf, 'f')

            # Save spectrum for later reference.
            # sz0 = get_intensity_on_axis(wfr);
            # wf.custom_fields['/misc/spectrum0'] = sz0

            # Propagate.
            bl.propagate(wfr)  # (exit wavefield saved in wf)

            # Save spectrum after propagation for later reference.
            # sz1 = get_intensity_on_axis(wfr);
            # wf.custom_fields['/misc/spectrum1'] = sz1

            # Switch back to time domain.
            # srwl.SetRepresElecField(wfr._srwl_wf, 't')

            if outFilePath is not None:
                incremental_filename = os.path.join(outFilePath, "%04d.h5" % (i))
                print("Saving propagated wavefront to " + incremental_filename)
                # mkdir_p(os.path.dirname(incremental_filename))
                wfr.store_hdf5(incremental_filename)

            if showSliceSummary == True:
                I = wfr.get_intensity().sum()  # calculated sum of the intensity
                print("Total intensity = %f" % I)
                # print('FWHMx[um]:', calculate_fwhm_x(wfr) * 1e6)
                # print('FWHMy [um]:', calculate_fwhm_y(wfr) * 1e6)

            if plot == True:
                # display the wavefield
                plotWavefront(wfr, "Wavefield")

            print("Done with propagation step %d." % (i))
            print("#" * 80)

    def propagate(self, wfr, outFilePath=None, plot=False):
        """
        Propagate wavefront through multislice object (beamline).

        :param wf: Input wavefront (rewritten after propagation)
        :type wf: wpg.wavefront.Wavefront
        """
        self.propagateSlices(wfr, outFilePath=outFilePath, showSliceSummary=False)


def loadFromCache(filename):
    # filename is a pickled SRWLOptT object
    with open(filename, "rb") as pickle_file:
        trOpt = pickle.load(pickle_file)

    return trOpt


def testgreyscaleToSlices(N):

    wf = constructTestWaveField(z=5)

    msobject = greyscaleToSlices(
        "/opt/wpg/wpg/extensions/in/ic.tif",
        slices=N,
        invert=False,
        resolution=[600.0e-9, 600.0e-9, None],
        thickness=657.85e-9,  # TOTAL thickness (m)
        delta=0.0031406,  # delta value of material in slices
        attenLength=0.277984e-6,
    )

    # msobject.tr.showTransmissionFunction()

    msobject.tr.propagate(wf)

    plotWavefront(
        wf, "Wavefield propagated through multislice object",
    )
    writeIntensity(
        wf, "AfterObject", path="./out", polarization="horizontal", imgType="tif"
    )


def testgreyscaleToSlicesForLTZP():

    import matplotlib.pyplot as plt
    import imageio

    N = 3
    outPath = "/opt/wpg/wpg/extensions/out/50Slices63D/"
    # outPath='/home/alaleh/Desktop/out/50Slices63D/'
    sourceWavefield = "/home/alaleh/Desktop/out/2019-09-12/50Slices63D/wfPreOptic.h5"
    thicknessMap = "/opt/wpg/wpg/extensions/in/3nmRes63D_rz_30um.tif"
    wavelength = 2.066e-9  # m

    wf = constructTestWaveField(
        z1=15
    )  # alaleh's choice of z1, but moved to argument here, change to function undone)

    # wf = Wavefront()
    # wf.load_hdf5( sourceWavefield )

    plotWavefront(wf, "Wavefield before optic")
    I0 = wf.get_intensity(polarization="horizontal")

    # from extensions.driftVol import greyscaleToSlices
    msobject = greyscaleToSlices(
        thicknessMap,
        Use_PP(
            semi_analytical_treatment=1.0, zoom=0.3, sampling=1
        ),  # should reconsider zoom and sampling here
        slices=N,
        invert=False,
        resolution=[3.0e-9, 3.0e-9, None],
        thickness=655.91e-9,  # TOTAL thickness (m)
        delta=0.0031406,  # delta value of material in slices
        attenLength=0.277984e-6,
    )

    # msobject.tr.showTransmissionFunction()

    msobject.tr.propagate(wf, outFilePath=outPath)

    # compare I0/I1 assuming no zoom or scaling
    I1 = wf.get_intensity(polarization="horizontal")
    I1_total = I1.sum()
    I0_total = I0[I1 != 0].sum()
    transmission = I1_total / I0_total
    absorption = 1 - transmission
    print(
        "I0={}, I1={}, Absorption={}, Transmission={}".format(
            I0_total, I1_total, absorption, transmission
        )
    )

    plotWavefront(
        wf, "Wavefield propagated through multislice object",
    )
    writeIntensity(
        wf, "AfterObject", path="./out", polarization="horizontal", imgType="tif"
    )

    transmission = msobject.tr.getAmplitudeTransmission()
    thickness = msobject.tr.getThickness()
    phaseShift = msobject.tr.getPhaseShift(wavelength)
    plot(transmission, thickness, phaseShift)

    plotThickness(tr)

    imageio.imwrite(outPath + "/phaseShift.tif", np.float32(phaseShift))
    imageio.imwrite(outPath + "/thickness.tif", np.float32(thickness))
    imageio.imwrite(outPath + "/transmission.tif", np.float32(transmission))


def test_multislice_from_images():
    slices = testblObjectSliced(5, trFilePath="./in/ic.tif")
    multiSliceLTZP_X = multisliceOptE(
        objectSlices=slices, description="LTZP (multislice)"
    )
    multiSliceLTZP_X.printSliceParams()
    multiSliceLTZP_X.description = "LTZP X"
    multiSliceLTZP_X.propagate(wf)
    cplx = multiSliceLTZP_X.getTransmissionFunction()
    phase = multiSliceLTZP_X.getTransmissionFunction(part="phase")
    amp = multiSliceLTZP_X.getTransmissionFunction(part="amplitude")

    multiSliceLTZP_X.showTransmissionFunction()

    # now save tr
    multiSliceLTZP_X.writeTransmissionFunction(
        "out.pickle", comments="this is my transmission function, picled SRWLOPTT type"
    )

    # tr=loadFromCache('out.pickle')
    # cplx=get_transmissionFunctionCplx(tr)
    # show_transmissionFunctionCplx(cplx)


def plotThickness(tr):
    from extensions.utilPlot import plotTrThickness

    plotTrThickness(tr)


def plot(transmission, thickness, phaseshift):
    cmap = plt.get_cmap("PiYG")

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)

    th = ax0.imshow(thickness, cmap=cmap)
    fig.colorbar(th, ax=ax0)
    ax0.set_title("Thickness")

    ps = ax1.imshow(phaseshift, cmap=cmap)
    fig.colorbar(ps, ax=ax1)
    ax1.set_title("Phase Shift")

    tm = ax2.imshow(transmission, cmap=cmap)
    fig.colorbar(tm, ax=ax2)
    ax2.set_title("Transmission")

    fig.tight_layout()
    plt.show()


def constructTestWaveField(npoints=256, ekeV=0.6, z1=5, show=False):
    # For testing ....
    # Generate a simple Gaussian wavefield for testing purposes

    from wpg.srwlib import srwl_uti_ph_en_conv
    from extensions.gvrutils import calculate_theta_fwhm_cdr
    from wpg.generators import build_gauss_wavefront_xy

    # set wavefield parameters
    qnC = 0.01
    wlambda = srwl_uti_ph_en_conv(ekeV, _in_u="keV", _out_u="nm")
    theta_fwhm = calculate_theta_fwhm_cdr(ekeV, qnC)
    k = 2 * np.sqrt(2 * np.log(2))
    range_xy = (theta_fwhm / k * z1 * 5.0) * 2.0
    sigX = 12.4e-10 * k / (ekeV * 4 * np.pi * theta_fwhm)

    # construct wavefield
    wf0 = build_gauss_wavefront_xy(
        nx=npoints,
        ny=npoints,
        ekev=ekeV,
        xMin=-range_xy / 2,
        xMax=range_xy / 2,
        yMin=-range_xy / 2,
        yMax=range_xy / 2,
        sigX=sigX,
        sigY=sigX,
        d2waist=z1,
        _mx=0,
        _my=0,
    )
    wfr = Wavefront(srwl_wavefront=wf0)
    print(
        "Wavelength=%f, theta FWWM=%f, range XY = %f, sig X = %f"
        % (wlambda, theta_fwhm, range_xy, sigX)
    )

    wfr._srwl_wf.unitElFld = 1  #'sqrt(Phot/s/0.1%bw/mm^2)'

    if show == True:  # display the wavefield
        plotWavefront(wfr, "Wavefield at source")

    return wfr


if __name__ == "__main__":

    # test_multislice_from_images()

    # testgreyscaleToSlices(3)

    """
    testblObjectSliced(N,trFilePath,orientation='x',
                       resolution = [10.0e-9,10.0e-9,10.0e-9],
                       thickness= 131.0e-9,       #thickness of ech slice (m) 657.85e-9 , 492.6e-9, 985.2e-9, 448.93e-9
                       delta= 0.0031406,           #delta value of material in slice
                       attenLength= 0.277984e-6,  #attenuation length in material in slice
                       trPreFilePath = None #'extensions/in/cs.tif'
                       )
    """
    pass
