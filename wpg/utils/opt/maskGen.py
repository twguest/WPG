# from extensions.multisliceOptE import *
from extensions.utils.mask_utils import pixels2rspace as p2r
from extensions.twpg_wavefront import Wavefront
from extensions.utils.mask_utils import Photomask

import numpy as np

from os import listdir


class maskGen:
    """
    generates a double slit mask given input spec. and wfr properties
    """

    def __init__(self, wfr):
        self.wfr = wfr
        self.px, self.py = self.wfr.pixelsize()
        self.nx, self.ny = self.wfr.params.Mesh.nx, self.wfr.params.Mesh.ny

        print("Mask Pixels Dimensions: {} x {}".format(self.nx, self.ny))
        print("Mask Pixels Dimensions (R-Space): {} x {}".format(self.px, self.py))

    def generate(self, x, y, sep, n=1, invert=True):

        print("Generating Mask")
        print("Slit Width: {}".format(x))
        print("Slit Height: {}".format(y))
        print("Slit Seperation: {}".format(sep))

        self.mask = Photomask(self.wfr.params.Mesh.nx, self.wfr.params.Mesh.ny)
        self.mask.nslits(p2r(self.px, x), p2r(self.py, y), p2r(self.px, sep), n=n)

        if invert == True:
            print("Inverting Mask")
            self.mask.mask = (self.mask.mask - 1) ** 2

    def save(self, outdir):
        print("Saving Mask as np Array")
        np.save(outdir, self.mask.mask)
        print(type(self.mask.mask))
        print("Mask Saved @: {}".format(outdir))


if __name__ == "__main__":
    print("Testing Operations")
    maskdir = r"./extensions/in/masks/"
    maskname = r"testMask.npy"

    if maskname in listdir(maskdir):
        print("Loading Mask from File")
        mask = np.load(maskdir + maskname)
    else:
        wfr = Wavefront()
        wfr.load_hdf5(
            r"/nfs/data/users/twg/single_electron/15_SRWLOptD/wfr_mode_se.hdf5"
        )

        mask = maskGen(wfr)
        mask.generate(100e-09, 500e-09, 1e-06)
        mask.save(maskdir + maskname)
