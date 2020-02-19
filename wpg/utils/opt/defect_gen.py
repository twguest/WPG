import numpy as np

from skimage import draw
from copy import copy
from matplotlib import pyplot as plt
from scipy import ndimage as nd
from scipy.ndimage import gaussian_filter as gsnblur
from numpy.random import randint

# from extensions.utils.mask_utils import pixels2rspace as p2r
import pandas as pd
from sklearn.preprocessing import minmax_scale as norm


class Defect:
    def __init__(self, nx, ny):
        self.nx = nx
        self.ny = ny

        self.base = np.zeros((nx, ny))
        self.surface = copy(self.base)

        self.px, self.py = None, None

        self.info = pd.DataFrame()
        self.locs = []
        self.heights = []
        self.rads = []

    def add_rspace(self, px, py):
        self.px = px
        self.py = py

    def bump_height(self, r):
        omega = r / (2.355 * 1.5)
        h = (16 * np.log(2) / 18) * r
        return h

    def add_defect_gsn(self, r=1, sevd=1, loc="random", rmax=15):
        """
        """

        layer = copy(self.base)

        if loc == "random":
            loc1 = randint(self.nx)
            loc2 = randint(self.ny)
            layer[loc1, loc2] = 1
            layer = gsnblur(layer, r / 2)
            height = norm(
                [0, self.bump_height(r / 2), self.bump_height(rmax / 2)], (0, 255)
            )[1]
            layer *= height
        elif type(loc) == tuple:
            loc1 = loc[0]
            loc2 = loc[1]
            layer[loc1, loc2] = 1
            layer = gsnblur(layer, r / 2)
            height = norm(
                [0, self.bump_height(r / 2), self.bump_height(rmax / 2)], (0, 255)
            )[1]
            layer *= height
        else:
            print("Unrecognised Location Sequence")

        self.surface += layer

        if self.px is not None:
            self.locs.append((loc1 * self.px, loc2 * self.py))
            self.heights.append((height / 255) * self.bump_height(rmax / 2) * self.px)
            self.rads.append(r * self.px)
            self.h0 = self.bump_height(rmax / 2) * self.px

        else:
            self.locs.append((loc1, loc2))
            self.heights.append(height)
            self.rads.append(r)

    def set_df(self):

        if self.px is not None:
            self.info["Location (m)"] = self.locs
            self.info["Radius (m)"] = self.rads
            self.info["Height (m)"] = self.heights
        else:
            self.info["Location (pix)"] = self.locs
            self.info["Radius (pix)"] = self.rads
            self.info["Height (pix)"] = self.heights

    def plot(self):
        plt.imshow(self.surface)

    def export(self, outfile):
        np.save(outfile, self.surface.astype(int))

        self.set_df()
        self.info.to_pickle(outfile + ".df")

        print("Maximum Defect Height: {}".format(self.h0))


def generate_surface(wfr):

    defect = Defect(wfr.params.Mesh.nx, wfr.params.mesh.ny)
    px, py = wfr.pixel_size()
    defect.add_rspace(px, py)

    return defect


def random_defect_distribution(defect, n, r=(5, 15)):

    itr = 0
    while itr in range(n):
        defect.add_defect_gsn(np.random.random() * np.random.randint(r[0], r[1]))
        itr += 1


if __name__ == "__main__":
    print("Testing Defect Generator")
    defects = Defect(250, 250)
    defects.add_rspace(1e-09, 1e-09)
    random_defect_distribution(defects, 20, r=(5, 60))
    defects.plot()
    defects.set_df()
    print(defects.info)
