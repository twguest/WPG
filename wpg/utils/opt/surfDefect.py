import numpy as np

from skimage import draw
from copy import copy
from matplotlib import pyplot as plt
from scipy import ndimage as nd


class Defect:
    def __init__(self, nx, ny):

        self.nx, self.ny = nx, ny
        self.mask = np.ones((self.nx, self.ny), dtype=np.uint8)

    def slit(self, dx, dy, _x=0):
        """
        slit dimensions dx,dy,_x are in m
        """

        maxi = int(len(self.mask))

        dx = int(dx)
        dy = int(dy)
        xshift = int(dx / 2)
        yshift = int(dy / 2)

        _x = int(_x)

        xc = int(maxi / 2)
        yc = int(len(self.mask) - (maxi / 2))

        start = (yc - yshift, _x + xc - xshift)
        end = (yc + yshift, _x + xc + xshift)

        r, c = draw.rectangle(start=start, end=end, shape=self.mask.shape)
        self.mask[r, c] = 0
