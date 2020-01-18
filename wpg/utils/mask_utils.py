import numpy as np

from skimage import draw
from copy import copy
from matplotlib import pyplot as plt
from scipy import ndimage as nd

class Photomask():
    def __init__(self, nx, ny):
        
        self.nx, self.ny = nx, ny
        self.mask = np.ones((self.nx, self.ny), dtype = np.uint8)    
            
        
    def slit(self, dx, dy, _x = 0):
        """
        slit dimensions dx,dy,_x are in m
        """

        maxi = int(len(self.mask))
        
        dx = int(dx)
        dy = int(dy)
        xshift = int(dx/2)  
        yshift = int(dy/2) 

        _x = int(_x)
        
        xc = int(maxi/2)
        yc = int(len(self.mask)-(maxi/2))

        start = (yc - yshift, _x + xc - xshift)
        end = (yc + yshift, _x + xc + xshift)
        
        r,c = draw.rectangle(start = start, end = end, shape = self.mask.shape)
        self.mask[r,c] = 0
        
        
    def nslits(self, dx, dy, _a = 100, n = 2):
    
        for N in range(n):
            _x = (dx*N*2)+_a//2
            self.slit(dx,dy,_x = _x)
            self.slit(dx,dy,_x = -_x)
            
    def save_array(self, filename):
        np.save(filename, self.mask)
        
    def fourbeam(self, dx, dy, _a = None, n = 1):
        
        if _a == None:
            _a == dx
        
        for N in range(n):
            _x = (dx*N*2)+_a//2
            self.slit(dx, dy, _x)
            
        self.img = 0
        
        self.img += np.rot90(self.mask)
        self.img += np.rot90(np.rot90(self.mask))
        self.img += np.rot90(np.rot90(np.rot90(self.mask)))
        self.img += np.rot90(np.rot90(np.rot90(np.rot90(self.mask))))

        self.mask = copy(self.img)            
        
    
    def fivebeam(self, dx, dy, _a = None, n = 1):
        
        if _a == None:
            _a == dx
        
        for N in range(n):
            _x = (dx*N*2)+_a//2
            self.slit(dx, dy, _x)
            
        self.img = 0
        
        for ang in range(5):
            self.img += nd.rotate(self.mask, angle = ang*360/5+360/4, reshape = False)
            
        self.mask = copy(self.img)            
        
        
def pixels2rspace(pixelsize, actualsize):
    npixels = actualsize/pixelsize
    return npixels



if __name__ == "__main__":
    
    Mask = Photomask(5000,5000)
    Mask.fivebeam(100,500,1000,10)
    plt.imshow(Mask.mask)
    
    Mask = Photomask(500,500)
    Mask.fourbeam(10,50,100,10)
 
    
            
