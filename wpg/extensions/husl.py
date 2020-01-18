import operator
import math
import numpy as np
import warnings

__version__ = "4.0.3"

# original scalar functions of husl-python are available with the
# the '_sc_' prefix. 

m = [
    [3.240969941904521, -1.537383177570093, -0.498610760293],
    [-0.96924363628087, 1.87596750150772, 0.041555057407175],
    [0.055630079696993, -0.20397695888897, 1.056971514242878],
]
m_inv = [
    [0.41239079926595, 0.35758433938387, 0.18048078840183],
    [0.21263900587151, 0.71516867876775, 0.072192315360733],
    [0.019330818715591, 0.11919477979462, 0.95053215224966],
]
refX = 0.95045592705167
refY = 1.0
refZ = 1.089057750759878
refU = 0.19783000664283
refV = 0.46831999493879
kappa = 903.2962962
epsilon = 0.0088564516

# a decorator to distinguish and handle triple-scalar or triple array input
def triplescalar(func):
    def inner(triple):
        # check for scalar
        was_all_scalar = np.all([np.isscalar(t) or np.asarray(t).ndim==0 for t in triple])
        
        if was_all_scalar:
            # all array
            triple = [np.array(t,ndmin=1) for t in triple]
            # all array single entry
            triple = [np.array(t[0],ndmin=1) for t in triple]
        else:
            # enlarge dimensions
            max_dims = np.max([np.asarray(t).ndim for t in triple])
            triple = [np.array(t,ndmin=max_dims)for t in triple]
            # enlarge shape
            max_shape =np.array([np.asarray(t).shape for t in triple]).max(0)
            triple = [np.resize(t,max_shape) for t in triple]
            
        ret = func(triple)
        if was_all_scalar:
            return [np.float(r[0]) for r in ret]
        else:
            return ret
    return inner
    
# Public API
def complex_to_rgb(z=None,amin=None,amax=None,mode='special',phstart=0.,sat=1.0,as_image=True):
    """\
    Complex to RGB (Red,Green,Blue) transformation.
    Attempts to image a complex array *z* in R, G, B coloration, where
    phase(z) is mapped to Hue and podulus(z) is mapped to a choice of
    perceived brightness.
    
    Parameters
    ----------
    z : array-like
        Input should be two-dimensional. If no input is given, a colorwheel
        is produced.
        
    amin : float
        All z-values with modulus below `amin` will appear black. Defaults
        to np.abs(z).min() if no value for `amin` was provided
    
    amix : float
        All z-values with modulus below `amax` will appear white. Defaults
        to 0.9*np.abs(z).max() if no value for `amax` was provided
        
    mode : str
        Choose between 'special','chroma' or 'pastell'. For 'special', 
        the individual channels will not receive gamma correction but 
        the luminance (Y-value) with the effect that the image will
        preserve the modulus when watched in grayscale. For 'chroma', 
        the chroma of the channels is distorted while preserving the 
        lightness, with the effect, that blue and red channel become 
        saturated. For 'pastell', lightness and chroma are preserved 
        by changing saturation to the appropriate value, with the effect, 
        that only pastell colors are available.
    
    phstart : float
        Starting (hue) of for the phase. Choose in the range [0-2*pi].
        
    sat : float
        Saturation value, defaults to 1.0 (maximum saturation when possible).
        
    as_image : bool
        See return
        
    Returns
    -------
    rgb : ndarray
        If `as-images` is True, returns 8bit array (m,n,3) where last 
        axis is RGB color.
        If `as-images` is False, returns all three channels concatenated,
        RGB (3,m,n), as float values in the range [0,1.].
        
    """
    if z is None:
        x,y=np.indices((200,200))-99.5
        z=x+1j*y
        amax = 100
        amin = 0
        
    H=np.degrees(np.angle(z) + np.pi + phstart) % 360. 
    A = np.abs(z)
    amin = A.min() if amin is None else amin
    amax = 1.2*A.max() if amax is None else amax
    if np.allclose(amin,amax):
        amin = 0
    if np.allclose(0,amax):
        amix = 1
        
    A=(A-amin)/(amax-amin)
    S= sat * 100
    if str(mode)=='special':
        XYZ = np.asarray(luv_to_xyz(lch_to_luv(huslp_to_lch([H,S,np.sqrt(A)*100]))))
        R,G,B = [np.sum(XYZ*np.array(mi).reshape((3,)+(XYZ.ndim-1)*(1,)),0) for mi in m]
    elif str(mode)=='chroma':
        R,G,B = husl_to_rgb(H,S, A*100 )
    else:
        R,G,B = huslp_to_rgb(H,S, A*100 )
        
    if as_image:
        return np.uint8(np.array([R,G,B]).swapaxes(0,2) * 255)
    else:
        return np.array([R,G,B])
    
def husl_to_rgb(h, s, l):
    return lch_to_rgb(*husl_to_lch([h, s, l]))


def husl_to_hex(h, s, l):
    return rgb_to_hex(husl_to_rgb(h, s, l))


def rgb_to_husl(r, g, b):
    return lch_to_husl(rgb_to_lch(r, g, b))


def hex_to_husl(hex):
    return rgb_to_husl(*hex_to_rgb(hex))


def huslp_to_rgb(h, s, l):
    return lch_to_rgb(*huslp_to_lch([h, s, l]))


def huslp_to_hex(h, s, l):
    return rgb_to_hex(huslp_to_rgb(h, s, l))


def rgb_to_huslp(r, g, b):
    return lch_to_huslp(rgb_to_lch(r, g, b))


def hex_to_huslp(hex):
    return rgb_to_huslp(*hex_to_rgb(hex))


def lch_to_rgb(l, c, h):
    return xyz_to_rgb(luv_to_xyz(lch_to_luv([l, c, h])))


def rgb_to_lch(r, g, b):
    return luv_to_lch(xyz_to_luv(rgb_to_xyz([r, g, b])))


def _sc_get_bounds(L):
    sub1 = ((L + 16.0) ** 3.0) / 1560896.0
    sub2 = sub1 if sub1 > epsilon else L / kappa
    ret = []
    for [m1, m2, m3] in m:
        for t in [0, 1]:
            top1 = (284517.0 * m1 - 94839.0 * m3) * sub2
            top2 = (838422.0 * m3 + 769860.0 * m2 + 731718.0 * m1) * L * sub2 - 769860.0 * t * L
            bottom = (632260.0 * m3 - 126452.0 * m2) * sub2 + 126452.0 * t
            ret.append((top1 / bottom, top2 / bottom))
    return ret


def _sc_intersect_line_line(line1, line2):
    return (line1[1] - line2[1]) / (line2[0] - line1[0])


def _sc_distance_from_pole(point):
    return math.sqrt(point[0] ** 2 + point[1] ** 2)


def _sc_length_of_ray_until_intersect(theta, line):
    m1, b1 = line
    length = b1 / (math.sin(theta) - m1 * math.cos(theta))
    if length < 0:
        return None
    return length


def _sc_max_safe_chroma_for_L(L):
    lengths = []
    for [m1, b1] in _sc_get_bounds(L):
        x = _sc_intersect_line_line((m1, b1), (-1.0 / m1, 0.0))
        lengths.append(_sc_distance_from_pole((x, b1 + x * m1)))
    return min(lengths)


def _sc_max_chroma_for_LH(L, H):
    hrad = H / 360.0 * math.pi * 2.0
    lengths = []
    for line in _sc_get_bounds(L):
        l = _sc_length_of_ray_until_intersect(hrad, line)
        if l is not None:
            lengths.append(l)
    return min(lengths)

def _get_bounds(L):
    sub1 = ((L + 16.0) ** 3.0) / 1560896.0
    sub2 = np.array(L) / kappa
    sub2[sub1 > epsilon] = sub1[sub1 > epsilon] 
    ret = []
    for [m1, m2, m3] in m:
        for t in [0, 1]:
            top1 = (284517.0 * m1 - 94839.0 * m3) * sub2
            top2 = (838422.0 * m3 + 769860.0 * m2 + 731718.0 * m1) * L * sub2 - 769860.0 * t * L
            bottom = (632260.0 * m3 - 126452.0 * m2) * sub2 + 126452.0 * t
            ret.append((top1 / bottom, top2 / bottom))
    return ret


def _intersect_line_line(line1, line2):
    return (line1[1] - line2[1]) / (line2[0] - line1[0])

def _length_of_ray_until_intersect(theta, line):
    m1, b1 = line
    length = b1 / (np.sin(theta) - m1 * np.cos(theta))
    """
    if length < 0:
        return None
    """
    length[length<0]=np.infty
    return length


def _max_safe_chroma_for_L(L):
    lengths = []
    for [m1, b1] in _get_bounds(L):
        x = _intersect_line_line((m1, b1), (-1.0 / m1, 0.0))
        lengths.append(np.sqrt(x**2 + (b1 + x * m1)**2))
    return np.asarray(lengths).min(0)


def _max_chroma_for_LH(L, H):
    hrad = H / 360.0 * np.pi * 2.0
    lengths = []
    for line in _get_bounds(L):
        l = _length_of_ray_until_intersect(hrad, line)
        """
        if l is not None:
            lengths.append(l)
        """
        lengths.append(l)
    return np.asarray(lengths).min(0)


def _sc_dot_product(a, b):
    return sum(map(operator.mul, a, b))

def _f(t):
    ret = np.array(t) / refY * kappa 
    ret[t > epsilon] =116 * np.power((t[t > epsilon] / refY), 1.0 / 3.0) - 16.0
    return ret  
    

def _sc_f(t):
    if t > epsilon:
        return 116 * math.pow((t / refY), 1.0 / 3.0) - 16.0
    else:
        return (t / refY) * kappa


def _f_inv(t):
    ret = np.array(t) * refY / kappa
    ret[t>8]=  refY * np.power((t[t>8] + 16.0) / 116.0, 3.0)
    return ret


def _sc_f_inv(t):
    if t > 8:
        return refY * math.pow((t + 16.0) / 116.0, 3.0)
    else:
        return refY * t / kappa


def _from_linear(c):
    c=np.array(c)
    ret = 12.92 * np.array(c)
    ret[c > 0.0031308]= 1.055 * np.power(c[c > 0.0031308], 1.0 / 2.4) - 0.055
    return ret


def _sc_from_linear(c):
    if c <= 0.0031308:
        return 12.92 * c
    else:
        return (1.055 * np.power(c, 1.0 / 2.4) - 0.055)


def _to_linear(c):
    a = 0.055
    c2= np.array(c)
    ret = np.array(c) / 12.92
    ret[c2 > 0.04045] = np.power((c2[c2 > 0.04045] + a) / (1.0 + a), 2.4)
    return ret  

def _sc_to_linear(c):
    a = 0.055

    if c > 0.04045:
        return (np.power((c + a) / (1.0 + a), 2.4))
    else:
        return (c / 12.92)




def _sc_rgb_prepare(triple):
    ret = []
    for ch in triple:
        ch = round(ch, 3)

        if ch < -0.0001 or ch > 1.0001:
            raise Exception("Illegal RGB value %f" % ch)

        if ch < 0:
            ch = 0
        if ch > 1:
            ch = 1

        # Fix for Python 3 which by default rounds 4.5 down to 4.0
        # instead of Python 2 which is rounded to 5.0 which caused
        # a couple off by one errors in the tests. Tests now all pass
        # in Python 2 and Python 3
        ret.append(round(ch * 255 + 0.001))

    return ret


def hex_to_rgb(hex):
    if hex.startswith('#'):
        hex = hex[1:]
    r = int(hex[0:2], 16) / 255.0
    g = int(hex[2:4], 16) / 255.0
    b = int(hex[4:6], 16) / 255.0
    return [r, g, b]


def rgb_to_hex(triple):
    [r, g, b] = triple
    return '#%02x%02x%02x' % tuple(_sc_rgb_prepare([r, g, b]))

@triplescalar
def xyz_to_rgb(triple):  
    XYZ = np.asarray(triple)
    RGB=[np.sum(XYZ*np.array(mi).reshape((3,)+(XYZ.ndim-1)*(1,)),0) for mi in m]
    ret = _from_linear(RGB)
    return list(ret)
    
@triplescalar
def rgb_to_xyz(triple):
    RGB = _to_linear(triple)
    XYZ = [np.sum(RGB*np.array(mi).reshape((3,)+(RGB.ndim-1)*(1,)),0) for mi in m_inv]
    return list(XYZ)


def _sc_xyz_to_rgb(triple):  
    xyz = map(lambda row: _sc_dot_product(row, triple), m)
    return list(map(_sc_from_linear, xyz))


def _sc_rgb_to_xyz(triple):
    rgbl = list(map(_sc_to_linear, triple))
    return list(map(lambda row: _sc_dot_product(row, rgbl), m_inv))

@triplescalar
def xyz_to_luv(triple):
    X, Y, Z = triple

    mask1 = (X != 0.0) & (X != 0.0) & (Z != 0.0)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        varU = (4.0 * X) / (X + (15.0 * Y) + (3.0 * Z))
        varV = (9.0 * Y) / (X + (15.0 * Y) + (3.0 * Z))
        L = _f(Y)

    # Black will create a divide-by-zero error
    mask = (L!=0) & mask1
    L[~mask]=0.0
    U = 13.0 * L * (varU - refU)
    V = 13.0 * L * (varV - refV)
    U[~mask]=0.0
    V[~mask]=0.0
    return [L, U, V]


def _sc_xyz_to_luv(triple):
    X, Y, Z = triple

    if X == Y == Z == 0.0:
        return [0.0, 0.0, 0.0]

    varU = (4.0 * X) / (X + (15.0 * Y) + (3.0 * Z))
    varV = (9.0 * Y) / (X + (15.0 * Y) + (3.0 * Z))
    L = _sc_f(Y)

    # Black will create a divide-by-zero error
    if L == 0.0:
        return [0.0, 0.0, 0.0]

    U = 13.0 * L * (varU - refU)
    V = 13.0 * L * (varV - refV)

    return [L, U, V]

@triplescalar
def luv_to_xyz(triple):
    L, U, V = triple

    mask = (L == 0)
    # copy
    L2 = np.array(L).astype(np.float)
    L2[mask]=np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        varY = _f_inv(L2)
        varU = U / (13.0 * L2) + refU
        varV = V / (13.0 * L2) + refV
        Y = varY * refY
        X = 0.0 - (9.0 * Y * varU) / ((varU - 4.0) * varV - varU * varV)
        Z = (9.0 * Y - (15.0 * varV * Y) - (varV * X)) / (3.0 * varV)
    
    X[mask]=0.
    Y[mask]=0.
    Z[mask]=0.
    return [X, Y, Z]

def _sc_luv_to_xyz(triple):
    L, U, V = triple

    if L == 0:
        return [0.0, 0.0, 0.0]

    varY = _f_inv(L)
    varU = U / (13.0 * L) + refU
    varV = V / (13.0 * L) + refV
    Y = varY * refY
    X = 0.0 - (9.0 * Y * varU) / ((varU - 4.0) * varV - varU * varV)
    Z = (9.0 * Y - (15.0 * varV * Y) - (varV * X)) / (3.0 * varV)

    return [X, Y, Z]

@triplescalar
def luv_to_lch(triple):
    L, U, V = triple

    #C = (math.pow(math.pow(U, 2) + math.pow(V, 2), (1.0 / 2.0)))
    C=np.sqrt(U**2+V**2)
    #hrad = np.arctan2(V, U)
    H = np.degrees(np.arctan2(V, U)) % 360.

    return [L, C, H]

@triplescalar
def lch_to_luv(triple):
    L, C, H = triple

    Hrad = np.radians(H)
    U = np.cos(Hrad) * C
    V = np.sin(Hrad) * C

    return [L, U, V]
    

def _sc_luv_to_lch(triple):
    L, U, V = triple

    C = (math.pow(math.pow(U, 2) + math.pow(V, 2), (1.0 / 2.0)))
    hrad = (math.atan2(V, U))
    H = math.degrees(hrad)
    if H < 0.0:
        H = 360.0 + H

    return [L, C, H]


def _sc_lch_to_luv(triple):
    L, C, H = triple

    Hrad = math.radians(H)
    U = (math.cos(Hrad) * C)
    V = (math.sin(Hrad) * C)

    return [L, U, V]

@triplescalar
def husl_to_lch(triple):
    H, S, L = triple
    L=np.array(L)
    M1 =  (L > 99.9999999)
    M2 =  (L < 0.00000001)
    L[M1]=100.
    L[M2]=0.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        C = _max_chroma_for_LH(L, H)/ 100.0 * S
    
    C[M1 | M2]=0.

    return [L, C, H]

@triplescalar
def lch_to_husl(triple):
    L, C, H = triple
    L=np.array(L)
    M1 =  (L > 99.9999999)
    M2 =  (L < 0.00000001)
    L[M1]=100.
    L[M2]=0.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        S = C / _max_chroma_for_LH(L, H) * 100.0
    S[M1 | M2]=0.

    return [H, S, L]
    

def _sc_husl_to_lch(triple):
    H, S, L = triple

    if L > 99.9999999:
        return [100, 0.0, H]
    if L < 0.00000001:
        return [0.0, 0.0, H]

    mx = _sc_max_chroma_for_LH(L, H)
    C = mx / 100.0 * S

    return [L, C, H]


def _sc_lch_to_husl(triple):
    L, C, H = triple

    if L > 99.9999999:
        return [H, 0.0, 100.0]
    if L < 0.00000001:
        return [H, 0.0, 0.0]

    mx = _sc_max_chroma_for_LH(L, H)
    S = C / mx * 100.0

    return [H, S, L]

@triplescalar
def huslp_to_lch(triple):
    H, S, L = triple

    M1 =  (L > 99.9999999)
    M2 =  (L < 0.00000001)
    L[M1]=100.
    L[M2]=0.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        C = _max_safe_chroma_for_L(L) / 100.0 * S
    C[M1 | M2]=0.
    return [L, C, H]

@triplescalar
def lch_to_huslp(triple):
    L, C, H = triple

    M1 =  (L > 99.9999999)
    M2 =  (L < 0.00000001)
    L[M1]=100.
    L[M2]=0.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        S = C / _max_safe_chroma_for_L(L) * 100.0
    S[M1 | M2] = 0.
    
    return [H, S, L]
    

def _sc_huslp_to_lch(triple):
    H, S, L = triple

    if L > 99.9999999:
        return [100, 0.0, H]
    if L < 0.00000001:
        return [0.0, 0.0, H]

    mx = _sc_max_safe_chroma_for_L(L)
    C = mx / 100.0 * S

    return [L, C, H]


def _sc_lch_to_huslp(triple):
    L, C, H = triple

    if L > 99.9999999:
        return [H, 0.0, 100.0]
    if L < 0.00000001:
        return [H, 0.0, 0.0]

    mx = _sc_max_safe_chroma_for_L(L)
    S = C / mx * 100.0

    return [H, S, L]
