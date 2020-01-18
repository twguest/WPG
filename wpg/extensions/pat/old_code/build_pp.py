import shelve

d = shelve.open('shelve_data/asp_source_w_kb_optics_v4')  # load variables from shelf

pp = []
# NOTE: number is -1 to the image it creates ( KBM_S is element 11 but makes image 12
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# 0,1 White_B
pp.append([0, 0, 1.0, 0, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 2,3 HFM
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 0, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 4 DCM1
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 5,6 DCM2
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 0, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

# 7,8 Exit_S
#           01    2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 0, 0,          2,   .5,   2,   .5,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) # .5,2,.5,2
pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 9 Sec_Src_Ap
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 1, 0,          1,   20,   1,   20,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 10 Sec_SRC_W
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) #.5,3,.5,3
# 11 KBM_S
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 12 KBM_S_W
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 13,14 KBMV
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 0, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 15, 16 KBMH
#          0  1   2   1  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 0, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0,          .05,   20,   .05,  20,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# UpStreamW
# shift x,y
pp.append([0, 0, 1.0, 0, 0,          0.25,   4,   0.25,   4,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

# [ 0]: Auto-Resize (1) or not (0) Before propagation
# [ 1]: Auto-Resize (1) or not (0) After propagation
# [ 2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
# [ 3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
# [ 4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
# [ 5]: Horizontal Range modification factor at Resizing (1. means no modification)
# [ 6]: Horizontal Resolution modification factor at Resizing
# [ 7]: Vertical Range modification factor at Resizing
# [ 8]: Vertical Resolution modification factor at Resizing
# [ 9]: Type of wavefront Shift before Resizing (not yet implemented)
# [10]: New Horizontal wavefront Center position after Shift (not yet implemented)
# [11]: New Vertical wavefront Center position after Shift (not yet implemented)
# [12]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Horizontal Coordinate
# [13]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Vertical Coordinate
# [14]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Longitudinal Coordinate
# [15]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Horizontal Coordinate
# [16]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Vertical Coordinate

d['pp'] = pp

d.close()
