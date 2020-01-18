# build_pp.py
# Build the propagation parameters for each optical elements
# define wavefeild scaling and resolution after each element

import shelve
linux = True

if linux:
    d = shelve.open('shelve_data/asp_source_w_kb_optics_v4')  # load variables from shelf
else:
    d = shelve.open('shelve_data\\asp_source_w_kb_optics_v4')  # load variables from shelf
pp = []

#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# 1,2 White_B:
pp.append([0, 0, 1.0, 0, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 3,4 HFM
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 0, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 5 DCM1
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 1, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 6,7 DCM2
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 0, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

# 8,9 Exit_S
#           01    2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
#pp.append([0, 0, 1.0, 0, 0, .5, 2, .5, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # .5,2,.5,2
pp.append([0, 0, 1.0, 1, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 10 Sec_Src_Ap
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
#pp.append([0, 0, 1.0, 1, 0, 0.2, 5, 0.2, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0, .25, 5, .25, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 11 Sec_SRC_W
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 1, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 12 KBM_S
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 1, 0, 0.25, 4, 0.25, 4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 13 KBM_S_W
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 1, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 14,15 KBMV
#          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 0, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
pp.append([0, 0, 1.0, 1, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# 16, 17 KBMH
#          0  1   2   1  4           5    6    7    8          9    10   11   12   13  14    15    16
pp.append([0, 0, 1.0, 0, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#pp.append([0, 0, 1.0, 1, 0,.1, 10, .1, 10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] )

pp.append([0, 0, 1.0, 1, 0,.25, 4, .25, 4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# UpStreamW

# #          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 0, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# # 2,3 HFM
# #          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 0, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# # 4 DCM1
# #          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# # 5,6 DCM2
# #          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 0, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#
# # 7,8 Exit_S
# #           01    2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 0, 0,          2,   .5,   2,   .5,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) # .5,2,.5,2
# pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# # 9 Sec_Src_Ap
# #          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 1, 0,          1,   2,   1,   2,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# # 10 Sec_SRC_W
# #          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) #.5,3,.5,3
# # 11 KBM_S
# #          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# # 12 KBM_S_W
# #          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# # 13,14 KBMV
# #          0  1   2   3  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 0, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# pp.append([0, 0, 1.0, 1, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# # 15, 16 KBMH
# #          0  1   2   1  4           5    6    7    8          9    10   11   12   13  14    15    16
# pp.append([0, 0, 1.0, 0, 0,          1,   1,   1,   1,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# pp.append([0, 0, 1.0, 1, 0,          .5,   2,   .5,  2,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# # UpStreamW


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
