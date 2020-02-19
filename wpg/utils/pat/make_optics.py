#!/usr/bin/env python


# make_optics.py
# This script initializes optical elements to be simulated in the XFM beam line
# change parameters such as drift distances, apperture dimensions etc.


from PA_srw_wpg_utils import *

#
#
# def make_section1_optics():
#     el = []
#     rescale_param = []
#
#     ##### White Beam Slits: aperture 14.0m
#     el.append(srwlib.SRWLOptA("r", "a", 0.0006, 0.0006, 0.0, 0.0))
#     rescale_param.append([1, 1, 1, 1])
#
#     ##### Drift between White Beam Slits and HFM
#     el.append(srwlib.SRWLOptD(1.0))
#
#     rescale_param.append([1, 1, 1, 1])
#
#     ##### HFM: sphericalMirror 15.0m
#     el.append(my_hfm())
#     rescale_param.append([1, 1, 1, 1])
#
#     ##### Drift between HFM and DCM1
#     el.append(srwlib.SRWLOptD(1.0))
#     rescale_param.append([1, 1, 1, 1])
#
#     ##### DCM1: crystal 16.0m
#     el.append(my_dcm(1))
#     rescale_param.append([1, 1, 1, 1])
#
#     ##### Drift between DCM1 and DCM2
#     el.append(srwlib.SRWLOptD(0.1))
#     rescale_param.append([1, 1, 1, 1])
#
#     ##### DCM2: crystal 16.1m
#     el.append(my_dcm(2))
#     rescale_param.append([1, 1, 1, 1])
#
#     # Drift between DCM2 and Exit Slits
#     el.append(srwlib.SRWLOptD(2.0))
#     rescale_param.append([1, 1, 1, 1])
#
#     # Exit Slits: aperture 18.1m
#     el.append(srwlib.SRWLOptA("r", "a", 0.0006, 0.0006, 0.0, 0.0))  # .00025
#     rescale_param.append([1, 1, 1, 1])
#
#     # Drift between Exit Slits and Secondary Source Apperture
#     el.append(srwlib.SRWLOptD(13.0))
#     rescale_param.append([1, 1, 1, 1])
#
#     # Secondary Source Apperture: aperture 31.1m
#     el.append(srwlib.SRWLOptA("r", "a", 2e-05, 2e-05, 0.0, 0.0))
#     rescale_param.append([1, 1, 1, 1])
#
#     pps = pp_from_rescale_params(rescale_param)
#
#     return srwlib.SRWLOptC(el, pps)
#
#
# def make_section2_optics():
#     el = []
#     rescale_param = []
#
#     # Drift between Secondary Source Apperture and KBM_S1
#     el.append(srwlib.SRWLOptD(4.02))
#     rescale_param.append([1, 1, 1, 1])
#
#     # KBM_V: ellipsoidMirror 35.12m
#     el.append(my_kbmV(35.12, 0.45))
#     rescale_param.append([1, 1, 1, 1])
#
#     # Drift between KBM_V and KBM_H
#     el.append(srwlib.SRWLOptD(0.22))
#     rescale_param.append([1, 1, 1, 1])
#
#     # KBM_H: ellipsoidMirror 35.34m
#     el.append(my_kbmH(4.24, 0.23))
#     rescale_param.append([1, 1, 1, 1])
#
#     # Drift between KBM_H and KBM_S1
#     el.append(srwlib.SRWLOptD(0.18))
#     rescale_param.append([1, 1, 1, 1])
#
#     # KBM_S1: aperture 35.52m
#     el.append(srwlib.SRWLOptA("c", "a", 100e-06, 100e-06, 0.0, 0.0))
#     rescale_param.append([1, 1, 1, 1])
#
#     #  Drift from KBM_S1 to Focus sample at  35.57m
#     el.append(srwlib.SRWLOptD(0.05))
#     rescale_param.append([1, 1, 1, 1])
#
#     pps = pp_from_rescale_params(rescale_param)
#
#     return srwlib.SRWLOptC(el, pps)
#


def make_optics():
    # [ 5]: Horizontal Range modification factor at Resizing
    # [ 6]: Horizontal Resolution modification factor at Resizing
    # [ 7]: Vertical Range modification factor at Resizing
    # [ 8]: Vertical Resolution modification factor at Resizing
    el = []
    rescale_param = []

    # 1#### White Beam Slits: aperture 14.0m
    el.append(srwlib.SRWLOptA("r", "a", 1.1e-3, 0.2e-3, 0.0, 0.0))

    rescale_param.append([1.5, 4, 1.5, 4])

    # 2#### Drift between White Beam Slits and HFM
    el.append(srwlib.SRWLOptD(1.0))
    rescale_param.append([1, 1, 1, 1])

    # 3#### HFM: sphericalMirror 15.0m
    el.append(my_hfm())
    rescale_param.append([1, 1, 1, 1])

    # 4#### Drift between HFM and DCM1
    el.append(srwlib.SRWLOptD(1.0))
    rescale_param.append([1, 1, 1, 1])

    # 5#### DCM1: crystal 16.0m
    el.append(my_dcm(1))
    rescale_param.append([1, 1, 1, 1])

    # 6#### Drift between DCM1 and DCM2
    el.append(srwlib.SRWLOptD(0.1))
    rescale_param.append([1, 1, 1, 1])

    # 7#### DCM2: crystal 16.1m
    el.append(my_dcm(2))
    rescale_param.append([1, 1, 1, 1])

    # 8#### Drift between DCM2 and Exit Slits
    el.append(srwlib.SRWLOptD(2.0))
    rescale_param.append([1, 1, 1, 1])

    # 9#### Exit Slits: aperture 18.1m
    el.append(srwlib.SRWLOptA("r", "a", 600e-6, 600e-6, 0.0, 0.0))
    # el.append(srwlib.SRWLOptA("r", "a", 0.0006, 1, 0.0, 0.0))
    rescale_param.append([1, 1, 1, 1])  # rescale_param.append([.5, 2, .5, 2])

    # 10#### Drift between Exit Slits and Secondary Source Apperture
    el.append(srwlib.SRWLOptD(13.0))
    rescale_param.append(rescale_by_num(4))

    # 11#### Secondary Source Apperture: aperture 31.1m
    # el.append(srwlib.SRWLOptA("r", "a", 2.5e-6, 10e-6, 0.0, 0.0))
    el.append(srwlib.SRWLOptA("r", "a", 2.5e-6, 10e-6, 0.0, 0.0))
    rescale_param.append(rescale_by_num(0.01))

    # 12#### shift to center
    # el.append(srwlib.SRWLOptShift(_shift_x=-670e-9, _shift_y=-678.32e-9))
    el.append(srwlib.SRWLOptShift(_shift_x=0, _shift_y=0))
    rescale_param.append([3, 1.5, 5, 0.5])

    # 13, 14
    drift_vol_pp = [[10, 0.1, 6, 0.1], [3, 1, 1, 1]]  #
    drift_in_steps(4.02, drift_vol_pp, el, rescale_param)

    #
    # # 15#### KBM_VH Thin Lens
    # vkbfoc = 1. / (1. / 0.45 + 1. / 35.12)
    # hkbfoc = 1. / (1. / 0.45 + 1. / 4.24)
    # el.append(srwlib.SRWLOptL(_Fx=hkbfoc, _Fy=vkbfoc))
    # rescale_param.append([1, 1, 1, 1])

    # # 15#### KBM_V Thin Lens
    # vkbfoc = 1. / (1. / 0.45 + 1. / 35.12)
    # el.append(srwlib.SRWLOptL(_Fx =1e23 , _Fy=vkbfoc))
    # rescale_param.append([1,1, 0.5,2])
    #
    #
    # #
    # # 16#### Drift between KBM_V and KBM_H
    # el.append(srwlib.SRWLOptD(0.22))
    # rescale_param.append([1,1,1,1])
    #
    # #
    # #
    # #
    # # 17#### KBM_H Thin Lens
    # hkbfoc = 1. / (1. / 0.23 + 1. / 4.24)
    # el.append(srwlib.SRWLOptL(_Fx=hkbfoc, _Fy=1e23))
    # #rescale_param.append([1,2.5,1,0.75])
    # rescale_param.append([2, 0.125, 1.5, 1.5])

    # 15#### KBM_V: ellipsoidMirror 35.12m
    el.append(my_kbmV(35.12, 0.45))
    # el.append(my_kbmV(35.12, 0.46))
    # rescale_param.append([0.145, 1/0.16, 0.2, 1/0.2])
    rescale_param.append([0.16, 1 / 0.16, 0.2, 1 / 0.2])
    # rescale_param.append(rescale_by_num(0.145))

    # 16#### Drift between KBM_V and KBM_H
    el.append(srwlib.SRWLOptD(0.22))
    rescale_param.append(rescale_by_num(1.25))

    # 17#### KBM_H: ellipsoidMirror 35.34m
    # el.append(my_kbmH(4.24, 0.233))
    el.append(my_kbmH(4.24, 0.233))
    rescale_param.append(rescale_by_num(0.4))

    # 18-25 onwards ###### drift to sample
    #    drift_vol_pp = [[1, 1, 1,1],
    #                    [1, 0.5, 1, 1],
    #                    [1, 0.5, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    rescale_by_num(1),
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    rescale_by_num(1),
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    [1, 1, 1, 1],
    #                    ]
    #    drift_vol_pp = [[1, 1, 1,1]]
    #    drift_in_steps(0.27, drift_vol_pp, el, rescale_param)

    el.append(srwlib.SRWLOptD(0.27))
    """rescale_param.append(Use_PP( 
                 semi_analytical_treatment=1,
                 zoom_h=0.5, zoom_v=0.5,
                 sampling_h=None, sampling_v=2))
    """
    rescale_param.append([1, 0.5, 0.5, 1, 2])

    #  rescale_param.append([1,0.5,1,2,2])

    ###0.27 for the focus in the thesis.
    # drift_in_steps(0.25, drift_vol_pp, el, rescale_param)

    # #
    # drift_vol_pp = [[1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1]
    #                 ]
    # #
    # drift_in_steps(0.5, drift_vol_pp, el, rescale_param)

    # drift_vol_pp = [[1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 [1, 1, 1, 1],
    #                 ]
    # drift_in_steps(0.44, drift_vol_pp, el, rescale_param)

    # # 17#### KBM_S1: aperture 35.52m
    # el.append(srwlib.SRWLOptA("c", "a", 100e-06, 100e-06, 0.0, 0.0))
    # rescale_param.append([0.25, 4, 0.25, 4])
    #
    # # 16#### Drift between KBM_H and KBM_S1
    # el.append(srwlib.SRWLOptD(0.18))
    # rescale_param.append([1, 1, 1, 1])
    #
    #
    # # 18####  Drift from KBM_S1 to Focus sample at  35.57m
    # el.append(srwlib.SRWLOptD(0.05))
    # rescale_param.append([1, 1, 1, 1])

    pps = pp_from_rescale_params(rescale_param)

    # TRY OPENING THE WBS to get full vertical beam!!
    return srwlib.SRWLOptC(el, pps)


def rescale_by_num(num):
    ls = [num, 1 / float(num), num, 1 / float(num)]
    return ls


def pp_from_rescale_params(rescale_params):

    a = [0, 0, 1, 0, 0]
    c = [0, 0, 0, 0, 0, 0, 0, 0]
    d = [0, 0, 1, 0]

    pps = []
    for rescale_param in rescale_params:
        if len(rescale_param) == 4:
            pp = PA_concat_list(a, rescale_param, c)
        elif len(rescale_param) == 5:
            pp = PA_concat_list(d, rescale_param, c)
        pps.append(pp)

    return pps


def drift_in_steps(d, rescale_params, el, pp):
    # 12#### Drift between Secondary Source Apperture and KBM_V
    x = d / len(rescale_params)

    for i, rescale_param in enumerate(rescale_params):
        el.append(srwlib.SRWLOptD(x))
        pp.append(rescale_param)


def my_hfm():
    HFM = srwlib.SRWLOptMirSph(
        _r=8871.45,
        _size_tang=0.95,
        _size_sag=0.005,
        _nvx=0.999902524484,
        _nvy=0.0,
        _nvz=-0.0139621463437,
        _tvx=0.0139621463437,
        _tvy=0.0,
        _x=0.0,
        _y=0.0,
    )
    return HFM


def my_dcm(orient):
    # opCr = srwlib.SRWLOptCryst(_d_sp=3.13557135638, _psi0r=-1.43534717388e-05, _psi0i=3.16676479565e-07,
    #                            _psi_hr=-7.5926201741e-06, _psi_hi=2.21095173107e-07, _psi_hbr=-7.5926201741e-06,
    #                            _psi_hbi=2.21095173107e-07, _tc=0.01, _ang_as=0.0)

    opCr = srwlib.SRWLOptCryst(
        _d_sp=3.13557135638,
        _psi0r=-1.80680829752e-05,
        _psi0i=4.94786384696e-07,
        _psi_hr=-9.56518891953e-06,
        _psi_hi=3.45446815392e-07,
        _psi_hbr=-9.56518891953e-06,
        _psi_hbi=3.45446815392e-07,
        _tc=0.01,
        _ang_as=0.0,
    )

    if orient == 1:
        # opCr.set_orient(-0.970946542231, 6.59748146207e-09, -0.239296494187, -0.239296494187, 1.62599496025e-09)
        opCr.set_orient(
            -0.258293438037,
            0.928106834743,
            -0.268145861745,
            -0.0718931646453,
            0.258328426713,
        )
    else:
        # opCr.set_orient(-3.48549712965e-09, -0.970946542231, -0.239296494187, -8.59024886897e-10, -0.239296494187)
        opCr.set_orient(
            -0.258293438037,
            0.928106834743,
            -0.268145861745,
            -0.0718931646453,
            0.258328426713,
        )

    return opCr


def my_kbmV(source_focus, image_focus):
    # mirror = srwlib.SRWLOptMirEl(_p=source_focus, _q=image_focus, _ang_graz=0.003, _size_tang=0.22, _size_sag=0.22,
    #                              _nvx=0.999995500003, _nvy=0.0, _nvz=-0.0029999955, _tvx=0.0029999955,
    #                              _tvy=0.0, _x=0.0, _y=0.0)

    mirror = srwlib.SRWLOptMirEl(
        _p=source_focus,
        _q=image_focus,
        _ang_graz=0.003,
        _size_tang=0.22,
        _size_sag=1e6,
        _nvx=0.999995500003,
        _nvy=0.0,
        _nvz=-0.0029999955,
        _tvx=0.0029999955,
        _tvy=0.0,
        _x=0.0,
        _y=0.0,
    )

    return mirror


def my_kbmH(source_focus, image_focus):
    # mirror = srwlib.SRWLOptMirEl(_p=source_focus, _q=image_focus, _ang_graz=0.003, _size_tang=0.22, _size_sag=0.22,
    #                              _nvx=0.0, _nvy=0.999995500003, _nvz=-0.0029999955, _tvx=0.0,
    #                              _tvy=0.0029999955, _x=0.0, _y=0.0)

    mirror = srwlib.SRWLOptMirEl(
        _p=source_focus,
        _q=image_focus,
        _ang_graz=0.003,
        _size_tang=0.22,
        _size_sag=1e6,
        _nvx=0.0,
        _nvy=0.999995500003,
        _nvz=-0.0029999955,
        _tvx=0.0,
        _tvy=0.0029999955,
        _x=0.0,
        _y=0.0,
    )
    return mirror


def main():
    op = make_optics()

    d = PA_open_shelf("optics")

    d["op"] = op  # save the optics

    d.close()  # close the shelf


if __name__ == "__main__":
    main()
