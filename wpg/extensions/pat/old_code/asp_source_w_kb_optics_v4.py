#!/usr/bin/env python


# aps_source_..v4.py
# This script initializes optical elements to be simulated in the XFM beam line
# change parameters such as drift distances, apperture dimensions etc.

# Imports (Some might not be used)
import os, sys
import shelve

sys.path.append(r'WPG')
sys.path.append(r'/opt/wpg/wpg')

linux = True


from wpg import srwl_bl, \
    srwlib, \
    srwl_uti_smp, \
    uti_io, \
    srwl_uti_src, \
    wavefront, \
    optical_elements, \
    beamline, \
    generators, \
    wpg_uti_wf, \
    srwlpy


def set_optics(v=None):
    el = []
    # 1,2 White_B: aperture 14.0m
    el.append(srwlib.SRWLOptA("r", "a", 0.0006, 0.0006, 0.0, 0.0))  # .0004
    el.append(srwlib.SRWLOptD(1.0))

    # 3,4 HFM: sphericalMirror 15.0m
    el.append(srwlib.SRWLOptMirSph(_r=8871.45, _size_tang=0.95, _size_sag=0.005, _nvx=0.999902524484, _nvy=0.0,
                                   _nvz=-0.0139621463437, _tvx=0.0139621463437, _tvy=0.0, _x=0.0, _y=0.0))

    el.append(srwlib.SRWLOptD(1.0))
    # 5 DCM1: crystal 16.0m
    opCr = srwlib.SRWLOptCryst(_d_sp=3.13557135638, _psi0r=-1.43534717388e-05, _psi0i=3.16676479565e-07,
                               _psi_hr=-7.5926201741e-06, _psi_hi=2.21095173107e-07, _psi_hbr=-7.5926201741e-06,
                               _psi_hbi=2.21095173107e-07, _tc=0.01, _ang_as=0.0)
    # Set crystal orientation:
    opCr.set_orient(-0.970946542231, 6.59748146207e-09, -0.239296494187, -0.239296494187, 1.62599496025e-09)
    el.append(opCr)

    # 6,7 DCM2: crystal 16.0m
    opCr = srwlib.SRWLOptCryst(_d_sp=3.13557135638, _psi0r=-1.43534717388e-05, _psi0i=3.16676479565e-07,
                               _psi_hr=-7.5926201741e-06, _psi_hi=2.21095173107e-07, _psi_hbr=-7.5926201741e-06,
                               _psi_hbi=2.21095173107e-07, _tc=0.01, _ang_as=0.0)
    #  Set crystal orientation:
    opCr.set_orient(-3.48549712965e-09, -0.970946542231, -0.239296494187, -8.59024886897e-10, -0.239296494187)
    el.append(opCr)

    el.append(srwlib.SRWLOptD(2.0))
    # 8,9  Exit_S: aperture 18.0m
    el.append(srwlib.SRWLOptA("r", "a", 0.0006, 0.0006, 0.0, 0.0))  # .00025
    el.append(srwlib.SRWLOptD(12.0))
    # 10  Sec_Src_Ap: aperture 30.0m
    el.append(srwlib.SRWLOptA("r", "a", 2e-05, 2e-05, 0.0, 0.0))
    # 11 Sec_SRC_W: watch 30.0m

    el.append(srwlib.SRWLOptD(4.99)) #4.9
    # 12 KBM_S: aperture 34.9m
    el.append(srwlib.SRWLOptA("r", "a", 2e-05, 2e-05, 0.0, 0.0))
    #el.append(srwlib.SRWLOptA("r", "a", 1000e-05, 1000e-05, 0.0, 0.0))
    # 13 KBM_S_W: watch 34.9m

    el.append(srwlib.SRWLOptD(0.01)) #0.1
    # 14, 15 KBMV: ellipsoidMirror 35.0m
    el.append(
        srwlib.SRWLOptMirEl(_p=35.0, _q=0.47, _ang_graz=0.003, _size_tang=0.22, _size_sag=0.22, _nvx=0.999995500003,
                            _nvy=0.0, _nvz=-0.0029999955, _tvx=0.0029999955, _tvy=0.0, _x=0.0, _y=0.0))

    el.append(srwlib.SRWLOptD(0.22))
    # 16,17 KBMH: ellipsoidMirror 35.22m
    el.append(srwlib.SRWLOptMirEl(_p=5.22, _q=0.25, _ang_graz=0.003, _size_tang=0.22, _size_sag=0.22, _nvx=0.0,
                                  _nvy=0.999995500003, _nvz=-0.0029999955, _tvx=0.0, _tvy=0.0029999955, _x=0.0, _y=0.0))

    el.append(srwlib.SRWLOptD(0.25))

    return srwlib.SRWLOptC(el)


varParam = srwl_bl.srwl_uti_ext_options([
    ['name', 's', 'asp_source_w_kb_optics_v4', 'simulation name'],

    # ---Data Folder
    ['fdir', 's', '', 'folder (directory) name for reading-in input and saving output data files'],

    # ---Electron Beam
    ['ebm_nm', 's', '', 'standard electron beam name'],
    ['ebm_nms', 's', '', 'standard electron beam name suffix: e.g. can be Day1, Final'],
    ['ebm_i', 'f', 0.2, 'electron beam current [A]'],
    ['ebm_e', 'f', 3.0, 'electron beam avarage energy [GeV]'],
    ['ebm_de', 'f', 0.0, 'electron beam average energy deviation [GeV]'],
    ['ebm_x', 'f', 0.0, 'electron beam initial average horizontal position [m]'],
    ['ebm_y', 'f', 0.0, 'electron beam initial average vertical position [m]'],
    ['ebm_xp', 'f', 0.0, 'electron beam initial average horizontal angle [rad]'],
    ['ebm_yp', 'f', 0.0, 'electron beam initial average vertical angle [rad]'],
    ['ebm_z', 'f', 0., 'electron beam initial average longitudinal position [m]'],
    ['ebm_dr', 'f', -1.034, 'electron beam longitudinal drift [m] to be performed before a required calculation'],
    ['ebm_ens', 'f', 0.001028, 'electron beam relative energy spread'],
    ['ebm_emx', 'f', 1.0198e-08, 'electron beam horizontal emittance [m]'],
    ['ebm_emy', 'f', 1.0198e-10, 'electron beam vertical emittance [m]'],
    # Definition of the beam through Twiss:
    ['ebm_betax', 'f', 2.452, 'horizontal beta-function [m]'],
    ['ebm_betay', 'f', 2.452, 'vertical beta-function [m]'],
    ['ebm_alphax', 'f', 0.0, 'horizontal alpha-function [rad]'],
    ['ebm_alphay', 'f', 0.0, 'vertical alpha-function [rad]'],
    ['ebm_etax', 'f', 0.1, 'horizontal dispersion function [m]'],
    ['ebm_etay', 'f', 0.0, 'vertical dispersion function [m]'],
    ['ebm_etaxp', 'f', 0.0, 'horizontal dispersion function derivative [rad]'],
    ['ebm_etayp', 'f', 0.0, 'vertical dispersion function derivative [rad]'],
    # Definition of the beam through Moments:
    ['ebm_sigx', 'f', 0.000188608949947, 'horizontal RMS size of electron beam [m]'],
    ['ebm_sigy', 'f', 1.58131261931e-05, 'vertical RMS size of electron beam [m]'],
    ['ebm_sigxp', 'f', 6.44907267257e-05, 'horizontal RMS angular divergence of electron beam [rad]'],
    ['ebm_sigyp', 'f', 6.44907267257e-06, 'vertical RMS angular divergence of electron beam [rad]'],
    ['ebm_mxxp', 'f', 0.0, 'horizontal position-angle mixed 2nd order moment of electron beam [m]'],
    ['ebm_myyp', 'f', 0.0, 'vertical position-angle mixed 2nd order moment of electron beam [m]'],

    # ---Undulator
    ['und_bx', 'f', 0.0, 'undulator horizontal peak magnetic field [T]'],
    ['und_by', 'f', 0.8, 'undulator vertical peak magnetic field [T]'],
    ['und_phx', 'f', 0.0, 'initial phase of the horizontal magnetic field [rad]'],
    ['und_phy', 'f', 0.0, 'initial phase of the vertical magnetic field [rad]'],
    ['und_b2e', '', '',
     'estimate undulator fundamental photon energy (in [eV]) for the amplitude of sinusoidal magnetic field defined by und_b or und_bx, und_by',
     'store_true'],
    ['und_e2b', '', '', 'estimate undulator field amplitude (in [T]) for the photon energy defined by w_e',
     'store_true'],
    ['und_per', 'f', 0.022, 'undulator period [m]'],
    ['und_len', 'f', 1.98, 'undulator length [m]'],
    ['und_zc', 'f', 0.0, 'undulator center longitudinal position [m]'],
    ['und_sx', 'i', 1, 'undulator horizontal magnetic field symmetry vs longitudinal position'],
    ['und_sy', 'i', -1, 'undulator vertical magnetic field symmetry vs longitudinal position'],
    ['und_g', 'f', 6.72, 'undulator gap [mm] (assumes availability of magnetic measurement or simulation data)'],
    ['und_ph', 'f', 0.0, 'shift of magnet arrays [mm] for which the field should be set up'],
    ['und_mdir', 's', '', 'name of magnetic measurements sub-folder'],
    ['und_mfs', 's', '', 'name of magnetic measurements for different gaps summary file'],

    # ---Calculation Types
    # Electron Trajectory
    ['tr', '', '', 'calculate electron trajectory', 'store_true'],
    ['tr_cti', 'f', 0.0, 'initial time moment (c*t) for electron trajectory calculation [m]'],
    ['tr_ctf', 'f', 0.0, 'final time moment (c*t) for electron trajectory calculation [m]'],
    ['tr_np', 'f', 10000, 'number of points for trajectory calculation'],
    ['tr_mag', 'i', 1, 'magnetic field to be used for trajectory calculation: 1- approximate, 2- accurate'],
    ['tr_fn', 's', 'res_trj.dat', 'file name for saving calculated trajectory data'],
    ['tr_pl', 's', '',
     'plot the resulting trajectiry in graph(s): ""- dont plot, otherwise the string should list the trajectory components to plot'],

    # Single-Electron Spectrum vs Photon Energy
    ['ss', '', '', 'calculate single-e spectrum vs photon energy', 'store_true'],
    ['ss_ei', 'f', 100.0, 'initial photon energy [eV] for single-e spectrum vs photon energy calculation'],
    ['ss_ef', 'f', 20000.0, 'final photon energy [eV] for single-e spectrum vs photon energy calculation'],
    ['ss_ne', 'i', 10000, 'number of points vs photon energy for single-e spectrum vs photon energy calculation'],
    ['ss_x', 'f', 0.0, 'horizontal position [m] for single-e spectrum vs photon energy calculation'],
    ['ss_y', 'f', 0.0, 'vertical position [m] for single-e spectrum vs photon energy calculation'],
    ['ss_meth', 'i', 1,
     'method to use for single-e spectrum vs photon energy calculation: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['ss_prec', 'f', 0.01,
     'relative precision for single-e spectrum vs photon energy calculation (nominal value is 0.01)'],
    ['ss_pol', 'i', 6,
     'polarization component to extract after spectrum vs photon energy calculation: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['ss_mag', 'i', 1,
     'magnetic field to be used for single-e spectrum vs photon energy calculation: 1- approximate, 2- accurate'],
    ['ss_ft', 's', 'f', 'presentation/domain: "f"- frequency (photon energy), "t"- time'],
    ['ss_u', 'i', 1,
     'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['ss_fn', 's', 'res_spec_se.dat', 'file name for saving calculated single-e spectrum vs photon energy'],
    ['ss_pl', 's', '',
     'plot the resulting single-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],

    # Multi-Electron Spectrum vs Photon Energy (taking into account e-beam emittance, energy spread and collection aperture size)
    ['sm', '', '', 'calculate multi-e spectrum vs photon energy', 'store_true'],
    ['sm_ei', 'f', 100.0, 'initial photon energy [eV] for multi-e spectrum vs photon energy calculation'],
    ['sm_ef', 'f', 20000.0, 'final photon energy [eV] for multi-e spectrum vs photon energy calculation'],
    ['sm_ne', 'i', 10000, 'number of points vs photon energy for multi-e spectrum vs photon energy calculation'],
    ['sm_x', 'f', 0.0, 'horizontal center position [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_rx', 'f', 0.001,
     'range of horizontal position / horizontal aperture size [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_nx', 'i', 1, 'number of points vs horizontal position for multi-e spectrum vs photon energy calculation'],
    ['sm_y', 'f', 0.0, 'vertical center position [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_ry', 'f', 0.001,
     'range of vertical position / vertical aperture size [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_ny', 'i', 1, 'number of points vs vertical position for multi-e spectrum vs photon energy calculation'],
    ['sm_mag', 'i', 1,
     'magnetic field to be used for calculation of multi-e spectrum spectrum or intensity distribution: 1- approximate, 2- accurate'],
    ['sm_hi', 'i', 1,
     'initial UR spectral harmonic to be taken into account for multi-e spectrum vs photon energy calculation'],
    ['sm_hf', 'i', 15,
     'final UR spectral harmonic to be taken into account for multi-e spectrum vs photon energy calculation'],
    ['sm_prl', 'f', 1.0,
     'longitudinal integration precision parameter for multi-e spectrum vs photon energy calculation'],
    ['sm_pra', 'f', 1.0, 'azimuthal integration precision parameter for multi-e spectrum vs photon energy calculation'],
    ['sm_meth', 'i', -1,
     'method to use for spectrum vs photon energy calculation in case of arbitrary input magnetic field: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler", -1- dont use this accurate integration method (rather use approximate if possible)'],
    ['sm_prec', 'f', 0.01,
     'relative precision for spectrum vs photon energy calculation in case of arbitrary input magnetic field (nominal value is 0.01)'],
    ['sm_nm', 'i', 1,
     'number of macro-electrons for calculation of spectrum in case of arbitrary input magnetic field'],
    ['sm_na', 'i', 5,
     'number of macro-electrons to average on each node at parallel (MPI-based) calculation of spectrum in case of arbitrary input magnetic field'],
    ['sm_ns', 'i', 5,
     'saving periodicity (in terms of macro-electrons) for intermediate intensity at calculation of multi-electron spectrum in case of arbitrary input magnetic field'],
    ['sm_type', 'i', 1, 'calculate flux (=1) or flux per unit surface (=2)'],
    ['sm_pol', 'i', 6,
     'polarization component to extract after calculation of multi-e flux or intensity: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['sm_rm', 'i', 1,
     'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
    ['sm_fn', 's', 'res_spec_me.dat', 'file name for saving calculated milti-e spectrum vs photon energy'],
    ['sm_pl', 's', '',
     'plot the resulting spectrum-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],
    # to add options for the multi-e calculation from "accurate" magnetic field

    # Power Density Distribution vs horizontal and vertical position
    ['pw', '', '', 'calculate SR power density distribution', 'store_true'],
    ['pw_x', 'f', 0.0,
     'central horizontal position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_rx', 'f', 0.015,
     'range of horizontal position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_nx', 'i', 100, 'number of points vs horizontal position for calculation of power density distribution'],
    ['pw_y', 'f', 0.0,
     'central vertical position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_ry', 'f', 0.015,
     'range of vertical position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_ny', 'i', 100, 'number of points vs vertical position for calculation of power density distribution'],
    ['pw_pr', 'f', 1.0, 'precision factor for calculation of power density distribution'],
    ['pw_meth', 'i', 1, 'power density computation method (1- "near field", 2- "far field")'],
    ['pw_zst', 'f', 0.,
     'initial longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
    ['pw_zfi', 'f', 0.,
     'final longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
    ['pw_mag', 'i', 1, 'magnetic field to be used for power density calculation: 1- approximate, 2- accurate'],
    ['pw_fn', 's', 'res_pow.dat', 'file name for saving calculated power density distribution'],
    ['pw_pl', 's', '',
     'plot the resulting power density distribution in a graph: ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    # Single-Electron Intensity distribution vs horizontal and vertical position
    ['si', '', '',
     'calculate single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position',
     'store_true'],
    # Single-Electron Wavefront Propagation
    ['ws', '', '', 'calculate single-electron (/ fully coherent) wavefront propagation', 'store_true'],
    # Multi-Electron (partially-coherent) Wavefront Propagation
    ['wm', '', '', 'calculate multi-electron (/ partially coherent) wavefront propagation', 'store_true'],

    ['w_e', 'f', 8263.0,
     'photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ef', 'f', -1.0,
     'final photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ne', 'i', 1, 'number of points vs photon energy for calculation of intensity distribution'],
    ['w_x', 'f', 0.0, 'central horizontal position [m] for calculation of intensity distribution'],
    ['w_rx', 'f', 0.0025, 'range of horizontal position [m] for calculation of intensity distribution'],
    ['w_nx', 'i', 100, 'number of points vs horizontal position for calculation of intensity distribution'],
    ['w_y', 'f', 0.0,
     'central vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ry', 'f', 0.0025,
     'range of vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ny', 'i', 100, 'number of points vs vertical position for calculation of intensity distribution'],
    ['w_smpf', 'f', 0.1,
     'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_meth', 'i', 1, 'method to use for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_prec', 'f', 0.01,
     'relative precision for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_u', 'i', 1,
     'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['si_pol', 'i', 6,
     'polarization component to extract after calculation of intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['si_type', 'i', 0,
     'type of a characteristic to be extracted after calculation of intensity distribution: 0- Single-Electron Intensity, 1- Multi-Electron Intensity, 2- Single-Electron Flux, 3- Multi-Electron Flux, 4- Single-Electron Radiation Phase, 5- Re(E): Real part of Single-Electron Electric Field, 6- Im(E): Imaginary part of Single-Electron Electric Field, 7- Single-Electron Intensity, integrated over Time or Photon Energy'],
    ['w_mag', 'i', 1,
     'magnetic field to be used for calculation of intensity distribution vs horizontal and vertical position: 1- approximate, 2- accurate'],

    ['si_fn', 's', 'res_int_se.dat',
     'file name for saving calculated single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position'],
    ['si_pl', 's', '',
     'plot the input intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],
    ['ws_fni', 's', 'res_int_pr_se.dat',
     'file name for saving propagated single-e intensity distribution vs horizontal and vertical position'],
    ['ws_pl', 's', '',
     'plot the resulting intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    ['wm_nm', 'i', 100000,
     'number of macro-electrons (coherent wavefronts) for calculation of multi-electron wavefront propagation'],
    ['wm_na', 'i', 5,
     'number of macro-electrons (coherent wavefronts) to average on each node for parallel (MPI-based) calculation of multi-electron wavefront propagation'],
    ['wm_ns', 'i', 5,
     'saving periodicity (in terms of macro-electrons / coherent wavefronts) for intermediate intensity at multi-electron wavefront propagation calculation'],
    ['wm_ch', 'i', 0,
     'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y'],
    ['wm_ap', 'i', 0,
     'switch specifying representation of the resulting Stokes parameters: coordinate (0) or angular (1)'],
    ['wm_x0', 'f', 0, 'horizontal center position for mutual intensity cut calculation'],
    ['wm_y0', 'f', 0, 'vertical center position for mutual intensity cut calculation'],
    ['wm_ei', 'i', 0,
     'integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from w_e, w_ef'],
    ['wm_rm', 'i', 1,
     'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
    ['wm_fni', 's', 'res_int_pr_me.dat',
     'file name for saving propagated multi-e intensity distribution vs horizontal and vertical position'],

    # to add options
    ['op_r', 'f', 14.0, 'longitudinal position of the first optical element [m]'],

    # Former appParam:
    ['source_type', 's', 'u',
     'source type, (u) idealized undulator, (t), tabulated undulator, (m) multipole, (g) gaussian beam'],
])


def main():
    v = srwl_bl.srwl_uti_parse_options(varParam,
                                       use_sys_argv=True)  # Make a namespace (similar to dictionary) for the parameters

    source_type, mag = srwl_bl.setup_source(v)  # Note sure if this is actually required?

    op = set_optics(v)  # Make the container or optics

    return v, source_type, mag, op


if __name__ == '__main__':
    v, source_type, mag, op = main()  # Get the things made in main()

    if linux:
        d = shelve.open('shelve_data/asp_source_w_kb_optics_v4')  # Open the shelf doc for this instance
    else:
        d = shelve.open('shelve_data\\asp_source_w_kb_optics_v4')  # Open the shelf doc for this instance

    d['v'] = v  # save the parameters

    d['op'] = op  # save the optics

    d.close()  # close the shelf
