# Imports
import os, sys
sys.path.append('/opt/wpg/wpg')
import matplotlib.pyplot as plt
import numpy as np
import time as t
import shelve
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
import array as ar



if __name__ == '__main__':

    pixels = [100, 250, 512, 720, 1024, 2048]  # Define the pixel dimensions of the initial wave field

    wfs = []  # initialize empty list of srw wavefront objects
    d = shelve.open('shelve_data/asp_source_w_kb_optics_v4')  # Open the shelf for an instance
    v = d['v']  # download the parameters and optics from the shelf
    
    

    elecBeam = srwlib.SRWLPartBeam()  # initialize a particle beam

    elecBeam.from_Twiss(_Iavg=v.ebm_i,  # twiss parameters are saved in the v object namespace
                        _e=v.ebm_e,
                        _sig_e=v.ebm_de,
                        _emit_x=v.ebm_emx,
                        _beta_x=v.ebm_betax,
                        _alpha_x=v.ebm_alphax,
                        _eta_x=v.ebm_etax,
                        _eta_x_pr=v.ebm_etaxp,
                        _emit_y=v.ebm_emy,
                        _beta_y=v.ebm_betay,
                        _alpha_y=v.ebm_alphay,
                        _eta_y=v.ebm_etay,
                        _eta_y_pr=v.ebm_etayp)

    # Some more undulator definitions
    numPer = v.und_len / v.und_per  # Number of ID Periods (without counting for terminations
    undPer = v.und_per  # Period Length [m]
    Bx = v.und_bx  # Peak Horizontal field [T]
    By = v.und_by  # Peak Vertical field [T]
    phBx = v.und_phx  # Initial Phase of the Horizontal field component
    phBy = v.und_phy  # Initial Phase of the Vertical field component
    sBx = v.und_sx  # Symmetry of the Horizontal field component vs Longitudinal position
    sBy = v.und_sy  # Symmetry of the Vertical field component vs Longitudinal position
    xcID = 0  # Transverse Coordinates of Undulator Center [m]
    ycID = 0

    # Z center of the insertion device
    zcID = -numPer * undPer  # Longitudinal Coordinate of Undulator Center wit hrespect to Straight Section Center [m]

    # Unclear how this works, should be - and on the order of the total size of the undulator ~2m
    # zcID = -2.012  #Longitudinal Coordinate of Undulator Center wit hrespect to Straight Section Center [m]
    # zcID = -numPer*undPer/2  #Longitudinal Coordinate of Undulator Center wit hrespect to Straight Section Center [m]

    # Make the undulator
    und = srwlib.SRWLMagFldU(
        [srwlib.SRWLMagFldH(1, 'v', By, phBy, sBy, 1), srwlib.SRWLMagFldH(1, 'h', Bx, phBx, sBx, 1)], undPer,
        numPer)  # Ellipsoidal Undulator
    magFldCnt = srwlib.SRWLMagFldC([und], ar.array('d', [xcID]), ar.array('d', [ycID]),
                                   ar.array('d', [zcID]))  # Container of all Field Elements

    # ***********Precision Parameters for SR calculation
    # Change here for calculation changes.
    meth = 1  # SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
    relPrec = 1  # relative precision
    zStartInteg = 0  # longitudinal position to start integration (effective if < zEndInteg)
    zEndInteg = 0  # longitudinal position to finish integration (effective if > zStartInteg)
    npTraj = 10000  # v.tr_np #Number of points for trajectory calculation
    useTermin = 1  # 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
    sampFactNxNyForProp = 0  # sampling factor for adjusting nx, ny (effective if > 0)
    arPrecPar = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, sampFactNxNyForProp]
    for i in pixels:  # create a wavefront for every option in pixels.

        time_start = t.time()  # Start a timer



        wf1 = srwlib.SRWLWfr()  # intialize wave feild object

        # ***********Initial Wavefront data placeholder
        wf1.allocate(1, i, i)  # Numbers of points vs Photon Energy, Horizontal and Vertical Positions
        # wf1.allocate(1, 720, 720)  # Numbers of points vs Photon Energy, Horizontal and Vertical Positions
        # wf1.mesh.zStart = 1
        wf1.mesh.zStart = 14  # Initial drift to first element
        # Longitudinal Position [m] from Center of Straight Section at which SR has to be calculated
        wf1.mesh.eStart = v.w_e  # Initial Photon Energy [eV]
        wf1.mesh.eFin = v.w_e  # Final Photon Energy [eV]
        wf1.mesh.xStart = -1e-3  # Initial Horizontal Position [m]
        wf1.mesh.xFin = 1e-3  # Final Horizontal Position [m]
        wf1.mesh.yStart = -1e-3  # Initial Vertical Position [m]
        wf1.mesh.yFin = 1e-3  # Final Vertical Position [m]

        wf1.partBeam = elecBeam  # append the undulatr beam we made to the wave front

        srwlpy.CalcElecFieldSR(wf1, 0, magFldCnt,
                               arPrecPar)  # calculate the electricfeild from the appended particle beam

        wfs.append(wf1)  # save the wf to the shelf

        time_finish = t.time()  # stop the timer

        # Print info about the wave field made
        print('\n\nMade WF with ' + str(i) + ' pixels,')
        print('Time taken: ' + str(time_finish - time_start) + ' seconds\n\n')

    d['wfs'] = wfs  # save the list of wave fields to the shelf

    d.close()  # close the shelf
