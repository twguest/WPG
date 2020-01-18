import wpg.srwlpy as srwl
from wpg.srwlib import SRWLWfr, SRWLRadMesh, SRWLGsnBm, SRWLStokes

import random
import numpy as np
from array import *
from math import *
from copy import *

from extensions.twpg_wavefront import Wavefront
from extensions.run.metrology_bl import sxri_bl 
 
from wpg.srwlib import SRWLRadMesh as Mesh
from joblib import load

def srwl_wfr_emit_prop_multi_e(output, _e_beam, _mag, _mesh, ne): #OC07092016 (added _wr)
    """
    Calculate Stokes Parameters of Emitted (and Propagated, if beamline is defined) Partially-Coherent SR
    :param _e_beam: Finite-Emittance e-beam (SRWLPartBeam type)
    :param _mag: Magnetic Field container (magFldCnt type)
    :param _mesh: mesh vs photon energy, horizontal and vertical positions (SRWLRadMesh type) on which initial SR should be calculated
    :param _sr_meth: SR Electric Field calculation method to be used (0- "manual", 1- "auto-undulator", 2- "auto-wiggler")
    :param _sr_rel_prec: relative precision for SR Electric Field calculation (usually 0.01 is OK, the smaller the more accurate)
    :param _n_part_tot: total number of "macro-electrons" to be used in the calculation
    :param _n_part_avg_proc: number of "macro-electrons" to be used in calculation at each "slave" before sending Stokes data to "master" (effective if the calculation is run via MPI)
    :param _n_save_per: periodicity of saving intermediate average Stokes data to file by master process
    :param _file_path: path to file for saving intermediate average Stokes data by master process
    :param _sr_samp_fact: oversampling factor for calculating of initial wavefront for subsequent propagation (effective if >0)
    :param _opt_bl: optical beamline (container) to propagate the radiation through (SRWLOptC type)
    :param _pres_ang: switch specifying presentation of the resulting Stokes parameters: coordinate (0) or angular (1)
    :param _char: radiation characteristic to calculate:
        0- Intensity (s0);
        1- Four Stokes components;
        2- Mutual Intensity Cut vs X;
        3- Mutual Intensity Cut vs Y;
        4- Mutual Intensity Cut vs X & Y;
        10- Flux
    :param _x0: horizontal center position for mutual intensity calculation
    :param _y0: vertical center position for mutual intensity calculation
    :param _e_ph_integ: integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from _mesh
    :param _rand_meth: method for generation of pseudo-random numbers for e-beam phase-space integration:
        1- standard pseudo-random number generator
        2- Halton sequences
        3- LPtau sequences (to be implemented)
    :param _tryToUseMPI: switch specifying whether MPI should be attempted to be used
    :param _wr: initial wavefront radius [m] to assume at wavefront propagation (is taken into account if != 0)
    """
    
    outpath = r"./extensions/out/me/"

    _sr_meth = 1,
    _sr_rel_prec = 0.1,
    _n_part_tot = 10,
    _n_part_avg_proc = 10,
    _n_save_per = 10,
    _sr_samp_fact = 10,
    _file_path = r"out.ascii"

    log = []
    
    wfr = SRWLWfr() #Wavefronts to be used in each process
    wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
    
    wfr.mesh.set_from_other(_mesh)    
    wfr.partBeam = deepcopy(_e_beam)
    arPrecParSR = [_sr_meth, _sr_rel_prec, 0, 0, 1000, 1, 0] #to add npTraj, useTermin ([4], [5]) terms as input parameters

    meshRes = SRWLRadMesh(_mesh.eStart, _mesh.eFin, _mesh.ne, _mesh.xStart, _mesh.xFin, _mesh.nx, _mesh.yStart, _mesh.yFin, _mesh.ny, _mesh.zStart) #to ensure correct final mesh if _opt_bl==None


    sigX = 389e-06
    sigXp = 19.7e-06
    sigY = 40.8e-06
    sigYp = 8.0e-06

    randAr = np.zeros(5)
    
    wfr.partBeam.arStatMom2[0] = 0 #<(x-x0)^2>
    wfr.partBeam.arStatMom2[1] = 0 #<(x-x0)*(xp-xp0)>
    wfr.partBeam.arStatMom2[2] = 0 #<(xp-xp0)^2>
    wfr.partBeam.arStatMom2[3] = 0 #<(y-y0)^2>
    wfr.partBeam.arStatMom2[4] = 0 #<(y-y0)*(yp-yp0)>
    wfr.partBeam.arStatMom2[5] = 0 #<(yp-yp0)^2>
    wfr.partBeam.arStatMom2[10] = 0
    for i in range(ne): #loop over macro-electrons
# =============================================================================
#         
#         for ir in range(5): #to expend to 6D eventually
#             randAr[ir] = random.gauss(0, 0)
#         wfr.partBeam.partStatMom1.x = randAr[1]*sigX
#         wfr.partBeam.partStatMom1.xp = randAr[1]*sigXp
#         wfr.partBeam.partStatMom1.y = randAr[2]*sigY
#         wfr.partBeam.partStatMom1.yp = randAr[3]*sigYp
#         
# =============================================================================
        wfr.presCA = 0 #presentation/domain: 0- coordinates, 1- angles
        wfr.presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time
        
          
        print('i=', i, 'Electron Coord.: x=', wfr.partBeam.partStatMom1.x, 'x\'=', wfr.partBeam.partStatMom1.xp, 'y=', wfr.partBeam.partStatMom1.y, 'y\'=', wfr.partBeam.partStatMom1.yp, 'E=',  wfr.partBeam.partStatMom1.gamma*0.51099890221e-03)
          
        data = [i, wfr.partBeam.partStatMom1.x, wfr.partBeam.partStatMom1.xp, wfr.partBeam.partStatMom1.y, wfr.partBeam.partStatMom1.yp ]
             
        srwl.CalcElecFieldSR(wfr, 0, _mag, arPrecParSR) #calculate Electric Field emitted by current electron
        cwfr = Wavefront(wfr)
        output += cwfr.data.arrEhor


    return "Finished Multi Electron Propagation"


if __name__ == "__main__":
    wfr = Wavefront()
    wfr.load_hdf5(r"./extensions/out/se/wfr_mode_se.hdf5")
    
    output = np.memmap(r"./extensions/out/tmp/memmap_me", dtype=wfr.data.arrEhor.dtype,
                   shape= wfr.data.arrEhor.shape, mode='w+')
    
    beamline = sxri_bl()
    beamline.setup_OE()
    _mesh = Mesh(_eStart = beamline.params.w_e, _eFin = beamline.params.w_e, _ne=1,
             _nx=2048, _ny=2048,
             _xStart = -20e-04, _xFin = 20e-04,
             _yStart = -20e-04, _yFin = 20e-04,
             _zStart = 1)
    beamline.setup_src()
    srwl_wfr_emit_prop_multi_e(output, beamline.elecBeam, beamline.magFldCnt,
          _mesh, 1)

    final = load(r"./extensions/out/tmp/memmap_me")