from extensions.twg.current import sxri_bl
#from gvr_utils import plotWavefront
#from scipy import ndimage



# =============================================================================
# #### PROPOGATE FULLY COHERENT
beamline = sxri_bl() 
beamline.und_run(read =  False) 
beamline.setup_OE()
beamline.setup_beamline(VERBOSE = False)
beamline.propagate_beamline(read = False)
 
beamline.setup_endstation()
beamline.propagate_endstation()
 
# =============================================================================

# =============================================================================
# 
# #### PROPOGATE PARTIALLY COHERENT
# beamline = sxri_bl() 
# beamline.setup_gsm(p_mu = .75)
# beamline.setup_OE()
# beamline.setup_beamline(VERBOSE = False)
# beamline.propagate_modes()
# beamline.collapse_modes()
# beamline.wfr.coherence()
# 
# ### STILL NEED TO SETUP AND PROPOGATE ENDSTATION
# 
# =============================================================================
