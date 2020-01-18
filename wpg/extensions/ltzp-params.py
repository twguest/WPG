
  # propagation parameter list & description:
    #[0]:  Auto-Resize (1) or not (0) Before propagation
    #[1]:  Auto-Resize (1) or not (0) After propagation
    #[2]:  Relative Precision for propagation with Auto-Resizing (1. is nominal)
    #[3]:  Allow (1) or not (0) for semi-analytical treatment of quadratic phase terms at propagation
    #[4]:  Do any Resizing on Fourier side, using FFT, (1) or not (0)
    #[5]:  Horizontal Range modification factor at Resizing (1. means no modification)
    #[6]:  Horizontal Resolution modification factor at Resizing
    #[7]:  Vertical Range modification factor at Resizing
    #[8]:  Vertical Resolution modification factor at Resizing
    #[9]:  Type of wavefront Shift before Resizing (not yet implemented)
    #[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
    #[11]: New Vertical wavefront Center position after Shift (not yet implemented)



wf0 created with 5000x5000 pixels



 el.append(SRWLOptD(5)) #length in m
    pp.append(propagationParameters(RangeX=2.0, RangeY=2.0, ResolutionX=2.0, ResolutionY=2.0))
PreOptics

SRWLOptD(5)   [0, 0, 1.0, 0, 0, 2.0, 2.0, 2.0, 2.0, 0, 0, 0]     propagationParameters(RangeX=2.0, RangeY=2.0, ResolutionX=2.0, ResolutionY=2.0))

             
Optical Element: Aperture / Obstacle
Prop. parameters = [0, 0, 1.0, 0, 0, 0.08, 30, 0.08, 30, 0, 0, 0]
        
Optical Element: Drift Space
Prop. parameters = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
        

slice   
Prop. parameters =       [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]

final drift =       [0, 0, 1.0, 0, 0, 10.0, 10.0, 0.1, 0.1, 0, 0, 0]