# build_optics.py
# using the premade optical elements and propagation perameters, build the optical beamline and
# propagate the wavefronts from the source to the sample plane

from PA_srw_wpg_utils import *
from extensions.beamlineExt import bl
from extensions.gvrutils import *
from wpg.srwl_uti_smp import srwl_opt_setup_transm_from_file
from wpg.optical_elements import *
import os
import pickle


def make_beamline():
    optics_shelf = PA_open_shelf("optics")

    op = optics_shelf["op"]  # load optical elements defined in make_optics5

    bls = []  # Init list of beamlines, one for each optical element

    for i, opt in enumerate(
        op.arOpt
    ):  # loop through each optical element (-2 to chop outlast two elements)
        bl = beamline.Beamline()  # make a seperate beamline
        bl.append(
            opt, op.arProp[i]
        )  # append just the single element to the new beamline
        bls.append(
            bl
        )  # append the beamline with a single element to the list of beamlines

    optics_shelf.close()

    return bls


def load_wf():
    source_shelf = PA_open_shelf("source")

    wfs = source_shelf["wfs"]  # load wavefronts from make_source

    wf1 = wfs[3]  # select a single wavefront (check 'pixels' list in make_und_source)

    wf_wpg = wavefront.Wavefront(wf1)

    source_shelf.close()

    return wf_wpg


def propagate_wavefront(wf_wpg, bls, fname="x"):

    wf_wpg_inten = wf_wpg.get_intensity(slice_number=0)  # get intensity

    wf_wpg_inten -= np.min(wf_wpg_inten)

    # PA_wf_data_to_csv(wf_wpg, fname+'_array_data.csv', ['optics', fname] )
    PA_save_tiff(wf_wpg_inten, fname + "000.tiff", ["optics", fname, "norm"])
    # PA_save_tiff(np.log10(wf_wpg_inten), fname + 'log10_000.tiff', ['optics', fname, 'log'])

    # Propagate the wave field through every optical element
    for i, bl in enumerate(bls):
        print(bl)
        bl.propagate(wf_wpg)  # propagate though the current optical element.

        wf_wpg_inten = wf_wpg.get_intensity(slice_number=0)  # get intensity
        wf_wpg_inten -= np.min(wf_wpg_inten)

        if (i + 1) / 10 == 0:
            fname_num = "00" + str(i + 1)
        else:
            fname_num = "0" + str(i + 1)

        # PA_wf_data_to_csv(wf_wpg, fname + '_array_data.csv', ['optics', fname], file_mode='a')
        PA_save_tiff(
            wf_wpg_inten, fname + fname_num + ".tiff", ["optics", fname, "norm"]
        )
        # PA_save_tiff(np.log10(wf_wpg_inten), fname + 'log10_'+ fname_num+'.tiff', ['optics', fname, 'log'])


OBJECTPLANEWFRFILE = "objPlane.h5"

if not (os.path.exists(OBJECTPLANEWFRFILE)):

    bls = make_beamline()

    wf_wpg = load_wf()

    propagate_wavefront(wf_wpg, bls, "drift_past_focus")

    plotWavefront(wf_wpg, "At object plane")

    imageio.imwrite("objPlane.tif", wf_wpg.get_intensity(polarization="horizontal"))
    wf_wpg.store_hdf5("objPlane.h5")


else:
    wf_wpg = Wavefront()
    wf_wpg.load_hdf5(OBJECTPLANEWFRFILE)


[nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wf_wpg)
dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
rx = xmax - xmin
ry = ymax - ymin

print(
    "WF at object Plane: Number={}, {}, Range={}, {}, Step={}, {}".format(
        nx, ny, rx, ry, dx, dy
    )
)


tr = pickle.load(open("io/sample/tr.p", "rb"))

# blSample = bl(SRWLOptC(_arOpt=[tr,],
#                       _arProp=[ propagationParameters(),]),
#                       description = 'Sample Transmission')
# blSample.propagate(wf_wpg)

srwl.PropagElecField(
    wf_wpg._srwl_wf, SRWLOptC(_arOpt=[tr,], _arProp=[propagationParameters(),])
)


"""
sample_bl = beamline.Beamline()

sample_bl.append(tr, [0, 0, 1.0, 0, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
sample_bl.propagate(wf_wpg)
"""


imageio.imwrite("esw.tif", wf_wpg.get_intensity(polarization="horizontal"))


# drift to detector
blDriftDetector = bl(
    SRWLOptC(
        _arOpt=[srwlib.SRWLOptD(3.8)],
        _arProp=[
            propagationParameters(
                SemiAnalyt=1, RangeX=10.0, ResolutionX=0.1, RangeY=10.0, ResolutionY=0.1
            )
        ],
    ),
    description="Drift to Detector",
)
blDriftDetector.propagate(wf_wpg)


# plotWavefront(wf_wpg,'At detector before resize/rescale')
imageio.imwrite("detPlane.tif", wf_wpg.get_intensity(polarization="horizontal"))

Nx, Ny = 512, 512
[nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wf_wpg)
dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
rx = xmax - xmin
ry = ymax - ymin

# Range and resolution that is required for next step:
Dx, Dy = 75e-6, 75e-6
Rx, Ry = Nx * Dx, Ny * Dy


print(
    "Before Rescale: Number={}, {}, Range={}, {}, Step={}, {}".format(
        nx, ny, rx, ry, dx, dy
    )
)
elResc = bl(
    SRWLOptC(
        _arOpt=[SRWLOptD(0)],
        _arProp=[propagationParameters(ResolutionX=dx / Dx, ResolutionY=dy / Dy)],
    ),
    description="Rescaling to match optic",
)
metrics = elResc.propagate(wfr=wf_wpg, describe=False)
# writeMetrics(strOutPath+'/'+p.LABEL+'-summary.csv', metrics, section.description)


[nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wf_wpg)
dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
rx = xmax - xmin
ry = ymax - ymin
print(
    "After Rescale: Number={}, {}, Range={}, {}, Step={}, {}".format(
        nx, ny, rx, ry, dx, dy
    )
)


print(
    "Before Resize: Number={}, {}, Range={}, {}, Step={}, {}".format(
        nx, ny, rx, ry, dx, dy
    )
)
elResc = bl(
    SRWLOptC(
        _arOpt=[SRWLOptD(0)],
        _arProp=[propagationParameters(RangeX=Rx / rx, RangeY=Ry / ry)],
    ),
    description="Resizing to match optic",
)
metrics = elResc.propagate(wfr=wf_wpg, describe=False)
# writeMetrics(strOutPath+'/'+p.LABEL+'-summary.csv', metrics, section.description)


[nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wf_wpg)
dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
rx = xmax - xmin
ry = ymax - ymin
print(
    "After Resize: Number={}, {}, Range={}, {}, Step={}, {}".format(
        nx, ny, rx, ry, dx, dy
    )
)


# plotWavefront(wf_wpg,'At detector after resize & rescale')
wf_wpg.store_hdf5("detPlane.h5")
imageio.imwrite("detPlaneR.tif", wf_wpg.get_intensity(polarization="horizontal"))
