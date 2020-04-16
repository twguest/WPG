import matplotlib

import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import minmax_scale as norm


from wpg.wavefront import Wavefront
from wpg.generators import build_gauss_wavefront_xy

J2EV = 6.24150934e18


def plotWavefront(wfr, item="both", save=False, savedir="", gsm=False, lims=None):

    if gsm == True:
        ii = wfr.II[:, :, 0]
    else:
        ii = wfr.get_intensity()[:, :, 0]

    II = ii * wfr.params.photonEnergy / J2EV

    phase = wfr.get_phase()[:, :, 0]
    phase = (phase + np.pi) % (2 * np.pi) - np.pi
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wfr)

    dx = (xmax - xmin) / (nx - 1)
    dy = (ymax - ymin) / (ny - 1)

    print("stepX, stepY [um]:", dx * 1e6, dy * 1e6, "\n")

    [x1, x2, y1, y2] = wfr.get_limits()
    matplotlib.rc("xtick", labelsize=28)
    matplotlib.rc("ytick", labelsize=28)
    if item == "both":

        fig = plt.figure(figsize=(21, 6))

        ax1 = fig.add_subplot(121)
        ax1.yaxis.set_major_formatter(FormatStrFormatter("%d"))
        ax1.xaxis.set_major_formatter(FormatStrFormatter("%d"))

        im1 = ax1.imshow(
            II, extent=[x1 * 1e6, x2 * 1e6, y1 * 1e6, y2 * 1e6], cmap=cm.gray
        )

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        cb1 = fig.colorbar(im1, cax=cax, orientation="vertical")
        cb1.ax.set_title("I (W/$mm^2$)", fontsize=20)
        ax1.set_xlabel("x ($\mu$m)", fontsize=20)
        ax1.set_ylabel("y ($\mu$m)", fontsize=20)
        ax1.axis("tight")

        ax2 = fig.add_subplot(122)
        im2 = ax2.imshow(
            phase, extent=[x1 * 1e6, x2 * 1e6, y1 * 1e6, y2 * 1e6], cmap=cm.phase
        )
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im2, cax=cax, orientation="vertical")

        cb2 = fig.colorbar(
            im2, cax=cax, orientation="vertical"
        )  # , ticks=[-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
        cb2.ax.set_ylabel("Phase (rad.)", fontsize=20)

        # cb2.ax.set_yticklabels(['$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'])

        ax2.set_xlabel("x ($\mu$m)", fontsize=20)
        ax2.set_ylabel("y ($\mu$m)", fontsize=20)
        ax2.axis("tight")
        ax2.yaxis.set_major_formatter(FormatStrFormatter("%d"))
        ax2.xaxis.set_major_formatter(FormatStrFormatter("%d"))

        fig.tight_layout()
    elif item == "phase":

        fig = plt.figure(figsize=(18, 18))

        ax2 = fig.add_subplot(111)
        im2 = ax2.imshow(
            phase, extent=[x1 * 1e6, x2 * 1e6, y1 * 1e6, y2 * 1e6], cmap=cm.phase
        )
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im2, cax=cax, orientation="vertical")

        cb2 = fig.colorbar(
            im2,
            cax=cax,
            orientation="vertical",
            ticks=[-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
        )
        cb2.ax.set_ylabel("Phase (rad.)", fontsize=24)

        cb2.ax.set_yticklabels(["$-\pi$", "$-\pi/2$", "0", "$\pi/2$", "$\pi$"])
        ax2.set_xlabel("x ($\mu$m)", fontsize=24)
        ax2.set_ylabel("y ($\mu$m)", fontsize=24)
        ax2.axis("tight")
        ax2.yaxis.set_major_formatter(FormatStrFormatter("%d"))
        ax2.xaxis.set_major_formatter(FormatStrFormatter("%d"))
        plt.subplots_adjust(top=0.9)

        if lims is not None:
            ax2.set_xlim(lims[0])
            ax2.set_ylim(lims[1])

    elif item == "intensity":

        fig = plt.figure(figsize=(18, 18))
        ax1 = fig.add_subplot(111)
        ax1.yaxis.set_major_formatter(FormatStrFormatter("%d"))
        ax1.xaxis.set_major_formatter(FormatStrFormatter("%d"))

        im1 = ax1.imshow(
            II, extent=[x1 * 1e6, x2 * 1e6, y1 * 1e6, y2 * 1e6], cmap=cm.gray
        )  # "CMRmap")

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        cb1 = fig.colorbar(im1, cax=cax, orientation="vertical")
        cb1.ax.set_title("I (W/$mm^2$)", fontsize=24)
        cb1.ax.yaxis.set_major_formatter(FormatStrFormatter("%.2e"))
        ax1.set_xlabel("x ($\mu$m)", fontsize=24)
        ax1.set_ylabel("y ($\mu$m)", fontsize=24)
        ax1.axis("tight")

        if lims is not None:
            ax1.set_xlim(lims[0])
            ax1.set_ylim(lims[1])
        plt.subplots_adjust(top=0.9, right=0.9)

    elif item == "image":
        matplotlib.rc("xtick", labelsize=20)
        matplotlib.rc("ytick", labelsize=20)
        fig = plt.figure(figsize=(9, 9))
        ax1 = fig.add_subplot(111)
        ax1.yaxis.set_major_formatter(FormatStrFormatter("%d"))
        ax1.xaxis.set_major_formatter(FormatStrFormatter("%d"))

        im1 = ax1.imshow(
            II, extent=[x1 * 1e6, x2 * 1e6, y1 * 1e6, y2 * 1e6], cmap=cm.gray
        )  # "CMRmap")
        ax1.axis("tight")
        plt.axis("off")

        if lims is not None:
            ax1.set_xlim(lims[0])
            ax1.set_ylim(lims[1])
        plt.subplots_adjust(top=0.9)

    plt.show()
    if save == True:
        fig.savefig(savedir + item + ".png", format="png")

    return fig


def plot3d(wfr, nplots=1):

    fig = plt.figure(figsize=[12, 12], dpi=100)
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    # ax._axis3don = False
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    
    import os

    os.path.join(".")
    ax.xaxis._axinfo["grid"]["color"] = (1, 1, 1, 0)
    ax.yaxis._axinfo["grid"]["color"] = (1, 1, 1, 0)
    # ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.w_xaxis.line.set_lw(0.0)
    ax.w_yaxis.line.set_lw(0.0)

    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("Irradiance (W/$mm^2$)")
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wfr)
    xa = np.linspace(xmin, xmax, nx)
    ya = np.linspace(ymin, ymax, ny)

    # Make data.
    ii = wfr.II[:, :, 0]
    phase = wfr.get_phase()[:, :, 0]
    xv, yv = np.meshgrid(xa, ya)

    surf = ax.plot_surface(
        xv,
        yv,
        ii,
        cmap="hot",
        linewidth=0,
        antialiased=False,
        facecolors=cm.hsv(norm(phase)),
    )

    normal = matplotlib.colors.Normalize(vmin=-2 * np.pi, vmax=2 * np.pi)

    m = cm.ScalarMappable(cmap=plt.cm.hsv, norm=normal)
    m.set_array([])

    cb = fig.colorbar(
        m, shrink=0.5, aspect=5, ticks=[-2 * np.pi, -np.pi, 0, np.pi, 2 * np.pi]
    )
    cb.ax.set_yticklabels(["-2$\pi$", "-$\pi$", "0", "$\pi$", "2$\pi$"])


def plotPhaseSpace(phasespace, title="", save=False, outdir="", fileid=""):
    x, x_, y, y_ = [], [], [], []

    for a in phasespace:
        x.append(a[1] * 1e6)
        x_.append(a[2] * 1e6)
        y.append(a[3] * 1e6)
        y_.append(a[4] * 1e6)

    ax1 = sns.jointplot(
        x, y, kind="hex", label=fileid, xlim=[-2000, 2000], ylim=[-20, 20]
    )
    ax1.set_axis_labels(
        "Horizontal Position ($\mu$m)", "Veritcal Position ($\mu$m)", fontsize=22
    )

    ax1.ax_joint.tick_params(axis="both", which="major", labelsize=16)

    ax2 = sns.jointplot(x_, y_, kind="kde", xlim=[-100, 100], ylim=[-6, 6])
    ax2.set_axis_labels(
        "Horizontal Trajectory ($\mu$rad)",
        "Vertical Trajectory ($\mu$rad)",
        fontsize=22,
    )
    ax2.ax_joint.tick_params(axis="both", which="major", labelsize=16)
    ax2.ax_joint.ticklabel_format(axis="both", style="sci")

    ax3 = sns.jointplot(x, x_, kind="reg", xlim=[-2000, 2000], ylim=[-100, 100])
    ax3.set_axis_labels(
        "Horizontal Position ($\mu$m)", "Horizontal Trajectory ($\mu$rad)", fontsize=22
    )
    ax3.ax_joint.tick_params(axis="both", which="major", labelsize=16)
    ax3.ax_joint.ticklabel_format(axis="both", style="sci")

    ax4 = sns.jointplot(y, y_, color="red", kind="reg", xlim=[-20, 20], ylim=[-6, 6])
    ax4.set_axis_labels(
        "Vertical Position ($\mu$m)", "Vertical Trajectory ($\mu$rad)", fontsize=22
    )
    ax4.ax_joint.tick_params(axis="both", which="major", labelsize=16)
    ax4.ax_joint.ticklabel_format(axis="both", style="sci")

    if save == True:
        ax1.savefig(outdir + "realspace" + str(phasespace.shape[0]) + ".png")
        ax2.savefig(outdir + "phasespace" + str(phasespace.shape[0]) + ".png")
        ax3.savefig(outdir + "x_ellipse" + str(phasespace.shape[0]) + ".png")
        ax4.savefig(outdir + "y_ellipse" + str(phasespace.shape[0]) + ".png")



flatui = ["#1663a6", "#f89406", "#14b258", "#bf110b"]
palettes = {"flatui": flatui}


def show_palettes():

    sns.palplot(palettes["flatui"])
    plt.title("Flat UI")


def rmsPlot(wfrs, axis="x", save=False, savedir=""):

    if type(wfrs) is not list:
        assert "input should be a list of wavefronts"

    sns.set()
    cmap = plt.get_cmap("tab10")
    cmap = [cmap(c) for c in range(100)]
    print(len(cmap))
    matplotlib.rc("xtick", labelsize=14)
    matplotlib.rc("ytick", labelsize=14)
    fig = plt.figure(figsize=(21, 6))
    ax = fig.add_subplot(111)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))

    itr = 0
    leg = []
    for wfr in wfrs:
        wfr.type = "se"
        wfr.electric_field()
        wfr.intensity()

        ii = wfr.II[:, :, 0]
        II = ii * wfr.params.photonEnergy / J2EV
        [x1, x2, y1, y2] = wfr.get_limits()

        axis = np.linspace(x1 * 1e6, x2 * 1e6, II.shape[1])
        plt.subplots_adjust(top=0.9)

        ax.plot(
            axis,
            II[II.shape[0] // 2, :],
            color=cmap[itr],
            label="$\sigma_x$ = " + "{:.2E}".format(get_rms(wfr)[4]) + " $\mu$m",
        )
        ax.plot(
            axis,
            II[:, II.shape[1] // 2,],
            color=cmap[itr],
            label="$\sigma_y$ = " + "{:.2E}".format(get_rms(wfr)[3]) + " $\mu$m",
        )
        itr += 1

    ax.legend(fontsize=16)
    ax.set_xlabel("x ($\mu$m)", fontsize=22)
    ax.set_ylabel("Intensity (W/$mm^2$)", fontsize=22)

    plt.show()

    if save == True:
        fig.savefig(savedir + "rsmplot")


def lineProf(wfr, show=True, lims=False, gsm=False, save=False, savedir=""):

    sns.set()

    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wfr)
    ax_x = np.linspace(xmin * 1e6, xmax * 1e6, nx)
    ax_y = np.linspace(ymin * 1e6, ymax * 1e6, ny)

    if gsm == False:
        II = wfr.get_intensity()
        ph = wfr.get_phase()
    elif gsm == True:
        II = wfr.II
        ph = wfr.get_phase()

    fig = plt.figure(figsize=[12, 8])

    Ix = II[:, II.shape[1] // 2, :]
    Iy = II[II.shape[0] // 2, :, :]

    px = ph[:, ph.shape[1] // 2, :]
    py = ph[ph.shape[0] // 2, :, :]

    matplotlib.rc("xtick", labelsize=18)
    matplotlib.rc("ytick", labelsize=18)

    ax1 = fig.add_subplot(111)
    ax1.plot(ax_x, Ix, color="blue", label="Horizontal Profile")
    ax1.set_ylabel("Intensity (W/$mm^2$)", fontsize=22)

    ax1.plot(ax_y, Iy, color="red", label="Vertical Profile")
    ax1.set_xlabel("Position ($\mu m$)", fontsize=22)

    plt.legend(fontsize=22)

    if show == True:
        plt.show()

def plotWavefront(wfr):
    
    fig = plt.figure(figsize = [12,8])
    
    ii = wfr.get_intensity()[:,:,0]
    ph = wfr.get_phase()[:,:,0]
    
    plt.imshow(ph*ii)
if __name__ == '__main__':
    
    ### CONSTRUCT ARB. WFR
    wfr = Wavefront(build_gauss_wavefront_xy(256, 256, 9.2, -5e-07, 5e-07, -5e-07, 5e-07, 1e-07, 1e-07, 1e-09))
    plotWavefront(wfr)