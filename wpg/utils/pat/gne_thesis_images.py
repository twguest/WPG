import numpy as np
from PIL import Image
import matplotlib.pyplot as plt


def imfile2array(fname):
    """input a filename (usally tiff) and return a numpy array of the data"""
    im = Image.open(fname)
    ar = np.array(im)
    return ar


def make_single_wf_plot(imfname, title, extent, plfname, label):
    plfname = r"C:\Users\patri\Documents\Uni\thesis\LaTeX\figs\exp_sim\\" + plfname
    plt.figure()
    plt.title(title)
    plt.imshow(imfile2array(imfname), cmap="jet", extent=extent)
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.colorbar()
    plt.savefig(plfname)


#
# make_single_wf_plot(r'C:\Users\patri\Documents\Uni\thesis\Synch\XFM_Sim\io\optics\base\norm\base000.tiff',
#           'Undulator Exit Wave Intensity\n 14m Upstream, 7.374keV',
#            [-1,1,-1,1],
#           'und_source_14m.png',
#           ('X Position [mm]','Y Position [mm]'))
#
#
# make_single_wf_plot(r'C:\Users\patri\Documents\Uni\thesis\Synch\XFM_Sim\io\optics\swap_kb_for_tl\norm\swap_kb_for_tl001.tiff',
#           'Wavefeild at WBS',
#            [-0.5,0.5,-0.5,0.5],
#           'wbs_14m.png',
#           ('X Position [mm]','Y Position [mm]'))
#

#
# make_single_wf_plot(r'C:\Users\patri\Documents\Uni\thesis\Synch\XFM_Sim\io\optics\swap_kb_for_tl\norm\swap_kb_for_tl002.tiff',
#           'Wavefeild 1m upstream of WBS',
#            [-0.5,0.5,-0.5,0.5],
#           'wbs_15m.png',
#           ('X Position [mm]','Y Position [mm]'))
#


# make_single_wf_plot(r'C:\Users\patri\Documents\Uni\thesis\Synch\XFM_Sim\io\optics\swap_kb_for_tl\norm\swap_kb_for_tl009.tiff',
#           'Wavefeild at ssa Slits',
#            [-0.5,0.5,-0.5,0.5],
#           'es_18m.png',
#           ('X Position [mm]','Y Position [mm]'))
#


# for i in range(20,36,4):
#     fname ='C:\\Users\\patri\\Documents\\Uni\\thesis\\Synch\\XFM_Sim\\io\\optics\\drift_past_focus\\norm\\drift_past_focus0'
#     fname += str(i)+'.tiff'
#     make_single_wf_plot(fname,'',
#                         [-132.219330484,130.8285099,-6.13762564551e1,5.9985435871e1],
#                         'drift_past_focus'+str(i)+'.png',
#                         ('X Position [um]','Y Position [um]'))
#
#


##Recon Comp


# vertical BIG PLOT
# plt.figure(dpi=100, figsize=(22.0/2,15.0/2))
#
#
# props = [30, 38, 40, 42, 44]
# prop_dist = ['35.46m','35.54m', '35.56m', '35.58m', '35.6m']
#
# im_margins = [[-120, 120, -60, 60],
#               [-60, 60, -30, 30],
#               [-35, 35, -17.5, 17.5],
#               [-35, 35, -17.5, 17.5],
#               [-35, 35, -17.5, 17.5]
#               ]
#
# plot_margins = [[-90, 90, -40, 40],
#               [-50, 50, -15, 15],
#               [-25, 25, -5, 5],
#               [-15, 15, -7.5, 7.5],
#               [-5, 5, -17.5, 17.5]
#               ]
# for i, prop in enumerate(props):
#     num = (2 * i + 1) + i
#     print(num)
#
#     plt.subplot(5, 3, num)
#     fname = 'C:\\Users\\patri\\Documents\\Uni\\thesis\\Synch\\XFM_Sim\\io\\optics\\drift_past_focus\\norm\\drift_past_focus0'
#     fname += str(prop) + '.tiff'
#     # print(fname)
#     im = imfile2array(fname)
#     print('Image Shape: ', im.shape)
#     plt.imshow(im, cmap='jet',
#                extent=[-132.219330484, 130.8285099, -6.13762564551e1, 5.9985435871e1],
#                aspect=0.5)
#     plt.title(prop_dist[i])
#     plt.colorbar()
#     plt.xlim(im_margins[i][0], im_margins[i][1])
#     plt.ylim(im_margins[i][2], im_margins[i][3])
#     if num == 7:
#         plt.ylabel('Y Position [um]')
#     elif num==13:
#         plt.xlabel('X Position [um]')
#
#     xpos = np.linspace(-132.219330484, 130.8285099, im.shape[1])
#     print(len(xpos))
#     ypos = np.linspace(-6.13762564551e1, 5.9985435871e1, im.shape[0])
#     print(len(ypos))
#
#     line_plot_x =im[int(im.shape[0] / 2),: ]
#     line_plot_y =im[:, int(im.shape[1] / 2)]
#
#
#     plt.subplot(5, 3, num + 1)
#     plt.plot(xpos, line_plot_x)#/ max(line_plot_x))
#     plt.xlim(plot_margins[i][0], plot_margins[i][1])
#
#     if num+1 ==8:
#         plt.ylabel('Intensity')
#
#     elif num+1==14:
#         plt.xlabel('X Position [um]')
#
#
#     plt.subplot(5, 3, num + 2)
#     plt.plot(ypos, line_plot_y)#/max(line_plot_y))
#     plt.xlim(plot_margins[i][2], plot_margins[i][3])
#
#     if num+2==9:
#         plt.ylabel('Intensity')
#
#     elif num+2==15:
#         plt.xlabel('Y Position [um]')
#
#
#
# plt.suptitle('Propagation Towards Focus')
# plt.subplots_adjust(left=0.07, right=0.97,bottom=0.08, top=0.93, hspace=0.5)
#
# plt.savefig(r'C:\Users\patri\Documents\Uni\thesis\LaTeX\figs\exp_sim\prop_arnd_foc.png')


plt.show()
