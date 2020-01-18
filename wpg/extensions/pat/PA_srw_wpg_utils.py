import os, sys
import shelve
import imageio
import numpy as np
sys.path.append(r'WPG')
sys.path.append(r'/opt/wpg/wpg')

OS_ENV = 'LINUX'

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



def PA_open_shelf(fname):

    if OS_ENV == 'LINUX':
        d = shelve.open('io/shelve_data/' + fname)  # Open the shelf doc for this instance
    else:
        d = shelve.open('io\\shelve_data\\'+fname)  # Open the shelf doc for this instance
    return d


# def PA_save_h5(wf, *args):
#     path = 'io'
#     if OS_ENV == 'LINUX':
#         for arg in args:
#             path = path + '/' + arg
#             if not os.path.exists(path):
#                 os.mkdir(path)
#         wf.store_hdf5(path + '/' + fname, )
#     else:
#         for arg in args:
#             path = path + '\\' + arg
#             if not os.path.exists(path):
#                 os.mkdir(path)
#         wf.store_hdf5(path+'\\'+fname,)






def PA_save_tiff(wf, fname, path = []):

    if OS_ENV == 'LINUX':

        save_path = 'io'
        for direc in path:
            save_path = save_path + '/' + direc
            if not os.path.exists(save_path):
                os.mkdir(save_path)
        
        try: 
            imageio.imwrite(save_path+'/'+fname, wf)
        except:
            pass
        else:
            print('Saved file to ' + save_path+'/'+fname)

    else:
        save_path = 'io'
        for direc in path:
            save_path = save_path + '\\' + direc
            if not os.path.exists(save_path):
                os.mkdir(save_path)


        try:
            imageio.imwrite(save_path +'\\'+fname, wf)
        except:
            pass
        else:
            print('Saved file to ' + save_path+'\\'+fname)


def PA_wf_data_to_csv(wf, fname, path = [], file_mode = 'w'):
    if OS_ENV == 'LINUX':
        save_path = 'io'
        for direc in path:
            save_path = save_path + '/' + direc
            if not os.path.exists(save_path):
                os.mkdir(save_path)


        start_x = str(wf._srwl_wf.mesh.xStart)
        finish_x = str(wf._srwl_wf.mesh.xFin)

        start_y = str(wf._srwl_wf.mesh.yStart)
        finish_y = str(wf._srwl_wf.mesh.yFin)

        npx_x = str(wf._srwl_wf.mesh.nx)
        npx_y = str(wf._srwl_wf.mesh.ny)

        px_size_x = str( (abs(wf._srwl_wf.mesh.xStart) + abs(wf._srwl_wf.mesh.xFin)) / abs(wf._srwl_wf.mesh.nx))
        px_size_y = str((abs(wf._srwl_wf.mesh.yStart) + abs(wf._srwl_wf.mesh.yFin)) / abs(wf._srwl_wf.mesh.ny))

        row = (start_x, finish_x, start_y, finish_y, npx_x, npx_y, px_size_x, px_size_y)
        content = ','.join(row)


        wf_file = open(save_path+'/'+fname, file_mode)

        if file_mode=='w':
            wf_file.write('Start_X,Finish_x,Start_Y,Finish_Y,Num_Px_X,Num_Px_y,Px_Size_X,Px_Size_X\n')

        wf_file.write(content+'\n')
        wf_file.close()



    else:
        save_path = 'io'
        for direc in path:
            save_path = save_path + '\\' + direc
            if not os.path.exists(save_path):
                os.mkdir(save_path)


        start_x = str(wf._srwl_wf.mesh.xStart)
        finish_x = str(wf._srwl_wf.mesh.xFin)

        start_y = str(wf._srwl_wf.mesh.yStart)
        finish_y = str(wf._srwl_wf.mesh.yFin)

        npx_x = str(wf._srwl_wf.mesh.nx)
        npx_y = str(wf._srwl_wf.mesh.ny)

        px_size_x = str( (abs(wf._srwl_wf.mesh.xStart) + abs(wf._srwl_wf.mesh.xFin)) / abs(wf._srwl_wf.mesh.nx))
        px_size_y = str((abs(wf._srwl_wf.mesh.yStart) + abs(wf._srwl_wf.mesh.yFin)) / abs(wf._srwl_wf.mesh.ny))

        row = (start_x, finish_x, start_y, finish_y, npx_x, npx_y, px_size_x, px_size_y)
        content = ','.join(row)


        wf_file = open(save_path+'\\'+fname, file_mode)

        if file_mode=='w':
            wf_file.write('Start_X,Finish_x,Start_Y,Finish_Y,Num_Px_X,Num_Px_y,Px_Size_X,Px_Size_X\n')

        wf_file.write(content+'\n')
        wf_file.close()









def PA_concat_list(*args):

    list1 = args[0]

    concatenated_list = []
    for arg in args:

        for arg_elm in arg:

            concatenated_list.append(arg_elm)

    return concatenated_list





