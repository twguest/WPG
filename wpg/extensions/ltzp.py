#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 09:57:28 2019

@author: gvanriessen
"""
import numpy as np
import math
import matplotlib.pyplot as plt

import pymesh
import meshio



# generate primitive shapes of LTZP oriented along the y axis that is perpendicular to the linear zones.




def rz(n,f,wavelength):
    return math.sqrt( (n*f*wavelength) + ((n*wavelength)**2.0)/4.0)



#def extrudePolygon2D(polygon):
#    """
#    Define polygon as points eg, [[-0.5, -0.3], [0.5, -0.3], [0.0, 0.5]]
#    """
#
#    p = pygalmesh.Polygon2D(polygon)
#    domain = pygalmesh.Extrude(p, [0.0, 0.3, 1.0])
#    pygalmesh.generate_mesh(
#        domain, "out.mesh", cell_size=0.1, edge_size=0.1, verbose=False
#    )
#
#    mesh = meshio.read("out.mesh")
#
#    #vol = sum(compute_volumes(mesh.points, mesh.cells["tetra"]))
#
#    return mesh

def extrudePolygon2D(polygon,thickness=100e-9):
    """
    Define polygon as list of points, either x-y 2-tuples 
    (e.g. [[-0.5, -0.3], [0.5, -0.3], [0.0, 0.5]]), or 
    x-y-z 3-tuples [[-0.5, -0.3,0], [0.5, -0.3,0], [0.0, 0.5,0]]),
    

    """
    
    if len(polygon[0]) == 2:
        poly3 = []
        for pt in polygon:
            poly3.append([pt[0],pt[1],0])
        polygon = poly3     
    print (polygon)
    
    
    import pygmsh
    geom = pygmsh.built_in.Geometry()
    
    poly = geom.add_polygon(polygon,
                            lcar=0.05)
    axis = [0, 0, thickness]
    
    geom.extrude(
        poly,
        translation_axis=axis,
        )
    
    #points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom, geo_filename="polygon.geo")

    mesh = pygmsh.generate_mesh(geom, geo_filename="polygon.geo")



   # from helpers import compute_volume
   # print ('Volume %s' % str(compute_volume(points, cells)) )


    #return points, cells, point_data, cell_data, field_data
    return mesh



def primitivePoly(n,f,wavelength,NZones,resolution,
                  width, 
                  shape = 'quadratic' # or 'kineform', or 'triangle'
                  ):

    zones = list(range(0,NZones))

    # zone widths
    zw = [rz(zone,f,wavelength) for zone in zones]

   
    #, x,y coordinate of points along edges of 2d shape
    x = np.linspace(0,width,width/resolution)

    if (shape=='quadratic'):
        # quadratic
        y = np.asarray([  xi**2 for xi in x[:int(len(x)/2)] ] + [ (width-xi)**2 for xi in x[int(len(x)/2):] ])
        y = zw[n]*(y/max(y))
    elif (shape=='kineform'):
        y = np.asarray([  xi**2/(2*wavelength*f) for xi in x[:int(len(x)/2)] ] + [ (width-xi)**2/(2*wavelength*f) for xi in x[int(len(x)/2):] ])
        y = zw[n]*(y/max(y))
    elif (shape=='triangle'):
        x,y = np.asarray([0.,width/2.0,width]),np.asarray([0.,zw[n],0.])
        
    y = y + np.sum(zw[1:n])    # translate along y

    #return list(zip(x,y))
    return [[p1, p2] for idx1, p1 in enumerate(x) for idx2, p2 in enumerate(y) if idx1==idx2]


# =============================================================================
# pymesh slicing
# 
# =============================================================================

def mergeMeshes(input_meshes, outname='merged.mesh'):

	out = pymesh.merge_meshes(input_meshes)
	if (outname != None):
		pymesh.save_mesh(outname, out)

	return out   

def sliceMesh(mesh,direction, N):

      

    slices = pymesh.slice_mesh(mesh, direction, N)
    
    for m,i in zip(slices,range(N)):

    
        pymesh.save_mesh('mesh_slice%d.stl' % i, m, ascii=True);
    

def testSinglePrimitive():

    AngleD = 18 # [degrees]
    Height=300.e-9
    Angle = math.radians(AngleD)
    proj = (Height * math.cos(Angle))+ Height;  # projected thickness at some angle,
    outerMostZone =  200e-9 # m
    
    Energy = 0.6   # keV
    wavelength = 2.06e-9  # m
    f = 2.904e-3   # focal length, 2.904e-3 for 30um ZP  % 4.84e-3 for 50 um ZP
    #outerMostZone = 200e-9  # %outer most zone, [m]
    width = 0.4e-6 # [m]
    NZones =  round( wavelength * f / 4.0 / outerMostZone**2.0)  # number of zones
    
    resolution = 1.e-9  # disrtance between points along x-axis
    
    fig = plt.figure()

    #for sh in ['quadratic','kineform', 'triangle']:
    for sh in ['quadratic']:
        n=1
        primitive = primitivePoly(n,f,wavelength,NZones,resolution,width, shape = sh)
        x,y=map(list, zip(*primitive))
        x,y=np.asarray(x), np.asarray(y)
        
        #mesh = extrudePolygon2D(primitive)
        points, cells, point_data, cell_data, field_data = extrudePolygon2D(primitive)
        meshio.write_points_cells('primitive.mesh', 
        	                      points, 
        	                      cells, 
        	                      cell_data=cell_data,
        	                      point_data=point_data,
        	                      field_data=field_data)
        #meshio.write_points_cells('primitive.stl', points, cells, cell_data=cell_data)
  
        plt.scatter(x*1e9, y*1e9,  marker='o')
    
    plt.show()
    

def testZoneOfPrimitives():
    
    
    AngleD = 18 # [degrees]
    Height=300.e-9
    Angle = math.radians(AngleD)
    proj = (Height * math.cos(Angle))+ Height;  # projected thickness at some angle,
    outerMostZone =  500e-9 # m
    
    Energy = 0.6   # keV
    wavelength = 2.06e-9  # m
    f = 2.904e-3   # focal length, 2.904e-3 for 30um ZP  % 4.84e-3 for 50 um ZP
    #outerMostZone = 200e-9  # %outer most zone, [m]
    width = 0.4e-6 # [m]
    NZones =  round( wavelength * f / 4.0 / outerMostZone**2.0)  # number of zones
    
    resolution = 10.e-9  # distance between points along x-axis
    
    fig = plt.figure()

    primitives = []
    for zoneIndex in range(1,NZones):
        
        primitive = primitivePoly(zoneIndex,f,wavelength, NZones,resolution, width, shape = 'quadratic')
        x,y=map(list, zip(*primitive))
        x,y=np.asarray(x), np.asarray(y)
            
        #points, cells, point_data, cell_data, field_data = extrudePolygon2D(primitive)
        mesh  = extrudePolygon2D(primitive)
        import meshio

        fname = 'primitive_%d.mesh' % zoneIndex
        meshio.write(fname, mesh)

        mesh  = pymesh.load_mesh(fname)
         
        #fname = 'primitive_%d.mesh' % zoneIndex
        #meshio.write_points_cells(fname, points, cells, cell_data=cell_data)
        #print ("Wrote %s" % fname)
#        primitives.append({'points' : points, 
#                           'cells' : cells, 
#                           'point_data' : point_data, 
#                           'cell_data' : cell_data, 
#                           'field_data' : field_data}) 
        #primitives.append(fname) 

        # primitiveMesh = meshio.Mesh(
        #     points=points,
        #     cells=cells,
        #     point_data=point_data,
        #     cell_data=cell_data,
        #     field_data=field_data,
        # )
        #     

        primitives.append(mesh)
        #primitives.append(fname)
        #slices = pymesh.slice_mesh(primitiveMesh, [0,0,1], 4)
    
    
        plt.scatter(x*1e9, y*1e9,  s=0.3, marker='o')
    
    plt.show()
    
    # meshList =[]
    # for p in primitives:
    #      pmesh = pymesh.load_mesh(p)
    #      #meshList.append({'mesh': pmesh})
    #      meshList.append(pmesh)

        
    print ("Merging mesh")    
    mmesh = pymesh.merge_meshes(primitives)


    #tree = pymesh.CSGTree({"union":  meshList });
    #mesh = tree.mesh
    #return mesh
    return mmesh
    

 #=============================================================================
# /pymesh slicing 
# =============================================================================



if __name__ == '__main__':	
    print ("tests")
    

    #print('Test generation of single primitive')
    #testSinglePrimitive()
    
    #print('Test generation of minimal primitive set comprising a LTZP')
    mesh = testZoneOfPrimitives()
    pymesh.save_mesh('merged.stl', mesh, ascii=True);

   # pymesh.Assembler(mesh, material=None)

# =============================================================================
#     print("Slice it")
#     
# =============================================================================
   
    slices = sliceMesh(mesh,[0,0,100],5)


    


