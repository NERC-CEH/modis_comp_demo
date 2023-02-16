# -*- coding: utf-8 -*-
"""
The shape module. 

Description
-----------

A couple of functions from the geospatial_learn.shape module (
https://github.com/Ciaran1981/geospatial-learn) for the purposes
of the demo of processing for part of SOC-D


"""
from scipy.stats import skew, kurtosis
from osgeo import gdal, ogr
from tqdm import tqdm
import numpy as np
from scipy.stats.mstats import mode
from pandas import DataFrame
gdal.UseExceptions()
ogr.UseExceptions()




        
def _bbox_to_pixel_offsets(rgt, geom):
    
    """ 
    Internal function to get pixel geo-locations of bbox of a polygon
    
    Parameters
    ----------
    
    rgt : array
          List of points defining polygon (?)
          
    geom : shapely.geometry
           Structure defining geometry
    
    Returns
    -------

    x offset: int
           
    y offset: int
           
    xcount: int
             rows of bounding box
             
    ycount: int
             columns of bounding box
    """
    
    xOrigin = rgt[0]
    yOrigin = rgt[3]
    pixelWidth = rgt[1]
    pixelHeight = rgt[5]
    ring = geom.GetGeometryRef(0)
    numpoints = ring.GetPointCount()
    pointsX = []; pointsY = []
    
    if (geom.GetGeometryName() == 'MULTIPOLYGON'):
        count = 0
        pointsX = []; pointsY = []
        for polygon in geom:
            geomInner = geom.GetGeometryRef(count)
            ring = geomInner.GetGeometryRef(0)
            numpoints = ring.GetPointCount()
            for p in range(numpoints):
                    lon, lat, z = ring.GetPoint(p)
                    pointsX.append(lon)
                    pointsY.append(lat)
            count += 1
    elif (geom.GetGeometryName() == 'POLYGON'):
        ring = geom.GetGeometryRef(0)
        numpoints = ring.GetPointCount()
        pointsX = []; pointsY = []
        for p in range(numpoints):
                lon, lat, z = ring.GetPoint(p)
                pointsX.append(lon)
                pointsY.append(lat)
            
    xmin = min(pointsX)
    xmax = max(pointsX)
    ymin = min(pointsY)
    ymax = max(pointsY)

    # Specify offset and rows and columns to read
    xoff = int((xmin - xOrigin)/pixelWidth)
    yoff = int((yOrigin - ymax)/pixelWidth)
    xcount = int((xmax - xmin)/pixelWidth)+1
    ycount = int((ymax - ymin)/pixelWidth)+1

    return (xoff, yoff, xcount, ycount)        

def _fieldexist(vlyr, field):
    """
    check a field exists
    """
    
    lyrdef = vlyr.GetLayerDefn()

    fieldz = []
    for i in range(lyrdef.GetFieldCount()):
        fieldz.append(lyrdef.GetFieldDefn(i).GetName())
    return field in fieldz

def zonal_stats(inShp, inRas, band, bandname, layer=None, stat = 'mean',
                write_stat=True, nodata_value=0, all_touched=True):
    
    """ 
    Calculate zonal stats for an OGR polygon file
    
    Parameters
    ----------
    
    inShp: string
                  input shapefile
        
    inRas: string
                  input raster

    band: int
           an integer val eg - 2

    bandname: string
               eg - blue
    layer: string
           if using a db type format with multi layers, specify the name of the
           layer in question
           
    stat: string
           string of a stat to calculate, if omitted it will be 'mean'
           others: 'mode', 'min','mean','max', 'std',' sum', 'count','var',
           skew', 'kurt (osis)'
                     
    write_stat: bool (optional)
                If True, stat will be written to OGR file, if false, dataframe
                only returned (bool)
        
    nodata_value: numerical
                   If used the no data val of the raster
    
    all_touched: bool
                    whether to use all touched when raterising the polygon
                    if the poly is smaller/comaparable to the pixel size, 
                    True is perhaps the best option
        
    """    
    # gdal/ogr-based zonal stats
    
    if all_touched == True:
        touch = "ALL_TOUCHED=TRUE"
    else:
        touch = "ALL_TOUCHED=FALSE"
        
    rds = gdal.Open(inRas, gdal.GA_ReadOnly)

    rb = rds.GetRasterBand(band)
    rgt = rds.GetGeoTransform()

    if nodata_value:
        nodata_value = float(nodata_value)
        rb.SetNoDataValue(nodata_value)

    vds = ogr.Open(inShp, 1) 
    
    # if we are using a db of some sort gpkg etc where we have to choose
    if layer !=None:
        vlyr = vds.GetLayerByName(layer)
    else:
        vlyr = vds.GetLayer()
    if write_stat != None:
        # if the field exists leave it as ogr is a pain with dropping it
        # plus can break the file
        if _fieldexist(vlyr, bandname) == False:
            vlyr.CreateField(ogr.FieldDefn(bandname, ogr.OFTReal))

    mem_drv = ogr.GetDriverByName('Memory')
    driver = gdal.GetDriverByName('MEM')

    # Loop through vectors
    stats = []
    feat = vlyr.GetNextFeature()
    features = np.arange(vlyr.GetFeatureCount())
    rejects = list()
    for label in tqdm(features):

        if feat is None:
            continue
#        debug
#        wkt=geom.ExportToWkt()
#        poly1 = loads(wkt)
        geom = feat.geometry()

        src_offset = _bbox_to_pixel_offsets(rgt, geom)
        src_array = rb.ReadAsArray(src_offset[0], src_offset[1], src_offset[2],
                               src_offset[3])
        if src_array is None:
            src_array = rb.ReadAsArray(src_offset[0]-1, src_offset[1], src_offset[2],
                               src_offset[3])
            if src_array is None:
                rejects.append(feat.GetFID())
                continue

        # calculate new geotransform of the feature subset
        new_gt = (
        (rgt[0] + (src_offset[0] * rgt[1])),
        rgt[1],
        0.0,
        (rgt[3] + (src_offset[1] * rgt[5])),
        0.0,
        rgt[5])

            
        # Create a temporary vector layer in memory
        mem_ds = mem_drv.CreateDataSource('out')
        mem_layer = mem_ds.CreateLayer('poly', None, ogr.wkbPolygon)
        mem_layer.CreateFeature(feat.Clone())

        # Rasterize it

        rvds = driver.Create('', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
     
        rvds.SetGeoTransform(new_gt)
        rvds.SetProjection(rds.GetProjectionRef())
        rvds.SetGeoTransform(new_gt)
        gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1], options=[touch])
        rv_array = rvds.ReadAsArray()
        
        # Mask the source data array with our current feature using np mask     

        masked = np.ma.MaskedArray(
            src_array,
            mask=np.logical_or(
                src_array == nodata_value,
                np.logical_not(rv_array)
            )
        )
        
        if stat == 'mode':
            feature_stats = mode(masked)[0]
        if stat == 'min':
            feature_stats = float(masked.min())
        if stat == 'mean':
            feature_stats = float(masked.mean())
        if stat == 'max':
            feature_stats = float(masked.max())
        if stat == 'median':
            feature_stats = float(np.median(masked[masked.nonzero()]))
        if stat == 'std':
            feature_stats = float(masked.std())
        if stat == 'sum':
            feature_stats = float(masked.sum())
#        elif stat is 'count':
#            feature_stats = int(masked.count())
        if stat == 'var':
            feature_stats = float(masked.var())
        if stat == 'skew':
            feature_stats = float(skew(masked[masked.nonzero()]))
        if stat == 'kurt':
            feature_stats = float(kurtosis(masked[masked.nonzero()]))
        
        stats.append(feature_stats)
        if write_stat != None:
            feat.SetField(bandname, feature_stats)
            vlyr.SetFeature(feat)
        feat = vlyr.GetNextFeature()
    if write_stat != None:
        vlyr.SyncToDisk()



    vds = None
    rds = None
    frame = DataFrame(stats)
    
    if write_stat != None:
        return frame
    
def zonal_stats_all(inShp, inRas, bandnames, 
                    statList = ['mean', 'min', 'max', 'median', 'std',
                                'var', 'skew', 'kurt']):
    """ 
    Calculate zonal stats for an OGR polygon file
    
    Parameters
    ----------
    
    inShp: string
                  input shapefile
        
    inRas: string
                  input raster

    band: int
           an integer val eg - 2

    bandnames: list
               eg - ['b','g','r','nir']
        
    nodata_value: numerical
                   If used the no data val of the raster
        
    """    

# zonal stats
    for bnd,name in enumerate(bandnames):
    
        [zonal_stats(inShp, inRas, bnd+1, name+st, stat=st, write_stat = True) for st in statList]

