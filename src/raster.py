
# -*- coding: utf-8 -*-
"""
The raster module. 

Description
-----------

A series of tools for the manipulation of geospatial imagery/rasters such as
masking or raster algebraic type functions and the conversion of Sentinel 2 
data to gdal compatible formats.  

"""
from osgeo import gdal
import numpy as np
from src.gdal_merge import _merge
from tqdm import tqdm

gdal.UseExceptions()


def array2raster(array, bands, inRaster, outRas, dtype, FMT=None):
    
    """
    Save a raster from a numpy array using the geoinfo from another.
    
    Parameters
    ----------      
    array: np array
            a numpy array.
    
    bands: int
            the no of bands. 
    
    inRaster: string
               the path of a raster.
    
    outRas: string
             the path of the output raster.
    
    dtype: int 
            though you need to know what the number represents!
            a GDAL datatype (see the GDAL website) e.g gdal.GDT_Int32
    
    FMT: string 
           (optional) a GDAL raster format (see the GDAL website) eg Gtiff, HFA, KEA.
        
    
    """

    if FMT == None:
        FMT = 'Gtiff'
        
    if FMT == 'HFA':
        fmt = '.img'
    if FMT == 'KEA':
        fmt = '.kea'
    if FMT == 'Gtiff':
        fmt = '.tif'    
    
    inras = gdal.Open(inRaster, gdal.GA_ReadOnly)    
    
    x_pixels = inras.RasterXSize  # number of pixels in x
    y_pixels = inras.RasterYSize  # number of pixels in y
    geotransform = inras.GetGeoTransform()
    PIXEL_SIZE = geotransform[1]  # size of the pixel...they are square so thats ok.
    #if not would need w x h
    x_min = geotransform[0]
    y_max = geotransform[3]
    # x_min & y_max are like the "top left" corner.
    projection = inras.GetProjection()
    geotransform = inras.GetGeoTransform()   

    driver = gdal.GetDriverByName(FMT)

    dataset = driver.Create(
        outRas, 
        x_pixels,
        y_pixels,
        bands,
        dtype)

    dataset.SetGeoTransform((
        x_min,    # 0
        PIXEL_SIZE,  # 1
        0,                      # 2
        y_max,    # 3
        0,                      # 4
        -PIXEL_SIZE))    

    dataset.SetProjection(projection)
    if bands == 1:
        dataset.GetRasterBand(1).WriteArray(array)
        dataset.FlushCache()  # Write to disk.
        dataset=None
        #print('Raster written to disk')
    else:
    # Here we loop through bands
        for band in range(1,bands+1):
            Arr = array[:,:,band-1]
            dataset.GetRasterBand(band).WriteArray(Arr)
        dataset.FlushCache()  # Write to disk.
        dataset=None
        #print('Raster written to disk')
        
def raster2array(inRas, bands=[1]):
    
    """
    Read a raster and return an array, either single or multiband

    
    Parameters
    ----------
    
    inRas: string
                  input  raster 
                  
    bands: list
                  a list of bands to return in the array
    
    """
    rds = gdal.Open(inRas)
   
   
    if len(bands) ==1:
        # then we needn't bother with all the crap below
        inArray = rds.GetRasterBand(bands[0]).ReadAsArray()
        
    else:
        #   The nump and gdal dtype (ints)
        #   {"uint8": 1,"int8": 1,"uint16": 2,"int16": 3,"uint32": 4,"int32": 5,
        #    "float32": 6, "float64": 7, "complex64": 10, "complex128": 11}
        
        # a numpy gdal conversion dict - this seems a bit long-winded
        dtypes = {"1": np.uint8, "2": np.uint16,
              "3": np.int16, "4": np.uint32,"5": np.int32,
              "6": np.float32,"7": np.float64,"10": np.complex64,
              "11": np.complex128}
        rdsDtype = rds.GetRasterBand(1).DataType
        inDt = dtypes[str(rdsDtype)]
        
        inArray = np.zeros((rds.RasterYSize, rds.RasterXSize, len(bands)), dtype=inDt) 
        for idx, band in enumerate(bands):  
            rA = rds.GetRasterBand(band).ReadAsArray()
            inArray[:, :, idx]=rA
   
   
    return inArray




def stack_ras(rasterList, outFile):
    """ 
    Stack some rasters 
    
    Parameters
    ----------- 
        
    rasterList: string
             the input image 
        
    outFile: string
              the output file path including file extension
        

    """
    _merge(names = rasterList, out_file = outFile)
    


def stat_comp(inRas, outMap, bandList = None,  stat = 'percentile', q = 95, 
                  blocksize=256,
                  FMT=None,  dtype = gdal.GDT_Float32):
    
    """
    Calculate depth wise stat on a multi band raster with selected or all bands
            
    Parameters 
    ---------- 
    
    inRas: string
               input Raster
    
    outMap: string
             the output raster calculated

    stat: string
           the statisitc to be calculated make sure there 
           are no nans as nan percentile is far too slow        

    blocksize: int
                the chunck processed 

    q: int
        the ith percentile if percentile is the stat used         
    
    FMT: string
          gdal compatible (optional) defaults is tif

    dtype: string
            gdal datatype (default gdal.GDT_Int32)
    """
    
    

    inDataset = gdal.Open(inRas)
    
    
    ootbands = len(bandList)
    
    bands = inDataset.RasterCount
    
    outDataset = _copy_dataset_config(inDataset, outMap = outMap,
                                     bands = 1)
        
    band = inDataset.GetRasterBand(1)
    cols = inDataset.RasterXSize
    rows = inDataset.RasterYSize
    
    # So with most datasets blocksize is a row scanline
    if blocksize == None:
        blocksize = band.GetBlockSize()
        blocksizeX = blocksize[0]
        blocksizeY = blocksize[1]
    else:
        blocksizeX = blocksize
        blocksizeY = blocksize

    def statChoose(X, stat, q):
        if stat == 'mean':
            stats = np.nanmean(X, axis=2)
        if stat == 'std':
            stats = np.nanstd(X, axis=2)
        if stat == 'percentile':
            # slow as feck
            stats = np.percentile(X, q, axis=2)    
        if stat == 'median':
            stats = np.nanmedian(X, axis=2)
            
        
        return stats
    
        
    if blocksizeY==1:
        rows = np.arange(outDataset.RasterYSize, dtype=np.int)                
        for row in tqdm(rows):
            i = int(row)
            j = 0
            X = np.zeros(shape = (blocksizeY , blocksizeX, ootbands))

            for ind, im in enumerate(bandList):
                array = inDataset.GetRasterBand(im).ReadAsArray(j,i,
                                            blocksizeX, j)
                array.shape = (1, blocksizeX)
                X[:,:,ind] = array
                    
            stArray = statChoose(X, stat, q)

            outDataset.GetRasterBand(band).WriteArray(stArray,j,i)


    # else it is a block            
    else:
        for i in tqdm(range(0, rows, blocksizeY)):
            if i + blocksizeY < rows:
                numRows = blocksizeY
            else:
                numRows = rows -i
        
            for j in range(0, cols, blocksizeX):
                if j + blocksizeX < cols:
                    numCols = blocksizeX
                else:
                    numCols = cols - j
                X = np.zeros(shape = (numRows, numCols, ootbands))

                for ind, im in enumerate(bandList):
                    array = inDataset.GetRasterBand(im).ReadAsArray(j,i,
                                            numCols, numRows)
                    
                    X[:,:,ind] = array
                    
                stArray = statChoose(X, stat, q)

                outDataset.GetRasterBand(1).WriteArray(stArray,j,i)

                    #print(i,j)
    outDataset.FlushCache()
    outDataset = None    



def _copy_dataset_config(inDataset, FMT = 'Gtiff', outMap = 'copy',
                         dtype = gdal.GDT_Int32, bands = 1):
    """Copies a dataset without the associated rasters.

    """

    
    x_pixels = inDataset.RasterXSize  # number of pixels in x
    y_pixels = inDataset.RasterYSize  # number of pixels in y
    geotransform = inDataset.GetGeoTransform()
    PIXEL_SIZE = geotransform[1]  # size of the pixel...they are square so thats ok.
    #if not would need w x h
    x_min = geotransform[0]
    y_max = geotransform[3]
    # x_min & y_max are like the "top left" corner.
    projection = inDataset.GetProjection()
    geotransform = inDataset.GetGeoTransform()   
    #dtype=gdal.GDT_Int32
    driver = gdal.GetDriverByName(FMT)
    
    # Set params for output raster
    outDataset = driver.Create(
        outMap, 
        x_pixels,
        y_pixels,
        bands,
        dtype)

    outDataset.SetGeoTransform((
        x_min,    # 0
        PIXEL_SIZE,  # 1
        0,                      # 2
        y_max,    # 3
        0,                      # 4
        -PIXEL_SIZE))
        
    outDataset.SetProjection(projection)
    
    return outDataset

def _quickwarp(inRas, outRas, proj='EPSG:27700'):
    
    """gdalwarp a dataset

    """
    ootRas = gdal.Warp(outRas, inRas, dstSRS=proj, format='Gtiff')
    ootRas.FlushCache()
    ootRas=None

