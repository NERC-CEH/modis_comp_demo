#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 13:17:20 2021

@author: ciaran

Use the AppEARS REST API to download some imagery with an aoi geojson

"""

import requests
import time
import os
import cgi
from joblib import Parallel, delayed
import geojson


#
#popular_layers = {"MOD13Q1.006": "_250m_16_days_NDVI", 
#               'MOD13A3.006': '_1_km_monthly_NDVI'}
#
#conus_lsat = ["CU_LC08.001", "CU_LE07.001", "CU_LT04.001", "CU_LT05.001"]

# e.g. 3 lsat7 bands
#prodLayer = [{
#        "layer": "SRB3",
#        "product": "CU_LE07.001"
#      },
#      {
#        "layer": "SRB4",
#        "product": "CU_LE07.001"
#      },
#      {
#        "layer": "SRB5",
#        "product": "CU_LE07.001"
#      } ]

def prodlist():
    
    """ 
    return a dict of APPears products
    
    Parameters
    ----------
    
    product: string
            you AppEARS password
            
    Returns
    -------
    
    a dict of APPears products
    
    Notes
    -----
    
    To look at the details do the following
    
    prodlist()
    
    products["MOD13Q1.006"]
    
    This can be used in conjunction with layer_search to ascertain which 
    product and layer you want

    """
    
    

    api = 'https://lpdaacsvc.cr.usgs.gov/appeears/api/'
    
    product_response = requests.get('{}product'.format(api)).json()

    # Print no. products available in AppEEARS
                         
    print('AρρEEARS currently supports {} products.'.format(len(product_response))) 
    

    
    prod_names = {p['ProductAndVersion'] for p in product_response} 
    
    return prod_names


# Ideally, get the product from above, then list the layers associated with it
# Then make a choice as to what you want. 
    

def layer_search(product):
    
    """ 
    return a dict of APPears layers associated with a product
    
    Parameters
    ----------
    
    product: string
            e.g. MOD13A3.006 (MODIS 1km NDVI) or CU_LC08.001 (conus lsat 8 SR)
            
    Returns
    -------
    
    a dict containing products that contain keyword

    """
    api = 'https://lpdaacsvc.cr.usgs.gov/appeears/api/'
    
    lyr_response = requests.get('{}product/{}'.format(api, product)).json()  
    
    return list(lyr_response.keys())
#    
 

def appears_download(aoi, user, password, start_date, end_date, out_dir, 
                     year_range=[2000, 2018], 
                     product_layer=[("MOD13Q1.006", "_250m_16_days_NDVI")],
                     recurring=True, proj='geographic', nt=8,
                     task='mytask', outformat='geotiff'):
    
    """ 
    Using the NASA AppEARS api to download a time series of imagery from a given
    sensor
    
    Parameters
    ----------
    
    aoi: string
         a geojson polygon/multipolygon demarcating the area of interest
          
    user: string
           your AppEARS username
           
    password: string
            you AppEARS password
    
    startdate: string
                month/day e.g. 05-01
                
    enddate: string
                month/day e.g. 05-31
    
    year_range: list of ints
                 e.g [2000, 2018]
    
    product_layer: list of tuples
             the product and corresponding layer to be downloaded
    
    recurring: bool
                whether the task recurs only in the demarcated period
    
    proj: string
             the projection of the download data 
                
    nt: int
           no of threads to use (parallel processing)
    
    task: string
            the name of the download task
    
    outformat: string
        the image format e.g. geotiff, netcdf4
        
    
    """
    


    # broadly based on the examples in the API guide
    
    inShp = aoi
    
    # get the prod layer for insertion into the json
    prod_layer = []
    
    for l in product_layer:
        prod_layer.append({
                "layer": l[1],
                "product": l[0]
              })
    
    response = requests.post('https://lpdaacsvc.cr.usgs.gov/appeears/api/login', 
                             auth=(user, password))
    token_response = response.json()
    
    # load in the json
    with open(inShp) as f:
        geoj = geojson.load(f)
   
    # create the task request json type structure
    task = {
        'task_type': 'area',
        'task_name': task,
        'params': {
             'dates': [{"endDate":  end_date,
                        "recurring": recurring, "startDate":  start_date, 
                        "yearRange": year_range}],
             'layers': prod_layer,
             'output': {
                     'format': {
                             'type': 'geotiff'}, 
                             'projection': proj},
             'geo': geoj
        }
    }
    
    # the task request
    
    # for convenience for insertion into the requests                 
    api = 'https://lpdaacsvc.cr.usgs.gov/appeears/api/'                 
    
    token = token_response['token']
    
    head = {'Authorization': 'Bearer {}'.format(token)} 
    
    response = requests.post(
        'https://lpdaacsvc.cr.usgs.gov/appeears/api/task', 
        json=task, 
        headers={'Authorization': 'Bearer {0}'.format(token)})
    task_response = response.json()
    print(task_response)
    
    # the task id which we use to track the progress or get stuff
    try:
        task_id = task_response['task_id']
    except:
        print('No taskID likely due to the product/layer combination being' 
              'incorrect')
    
    #for ref a response
    #status_response = requests.get("{}/status/{}".format(api, task_id), headers=head).json()
    
    
    # This will do for now - wait for response to change to done
    # TODO a progress spinny thing or something...
    starttime = time.time()

    while requests.get('{}task/{}'.format(api, task_id),
                       headers=head).json()['status'] != 'done':
        print(requests.get('{}task/{}'.format(api, task_id),
                           headers=head).json()['status'])
        time.sleep(20.0 - ((time.time() - starttime) % 20.0))
        
    print(requests.get('{}task/{}'.format(api, task_id), headers=head).json()['status'])
         
    bundle = requests.get("{}/bundle/{}".format(api, task_id), headers=head).json()
    
    files = {}
    for f in bundle['files']: 
        files[f['file_id']] = f['file_name'] 
        
        
    outDir = out_dir
    
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    #Why doesn't my formatting work???? - its the same!
    #    download_response = requests.get(requests.get("{}/bundle/{}/{}".format(api,
    #                                                  task_id, file), stream=True)) 
    # Anyway the f{...} seems to...
    
    # This is adapted from appears example but is far too slow. 
    # Need to get quicker tool to download the files 
    # perhaps look at sentinelsat or landsat tools see what they do. 
    def _downloader(file):
        """
        a quick func for use in parallel 
        """
        # Get a stream to the bundle file
        download_response = requests.get(f"{api}/bundle/{task_id}/{file}", 
                                         stream=True)
        # Parse the name from Content-Disposition header                                    
        filename = os.path.basename(cgi.parse_header(download_response.headers
                                                     ['Content-Disposition'])[1]['filename'])    
        # Create output file path
        filepath = os.path.join(outDir, filename) 
        # Write file to dest dir                                                                       
        with open(filepath, 'wb') as file:                                                                               
            for data in download_response.iter_content(chunk_size=8192): 
                file.write(data)
    
    # Use joblib to spread the task for speed up
    # Beyond 8 there seems to no speedup - the download libs may already 
    # muti-thread
    Parallel(n_jobs=nt, verbose=3)(delayed(_downloader)(file) for file in files)
    
    
 
    
    
    
    
    
    
    
    
    
