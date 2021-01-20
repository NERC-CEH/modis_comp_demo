A brief demo of the process to derive May composites of MODIS NDVI and attribute
~~~~~~~~~~~~~~

It should be noted that the input rasters have already been averaged by NASA/USGS prior to the averaging below.

**Please register with APPears so you can use your username and password in the workflow**

https://lpdaacsvc.cr.usgs.gov/appeears/

**The first part utilises a function I have written to access the APPears API to download the desired data within an AOI geojson.**

**The second stacks and composites rasters then attributes the shapefiles using functions from my own lib geospatial_learn.** 

https://github.com/Ciaran1981/geospatial-learn

The dependencies for this are numerous and it will take a long time to install hence, for brevity, I have just provided the necessary functions locally, so provided you work in this folder (which will be the case) you can import the function. The files are appears_down.py, raster.py and shape.py and reside in src. 

Installation and Use
~~~~~~~~~~~~~~~~~~~~

To install the required libs (there are not many) uses the conda system so ensure you have this first. The libs are all just standard python science and gdal so shouldn't cause problems in your base env, but to be safe just create an env. name something like below:

.. code-block:: bash

conda env create -f ndvi_demo.yml

conda activate ndvi_demo

jupyter notebook

Then open the ipynb and cycle through the cells.

