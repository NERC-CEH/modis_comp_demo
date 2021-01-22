A brief demo of the process to derive May composites of MODIS NDVI and attribute landcover polygons with FOSS libs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is an ipynb demo for both CEH land cover and 16 day data and one for Corine and monthly data. Both of course use gpkg subsets as the entire datasets would result in a fair wait (e.g. 6737558 features on CEH_lcm). 

It should be noted that the input rasters have already been averaged by NASA/USGS prior to the averaging below.

**Please register with APPears so you can use your username and password in the workflow**

https://lpdaacsvc.cr.usgs.gov/appeears/

**The first part utilises a function I have written to access the APPears API to download the desired data within an AOI geojson.**

**The second stacks and composites rasters then attributes the shapefiles using functions from my own lib geospatial_learn.** 

https://github.com/Ciaran1981/geospatial-learn

The dependencies for the above lib are numerous and it will take a while to install hence, for brevity, I have just provided cut down modules locally. Provided you are using the notebook herein you can import the functions from the local files. The files are appears_down.py, raster.py and shape.py and reside in src. I have left the original lib imports commented out if you decide to install/use it.

Installation and Use
~~~~~~~~~~~~~~~~~~~~

Installing the required libs (there are not many) uses the conda system so ensure you have this first. Clone this repo and cd into the directory then...

.. code-block:: bash

conda env create -f ndvi_demo.yml

conda activate ndvi_demo

jupyter notebook

Then open the ipynb and cycle through the cells.

