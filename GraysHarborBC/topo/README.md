# GraysHarbor/topo

Use [GraysHarborMTL.ipynb](topo/GraysHarborMTL.html) to download 1/3 arcsecond
topography data that is referenced to MHW and adjust it so that it is 
more closely aligned with Mean Tide Level (MTL). Since we are modeling
the tide, we want to use topography referenced to MTL if possible.

Use `fetch_etopo.py` to download 1 arcminute topography of the coastal
region from the etopo1 database.

Creat`topofiles` subdirectory and move the resulting .asc files into it.
