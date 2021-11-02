## `topo` directory

Use [GraysHarborMTL.ipynb](topo/GraysHarborMTL.html) to download 1/3 arcsecond
topography data that is referenced to MHW and adjust it so that it is 
more closely aligned with Mean Tide Level (MTL). Since we are modeling
the tide, we want to use topography referenced to MTL if possible.

Use `fetch_etopo.py` to download 1 arcminute topography of the coastal
region from the etopo1 database.

## `sine` directory

In this directory boundary condition is applied using
pure sine wave oscillation with amplitude 1. Examining the resulting gauge
time series at the Westport gauge allows us to determine how much this wave
is amplified and the time lag between the boundary signal and the response
at Westport.  

## `kingtide2015` directory

Simulates the King Tide event from December, 2015.  
The notebook downloads observations at the Westport and Aberdeen tide gauges and adjusts 
the Westport data to obtain the tidal signal to be applied on the left boundary.
After running GeoClaw with this input, the same notebook shows the
resulting simulated tides compared to the observations.

