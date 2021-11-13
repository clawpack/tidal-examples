## This directory simulate tides in Grays Harbor, WA with 3 methods. 

### First, a sine tidal signal is implemented to get 

1.The high tide time lag between the sine wave and NOAA station 1102 at Westport   

2.The high tide amplitude at NOAA station 1102

This is done in ./sine

### Second, simulates the King Tide event from December, 2015

The tidal condition is implemented using the existing NOAA 1102 tides prediction data, but the data is shifted by the time lag and divided by the amplitude observed in the sine example.

This is done in ./kingtide2015

## 3Methods.ipynb

It explains "Simple", "Riemann" and ocean forcing method.

## `topo` directory

It creates Grays Harbor topography. 


## `comparison` directory

It compares tidal predictions of the 3 methods (Simple, Riemann and ocean forcing) used for King Tide event from December, 2015


Feel free to contact yz3889@columbia.edu if you have any questions about the code.
