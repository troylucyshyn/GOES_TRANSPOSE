# GOES_TRANSPOSE

Note. Work in progress. My first time using Github.

Transpose GOES16 imagery to a different view.

## Introduction / Contents
The idea for this code originated from "https://github.com/lanceberc/GOES" who did some significant work deciphering GDAL and making it work for images files from the GOES16 and GOES17 weather satellites. The output of this code is NetCDF and/or GeoTiff that has been reprojected.

One the NetCDF has been created from jpg/png file downloaded from GOESTOOLS, there are numerous examples of how to manipulate these files to add borders, coastlines, and river as described from "https://unidata.github.io/python-gallery/examples/mapping_GOES16_TrueColor.html"

The collaboration between these two code sources allows the geostationary view to be reprojected/transposed. For example the attached files show the perspective as if the viewer is directly over Canada rather than viewing Canada from a location over the equator.



## Reprojecting GOES Images

## Adding Map Overlays

## References
CIRA / RAMMB SLIDER - Source of GeoColored PNGs
NOAA NESDIS - Source of GeoColored JPGs (and other archived data)
NASA GOES-R Series Product Users Guide Volume 4 - GOES image spatial specifications
NOAA OPC - Ocean Prediction Center weather charts
Proj.4 - Geospatial coordinate transformation library
GDAL - Geospatial Data Abstraction Library
numpy - Python scientific computing library
Pillow - A Python image manipulation library

## Done
