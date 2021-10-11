# -*- coding: utf-8 -*-

#ReTranspose 2 but with coordinates for all of Canada.

def gdal_error_handler(err_class, err_num, err_msg):
    errtype = {
            gdal.CE_None:'None',
            gdal.CE_Debug:'Debug',
            gdal.CE_Warning:'Warning',
            gdal.CE_Failure:'Failure',
            gdal.CE_Fatal:'Fatal'
    }
    err_msg = err_msg.replace('\n',' ')
    err_class = errtype.get(err_class, 'None')
    print('Error Number: %s' % (err_num))
    print('Error Type: %s' % (err_class))
    print('Error Message: %s' % (err_msg))



def transpose(sourcefn,outfilename):
    import sys
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import metpy  # noqa: F401
    import numpy as np
    from PIL import Image
    from osgeo import gdal #, ogr, osr
    import matplotlib
    import matplotlib.pyplot as plt
    import xarray
    import subprocess
    import gc
    #https://gdal.org/tutorials/raster_api_tut.html
    
    version_num = int(gdal.VersionInfo('VERSION_NUM'))
    if version_num<1100000:
        sys.exit("ERROR: Python bindings of GDAL 1.10 or later required")
    
    # Enable GDAL/OGR exceptions
    #gdal.UseExceptions()
    gdal.DontUseExceptions() # To disable GDAL/OGR exceptions at any point
    
    #Set variables for 
    #proj_lcc = "+proj=lcc +lon_0=-113 +lat_0=55 +lat_1=45 +lat_2=62" # This is for West. Canada
    proj_lcc = "+proj=lcc +lon_0=-90 +lat_0=60 +lat_1=40 +lat_2=65" # This is for all of Canada
   
    (w, h) = (1280, 960) # 16:9 aspect ratio HD-sized (not FullHD which is 1920x1280)
    #interestArea = [-130, 45, -85, 62] # This is is for West. Canada
    interestArea = [-130, 30, 10, 70] # This is for all of Canada
    proj_anti_mercator = "+proj=merc +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +over +lon_0=-180"
        
    goes16 = { # GOES-16, aka GOES-R, GOES-EAST
        "p_height": 35786023.0,           # perspective height from the ellipsoid
        "height": 42164160.0,             # from center of the earth
        "longitude": -75.0,
        "sweep_axis": 'x',
        "semi_major": 6378137.0,          # GRS-80 ellipsoid
        "semi_minor": 6356752.31414,      # GRS-80 ellipsoid
        "flattening": 298.257222096,
        "eccentricity": 0.0818191910435,
        # The other resolution (.5k, 4k, etc) can be added here
        "1k": {
            "resolution": 0.000028,       # radians per pixel
            "FD": {
                "x_offset": -0.151858,    # radians from nadir
                "y_offset":  0.151858,
                "shape": (10848, 10848),  # pixels in image
            },
            "CONUS": {
                "x_offset": -0.101346,
                "y_offset":  0.128226,
                "shape": (5000, 3000),
            }
        },
        "2k": {
            "resolution": 0.000056,
            "FD": {
                "x_offset": -0.151844,
                "y_offset":  0.151844,
                "shape": (5424, 5424),
            },
            "CONUS": {
                "x_offset": -0.101332,
                "y_offset":  0.128212,
                "shape": (2500, 1500),
            }
        }
    }
    
    goes17 = { # GOES-17, aka GOES-S, GOES-WEST
        "p_height": 35786023.0,
        "height": 42164160.0,
        "longitude": -137.0,
        "sweep_axis": 'x',
        "semi_major": 6378137.0,
        "semi_minor": 6356752.31414,
        "flattening": 298.257222096,
        "eccentricity": 0.0818191910435,
        "1k": {
            "resolution": 0.000028,
            "FD": {
                "x_offset": -0.151858,
                "y_offset":  0.151858,
                "shape": (10848, 10848),
            },
            "CONUS": { # aka PACUS because it's mostly Pacific Ocean
                "x_offset": -0.069986,
                "y_offset":  0.128226,
                "shape": (5000, 3000),
            }
        },
        "2k": {
            "resolution": 0.000056,
            "FD": {
                "x_offset": -0.151844,
                "y_offset":  0.151844,
                "shape": (5424, 5424),
            },
            "CONUS": {
                "x_offset": -0.069972,
                "y_offset":  0.128212,
                "shape": (2500, 1500),
            }
        },
    }
    
    if "GOES16" in sourcefn: #this changes upon satellite
        s = goes16
    elif "GOES17" in sourcefn: #this changes upon satellite:
        s = goes17
    else:
        sys.exit("Incompatiable filename.")
    res = s["2k"]
    sector = res["FD"]
    print("The source file is: " + sourcefn)
    destfn = sourcefn[0:len(sourcefn)-4] + "TRNSP.nc"
    print("The destination file is: " + destfn)
     
    #print ("Sector X offset is: ", sector['x_offset'])
    upper_left_x = sector['x_offset'] * s['p_height']
    upper_left_y = sector['y_offset'] * s['p_height']
    resolution_m = res['resolution'] * s['p_height']
    
    geotransform = [upper_left_x, resolution_m, 0, upper_left_y, 0, -resolution_m]
    proj = "+proj=geos +lon_0=%f +h=%f +a=%f +b=%f +f=%f +units=m +no_defs -ellps=GRS80 +sweep=%s +over" % (
        s['longitude'],
        s['p_height'],
        s['semi_major'],
        s['semi_minor'],
        1/s['flattening'], # Inverese flattening
        s['sweep_axis'])
    WKT = 'PROJCS["unnamed",GEOGCS["unnamed ellipse",DATUM["unknown",SPHEROID["unnamed",%f,%f]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]\
    ,PROJECTION["Geostationary_Satellite"],PARAMETER["central_meridian",%f],PARAMETER["satellite_height",%f],PARAMETER["false_easting",0],PARAMETER["fal\
    se_northing",0],UNIT["Meter",1],EXTENSION["PROJ4","%s"]]' % (s['semi_major'], s['flattening'], s['longitude'], s['p_height'], proj)
    
    #print("This what WKT looks like: ", WKT )
    #Setup for GTiff Format
    #prewarpOptions = gdal.WarpOptions( #https://gdal.org/drivers/raster/index.html#raster-drivers
        #format="netCDF",  #https://svn.osgeo.org/gdal/tags/gdal_1_2_5/frmts/formats_list.html
        #format="MEM",  #https://svn.osgeo.org/gdal/tags/gdal_1_2_5/frmts/formats_list.html
        #format="PNG",  #https://svn.osgeo.org/gdal/tags/gdal_1_2_5/frmts/formats_list.html
    #    format="GTiff",  #https://svn.osgeo.org/gdal/tags/gdal_1_2_5/frmts/formats_list.html
    #    width=w,
    #    height=h,
    #    outputBoundsSRS="EPSG:4326", # WGS84 - Allows use of lat/lon outputBounds
    #    outputBounds=interestArea,
    #    dstSRS=proj_lcc,
        #dstSRS=proj_anti_mercator,
    #    warpOptions=["SOURCE_EXTRA=500","WRITE_GDAL_TAGS=YES"], # Magic from The Internet
    #    multithread = True,
        #creationOption = listofcreationOptions
        # Not sure if this buys anything on my setup
    #)
    
    # Setup for netCDF options
    warpOptions = gdal.WarpOptions( #https://gdal.org/drivers/raster/index.html#raster-drivers
        format="netCDF",  #https://svn.osgeo.org/gdal/tags/gdal_1_2_5/frmts/formats_list.html
        #format="MEM",  #https://svn.osgeo.org/gdal/tags/gdal_1_2_5/frmts/formats_list.html
        width=w,
        height=h,
        outputBoundsSRS="EPSG:4326", # WGS84 - Allows use of lat/lon outputBounds
        outputBounds=interestArea,
        dstSRS=proj_lcc,
        warpOptions=["SOURCE_EXTRA=500","WRITE_GDAL_TAGS=YES"], # Magic from The Internet
        multithread = True,
        #creationOption = listofcreationOptions
        # Not sure if this buys anything on my setup
    )
    
    src = gdal.Open(sourcefn, gdal.GA_ReadOnly)
    #src.RasterCount #Identifies number of bands
    #print(str(src.RasterCount) + " Band(s) in: " + sourcefn)
    
    src.SetProjection(WKT)            # Projection and Georeferencing (from satellite longitude)
    src.SetGeoTransform(geotransform) # Scan angle and pixel width ...Need to look at this
    
    #Creates TIF from the prewarpOptions
    #destfnTIFF = destfn[0:len(destfn)-4] +".tif"
    #3predst = gdal.Warp(destfnTIFF, src, options=prewarpOptions)
    #destfnTIFF = None
    #if not predst:
    #    print("Warp failed %s" % (fn))
    #    sys.exit("GTiff Warp Option Failed.")
    
    #Creats netCDF file from warpOption
    #print("The destination filename is: " + destfn)
    dst = gdal.Warp(destfn, src, options=warpOptions)
    #Currently you would need for example to create a MEM dataset, and set the NETCDF_VARNAME metadata item on its bands,
    #and finally use CreateCopy() to generate the netCDF file.
    if not dst:
        print("Warp failed %s")
        sys.exit("netCDF Warp Option Failed.")
    
    # dsta = dst.ReadAsArray()          # Array shape is [band, row, col] - At this point the file is a numpndarray
    # print("dsta is: ", dsta)
    # if src.RasterCount == 1:
    #     arr = dsta # okay for a single data change 13, 14, etc.
    #     img = Image.fromarray(arr, 'P') # Convert to PIL image. The L should work for a single colour
    #     print("Single band")
    # elif src.RasterCount == 3:
    #     arr = dsta.transpose (1,2,0)   # Virtually change the shape to [row, col, band] - Works for colour but that's all
    #     img = Image.fromarray(arr, 'RGB') # Convert to PIL image. This works good for color but not single color.
    #     print("Three bands")
    # else:
    #     sys.exit("Error with number of layers.")
     
    # moddestfn = destfn[0:len(destfn)-4] +"-gdalPIL.png"
    # img.save(moddestfn, "PNG")
    # img.show()
    # dsta = None
    # arr = None
    # Clean up - important if running in a loop
    #src = None
    
    dst = None
    predst = None
    #src = None
    moddestfn = None
    #destfn=None
    
    print("Done RETRANSPOSE")
    
    ##Start the mapping process
    
    #https://www.naturalearthdata.com/downloads/
    
    D = xarray.open_dataset(destfn)
    
    # Load the three channels into appropriate R, G, and B variables
    if src.RasterCount == 1:
        single = D['Band1'].data
    elif src.RasterCount == 3:
        R = D['Band1'].data
        G = D['Band2'].data
        B = D['Band3'].data
        RGB = np.dstack([R, G, B]) # The RGB array for the true color image
    
    
    #sys.exit()
    ######################################################################
    # Plot with `Cartopy` Geostationary Projection
    # ----------------------------------------------
    #
    # The image above is not georeferenced. You can see the land and oceans, but we
    # do have enough information to draw state and country boundaries. Use the
    # `metpy.io` package to obtain the projection information from the file.  Then
    # use `Cartopy` to plot the image on a map. The GOES data and image is on a
    # [geostationary projection
    # ](https://proj4.org/operations/projections/geos.html?highlight=geostationary).
    
    # We'll use the `CMI_C02` variable as a 'hook' to get the CF metadata.
    dat = D.metpy.parse_cf('Band1')
    D = None
    geos = dat.metpy.cartopy_crs
    
    # We also need the x (north/south) and y (east/west) axis sweep of the ABI data
    x = dat.x
    y = dat.y
    
    fig = plt.figure(figsize=(18, 15),dpi=400)
    # Generate an Cartopy projection
    #lc = ccrs.LambertConformal(central_longitude=-107.5, standard_parallels=(49,60.0)) # For. W.Canada
    lc = ccrs.LambertConformal(central_longitude=-90, standard_parallels=(49,60.0)) # For all of Canada 
    
    ax = fig.add_subplot(1, 1, 1, projection=lc)
    #ax.set_extent([-130, -94, 48, 60], crs=ccrs.PlateCarree()) # sets the mapping range but not the image range. for W. Canada
    ax.set_extent([-130, -50, 40, 70], crs=ccrs.PlateCarree()) # sets the mapping range but not the image range. For all of Canada.
    

    if src.RasterCount == 1:
        img = ax.imshow(single, origin='lower', extent=(x.min(), x.max(), y.min(),
                        y.max()), transform=geos, cmap="jet") #note added jet in change color
    elif src.RasterCount == 3:
        img = ax.imshow(RGB, origin='lower',extent=(x.min(), x.max(),
                        y.min(), y.max()), transform=geos)
    
    ax.coastlines(resolution='50m', color='black', linewidth=0.5)
    hum_lon, hum_lat = -105.1229, 52.2020
    ax.plot(hum_lon, hum_lat, marker="o",markersize=5, linewidth =1.0,color='black',transform=ccrs.Geodetic())
    ax.text(hum_lon+.2, hum_lat+.2, 'Humboldt', horizontalalignment='left', transform=ccrs.Geodetic())
    
    swft_lon, swft_lat = -107.7972, 50.2851
    ax.plot(swft_lon, swft_lat, marker="o",markersize=5, linewidth =1.0,color='black',transform=ccrs.Geodetic())
    ax.text(swft_lon+.2, swft_lat+.2, 'SwiftCurrent', horizontalalignment='left', transform=ccrs.Geodetic())
    
    resol = '50m'
    country_bodr = cfeature.NaturalEarthFeature(category='cultural',
        name='admin_0_boundary_lines_land', scale=resol, facecolor='none', edgecolor='k')
    provinc_bodr = cfeature.NaturalEarthFeature(category='cultural', 
        name='admin_1_states_provinces_lines', scale=resol, facecolor='none', edgecolor='k')
    
    lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='b', facecolor='none')
    rivers = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', scale=resol, edgecolor='b', facecolor='none')
    
    
    # Add all features to the map
    ax.add_feature(lakes, linewidth=0.5)
    ax.add_feature(rivers, linewidth=0.5)
    ax.add_feature(country_bodr, linestyle='--', linewidth=0.8, edgecolor="k")  #USA/Canada
    ax.add_feature(provinc_bodr, linestyle='--', linewidth=0.6, edgecolor="k")
    
    #plt.colorbar(img, label='Brightness Temperatures (Â°C)', extend='both', orientation='horizontal', pad=0.05, fraction=0.05)
    mapfilename = destfn[0:len(destfn)-3]+"+mapped.jpg"
    plt.title(mapfilename, loc='left', fontweight='bold', fontsize=12)
    plt.savefig(outfilename, bbox_inches="tight")
    #plt.show()
    print("The jpg desgination is: ", destfn[0:len(destfn)-4])
    #subprocess.call(['rm', destfn]) # Can't do in windows
    #subprocess.call(['rm', sourcefn + ".aux.xml"])
    destfn = None
    src = None
    gc.collect(generation=2)