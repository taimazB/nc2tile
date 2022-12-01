import numpy as np
from netCDF4 import Dataset
import sys
import matplotlib.pyplot as plt
from scipy import interpolate
import multiprocessing
import os
import imageio
from datetime import datetime, timedelta
import math


field = 'SSH'
model = 'CIOPS'

maxTileLat = 85.0511287798066
tileSize = 512  # px
minZoom = 2
maxZoom = 11


modelDate = sys.argv[1]
modelHr = int(sys.argv[2])
forecastHr = int(sys.argv[3])


nc = Dataset("nc/%sT%02dZ_MSC_CIOPS-SalishSea_SeaSfcHeight_Sfc_LatLon0.008x0.005_PT%03dH.nc" % (modelDate, modelHr, forecastHr), 'r')
var = nc.variables['sossheig'][0].data
missingValue = nc['sossheig'].missing_value

##  latitude, longitude, depth
lonNC = nc.variables['longitude'][:].data
latNC = nc.variables['latitude'][:].data
lonNC[lonNC >= 180] -= 360


# Mercator
R = 6378137
xNC = R * lonNC * np.pi/180.
yNC = R * np.log(np.tan(np.pi/4 + latNC*np.pi/180/2))


# Fixed min/max values for all levels and times
varMin = -3
varMax = 3
step = 0.01

var[var==missingValue] = -9999
f = interpolate.interp2d(xNC, yNC, var, kind='linear')


###################################################################
###########################  FUNCTIONS  ###########################


def xMercator(lon):
    return R * lon * np.pi/180.


def yMercator(lat):
    return R * np.log(np.tan(np.pi/4 + lat*np.pi/180/2))


def saveImg(i, j):
    try:
        if(yTile[j*tileSize: (j+1)*tileSize].min() > yNC.max() or
           yTile[j*tileSize: (j+1)*tileSize].max() < yNC.min()):
            print('Exit 1')
            return
    except:
        print('Exit 2')
        return
    dateSave = (datetime.strptime("%s %d" % (modelDate,modelHr),'%Y%m%d %H') + timedelta(hours=forecastHr)).strftime('%Y%m%d_%H')
    devNull = os.system('mkdir -p tiles/%s_%s_%s/%d/%d' %
                        (model, field, dateSave, zoom, i))
    varNew = f(xTile[i*tileSize: (i+1)*tileSize],
                  yTile[j * tileSize:(j+1) * tileSize])
    varNew[varNew < varMin] = np.nan
    #
    # To trim the interpolation tail from the right side
    iLonMax = np.argmin(np.abs(xTile-xMercator(lonNC[-1])))
    if((i+1)*tileSize > iLonMax):
        if(i*tileSize > iLonMax):
            varNew[:, :] = np.nan
        else:
            varNew[:, iLonMax % tileSize:] = np.nan
    #
    # To trim the interpolation tail from the left side
    iLonMin = np.argmin(np.abs(xTile-xMercator(lonNC[0])))
    if(i*tileSize < iLonMin):
        if((i+1)*tileSize < iLonMin):
            varNew[:, :] = np.nan
        else:
            varNew[:, :iLonMin % tileSize] = np.nan
    #
    # To trim the interpolation tail from the top side
    jLatMax = np.argmin(np.abs(yTile-yMercator(latNC[-1])))
    if((j+1)*tileSize > jLatMax):
        if(j*tileSize > jLatMax):
            varNew[:, :] = np.nan
        else:
            varNew[jLatMax % tileSize:, :] = np.nan
    #
    # To trim the interpolation tail from the bottom side
    jLatMin = np.argmin(np.abs(yTile-yMercator(latNC[0])))
    if(j*tileSize < jLatMin):
        if((j+1)*tileSize < jLatMin):
            varNew[:, :] = np.nan
        else:
            varNew[:jLatMin % tileSize, :] = np.nan
    #
    if(np.isnan(np.nanmax(varNew))):
        return
    else:
        # Coloring
        varNewRounded = np.round(varNew, int(-math.log10(step)))
        varNewInt = ((varNewRounded-varMin)/step).astype(np.uint16)
        varNewInt[varNewInt < 0] = 0
        varNewInt[varNewInt > 100 * (varMax-varMin)] = 100*(varMax-varMin)
        varColored = colors[varNewInt]
        #
        # Saving
        imageio.imwrite('tiles/%s_%s_%s/%d/%d/%d.png' % (model, field, dateSave,
                                                            zoom, i, 2**zoom-j-1), np.flipud(varColored).astype(np.uint8))


def colorRange(color1, color2, n):
    colors = []
    for r, g, b, a in zip(np.linspace(color1[0], color2[0], n), np.linspace(color1[1], color2[1], n), np.linspace(color1[2], color2[2], n), np.linspace(color1[3], color2[3], n)):
        colors.append((r, g, b, a))
    return colors


##  Continents = transparent
color0 = [(0, 0, 0, 0)]

colors1 = colorRange([0, 51, 153, 255], [0, 255, 204, 255], 120)
colors2 = colorRange([0, 255, 204, 255], [51, 204, 51, 255], 120)
colors3 = colorRange([51, 204, 51, 255], [204, 204, 0, 255], 120)
colors4 = colorRange([204, 204, 0, 255], [204, 0, 0, 255], 120)
colors5 = colorRange([204, 0, 0, 255], [255, 204, 204, 255], 120)

# All ranges
colors = np.array(color0+colors1+colors2+colors3+colors4+colors5)


def saveTile():
    global zoom
    global xTile, yTile
    #
    for zoom in np.arange(minZoom, maxZoom+1):
        print("--  Start zoom %d" % zoom)
        noOfPoints = 2**zoom*tileSize
        #
        xTile = np.linspace(xMercator(-180),
                            xMercator(180), noOfPoints)
        yTile = np.linspace(yMercator(-maxTileLat),
                            yMercator(maxTileLat), noOfPoints)
        iStart = math.floor(np.abs(xTile-xNC[0]).argmin()/tileSize)
        iEnd = math.floor(np.abs(xTile-xNC[-1]).argmin()/tileSize)+1
        jStart = math.floor(np.abs(yTile-yNC[0]).argmin()/tileSize)
        jEnd = math.floor(np.abs(yTile-yNC[-1]).argmin()/tileSize)+1
        iters = np.array(np.meshgrid(np.arange(iStart, iEnd),
                                     np.arange(jStart, jEnd))).T.reshape(-1, 2)
        #
        with multiprocessing.Pool() as p:
            p.starmap(saveImg, iters)


saveTile()
