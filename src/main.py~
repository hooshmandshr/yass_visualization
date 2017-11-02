
from bokeh.layouts import row, widgetbox
from bokeh.models import Select, HoverTool
from bokeh import palettes
from bokeh.plotting import curdoc, ColumnDataSource
from bokeh.plotting import figure, show, output_file, gridplot
from bokeh.sampledata.autompg import autompg

import numpy as np
import scipy.io

## matlab file inputs
dic = scipy.io.loadmat('test.mat')
TEMP_GT = dic['templates_gt']
TEMP_DT = dic['templates']
GEOM = dic['geometry']
K_MATCH = dic['k_match'][0]
TP = dic['tp'][0]

nchan = GEOM.shape[0]
nUnit = TEMP_GT.shape[2]

## interface global parameters
COLORS = palettes.Spectral5
COLORS = palettes.Inferno256
COLORS = palettes.Plasma256

OPT = [str(i) for i in range(0, nUnit)]
Y_RANGE = (np.min(TEMP_GT), np.max(TEMP_GT))

#OPT = ['unit %s' %i for i in range(0, nUnit)]
def covEllipse(X, weight):
    mu = np.average(GEOM, axis=0, weights=weight)
    covar = np.cov(GEOM.T, aweights=weight)
    eig = np.linalg.eig(covar)
    width = np.sqrt(eig[0][0])
    height = np.sqrt(eig[0][1])
    v = eig[1][:, np.argmax(eig[0])]
    angle = np.arctan(v[1]/(v[0]+1e-6))
    return mu, width, height, angle

## static overlayn figure
def getSource(temps, tp):
    _nUnit = temps.shape[2]
    maxEnergy = []
    for i in range(0, _nUnit):
        template = temps[:, :, i]
        energy = np.linalg.norm(template, axis=0)
        maxEnergy.append(np.max(energy))
    maxEnergy = max(maxEnergy) 

    MU = [[],[]]
    WIDT = []
    HGHT = []
    ANGL = []
    COLS = []
    COLSTP = []
    for i in range(0, _nUnit):
        template = temps[:, :, i]
        energy = np.linalg.norm(template, axis=0)
        energyS = np.power(energy, 2)
        mainChan = np.argmax(energy)

        mu, width, height, angle = covEllipse(GEOM, energyS)

        MU[0].append(mu[0])
        MU[1].append(mu[1])
        WIDT.append(width)
        HGHT.append(height)
        ANGL.append(angle)
        COLS.append(COLORS[int(np.max(energy)/maxEnergy*255)])
        COLSTP.append(COLORS[int(tp[i]*255)])

    covSource = ColumnDataSource(data=dict(\
        x=MU[0],\
        y=MU[1],\
        width=WIDT,\
        height=HGHT,\
        angle=ANGL,\
        color=COLS,\
        colortp=COLSTP,\
        desc=['unit %s' %i for i in range(0, nUnit)],\
        ))

    return covSource
##
gtSource = getSource(TEMP_GT, TP)

TPD = np.zeros(TEMP_DT.shape[2])
for i, j in zip(range(0, len(K_MATCH)), K_MATCH):
    TPD[j-1] = TP[i]
dtSource = getSource(TEMP_DT, TPD)

def spatialTrace(I, TEMP, dSource, plth, pltw):
    temp = TEMP[:, :, I]
    energy = np.linalg.norm(temp, axis=0)
    energyS = np.power(energy, 2)
    mainChan = np.argmax(energy)
    p = figure(plot_height=plth, plot_width=pltw, title='Unit Spatial Trace', tools='pan,box_zoom,reset,hover')
    # draw the electrodes
    p.scatter(GEOM[:, 0], GEOM[:, 1], radius=1, color=COLORS[0], line_color=None)
    # draw the spatial trace
    p.scatter(GEOM[:, 0], GEOM[:, 1],\
            radius=energy,\
            fill_alpha=0.6,\
            color=COLORS[-100],\
            line_color=None)
    # draw the covariance of spatial trace
    for i in range(1, 4):
        p.oval(x=dSource.data['x'][I],\
                y=dSource.data['y'][I],\
                width=i*dSource.data['width'][I],\
                height=i*dSource.data['height'][I],\
                angle=dSource.data['angle'][I],
                line_color=COLORS[-100],
                fill_color=None)
    p.xaxis.axis_label = 'x Coordinate'
    p.yaxis.axis_label = 'y Coordinate'
    return p

def temporalTrace(I, TEMP, plth, pltw):
    temp = TEMP[:, :, I]
    p = figure(plot_height=plth, plot_width=pltw, title='Unit Mean Waveform', tools='pan,box_zoom,reset',\
            y_range=Y_RANGE)
    for i in range(0, nchan):
        p.line(range(0, nchan),\
                temp[:, i],\
                line_color=COLORS[4*i])
    p.xaxis.axis_label = 'time sample'
    p.yaxis.axis_label = 'voltage' 
    return p

def fullTrace(dSource, plth, pltw, _con):
    hover = HoverTool(\
            tooltips=[("index", "$index"),\
            ("desc", "@desc"),\
            ]\
            )
    sP = figure(plot_height=plth, plot_width=pltw, title='Spatial Covariance of all Units', tools = [hover, 'tap', 'reset']) 
    sP.scatter(GEOM[:, 0], GEOM[:, 1], radius=1, color=COLORS[0], line_color=None)
    sP.circle('x', 'y', color='color', size=5, source=dSource)
    if _con == 'energy':
        sP.oval('x',\
                'y',\
                width='width',
                height='height',
                angle='angle',
                color='color',
                fill_alpha=0.4,\
                source=dSource
                )
    else:
        sP.oval('x',\
                'y',\
                width='width',
                height='height',
                angle='angle',
                color='colortp',
                fill_alpha=0.4,\
                source=dSource
                )

    sP.xaxis.axis_label = 'x Coordinate'
    sP.yaxis.axis_label = 'y Coordinate'
    return sP


######
def create_figure(plth=450, pltw=450):
    I = int(unit.value)
    _con = contrast.value

    s1 = spatialTrace(I, TEMP_GT, gtSource, plth, pltw)
    ## draw temporal voltage trace
    s2 = temporalTrace(I, TEMP_GT, plth, pltw) 
    ## all spatial traces 
    s3 = fullTrace(gtSource, plth, pltw, _con) 

    s4 = spatialTrace(K_MATCH[I]-1, TEMP_DT, dtSource, plth, pltw)
    ## draw temporal voltage trace
    s5 = temporalTrace(K_MATCH[I]-1, TEMP_DT, plth, pltw) 
    ## all spatial traces 
    s6 = fullTrace(dtSource, plth, pltw, _con)

    output_file("test.html", title="example")
    p = gridplot([[s1, s2, s3], [s4, s5, s6]])
    #p = gridplot([[s1, s2, s3]])
    return p


def update(attr, old, new):
    layout.children[1] = create_figure()

unit = Select(title='Unit', value='0', options=OPT)
unit.on_change('value', update)

contrast = Select(title='Contrast', value='energy', options=['energy', 'accuracy'])
contrast.on_change('value', update)

#controls = widgetbox([x, y, color, size], width=200)
controls = widgetbox([unit, contrast], width=200)

layout = row(controls, create_figure())

curdoc().add_root(layout)
curdoc().title = "Electrode Array and Cells Geometry"
