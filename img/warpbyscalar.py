"""
paraview script to visualize the water surface in vtk format
file paths are hard-coded specially for the lab working machine
warp by scalar is stretched x300 times in z direction

"""

# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

# import the simple module from the paraview
from paraview.simple import *
# import numpy as np
import os
import time

ResetSession()
tStart = time.time()
LoadPalette(paletteName='WhiteBackground')

filePathAbsolute = '/Users/prujaka/work/codes/imex2d-extsgn/res.vtk'

print(filePathAbsolute)

# create a new 'Legacy VTK Reader'
t1 = time.time()
print('creating a new Legacy VTK Reader...')
globvtk = LegacyVTKReader(FileNames=[filePathAbsolute])


# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
## uncomment following to set a specific view size
## renderView1.ViewSize = [1257, 461]

# update the view to ensure updated data information
renderView1.Update()
t2 = time.time()
timeString = '{:.3f}'.format(t2 - t1)
print('Legacy VTK Reader: done\n')


# create a new 'Cell Data to Point Data'
t1 = time.time()
print('creating a new Cell Data to Point Data...')
cellDatatoPointData1 = CellDatatoPointData(Input=globvtk)

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1)

# hide data in view
Hide(globvtk, renderView1)

# show color bar/color legend
#cellDatatoPointData1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()
t2 = time.time()
timeString = '{:.3f}'.format(t2 - t1)
print('Cell Data to Point Data: done\n')

# find source
cellDatatoPointData1 = FindSource('CellDatatoPointData1')

# create a new 'Warp By Scalar'
t1 = time.time()
print('creating a new Warp By Scalar...')
warpByScalar1 = WarpByScalar(Input=cellDatatoPointData1)
warpByScalar1.Scalars = ['POINTS', 'depth']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1311, 506]

# show data in view
warpByScalar1Display = Show(warpByScalar1, renderView1)

# trace defaults for the display properties.
warpByScalar1Display.Representation = 'Outline'
warpByScalar1Display.AmbientColor = [0.0, 0.0, 0.0]
warpByScalar1Display.ColorArrayName = ['POINTS', '']
warpByScalar1Display.OSPRayScaleArray = 'depth'
warpByScalar1Display.OSPRayScaleFunction = 'PiecewiseFunction'
warpByScalar1Display.SelectOrientationVectors = 'velocity(m/s)'
warpByScalar1Display.ScaleFactor = 60.0
warpByScalar1Display.SelectScaleArray = 'depth'
warpByScalar1Display.GlyphType = 'Arrow'
warpByScalar1Display.GlyphTableIndexArray = 'depth'
warpByScalar1Display.GaussianRadius = 3.0
warpByScalar1Display.SetScaleArray = ['POINTS', 'depth']
warpByScalar1Display.ScaleTransferFunction = 'PiecewiseFunction'
warpByScalar1Display.OpacityArray = ['POINTS', 'depth']
warpByScalar1Display.OpacityTransferFunction = 'PiecewiseFunction'
warpByScalar1Display.DataAxesGrid = 'GridAxesRepresentation'
warpByScalar1Display.SelectionCellLabelFontFile = ''
warpByScalar1Display.SelectionPointLabelFontFile = ''
warpByScalar1Display.PolarAxes = 'PolarAxesRepresentation'
warpByScalar1Display.ScalarOpacityUnitDistance = 8.508363702731133

# get active source.
warpByScalar1 = GetActiveSource()

# Properties modified on warpByScalar1
warpByScalar1.ScaleFactor = 25.0


# hide data in view
Hide(cellDatatoPointData1, renderView1)

# find source
legacyVTKReader1 = FindSource('LegacyVTKReader1')

# update the view to ensure updated data information
renderView1.Update()

# change representation type
warpByScalar1Display.SetRepresentationType('Surface')
t2 = time.time()
timeString = '{:.3f}'.format(t2 - t1)
print('Warp By Scalar: done\n')


# find source
cellDatatoPointData1 = FindSource('CellDatatoPointData1')

# set active source
SetActiveSource(cellDatatoPointData1)

# get color transfer function/color map for 'depth'
depthLUT = GetColorTransferFunction('depth')

# get opacity transfer function/opacity map for 'depth'
depthPWF = GetOpacityTransferFunction('depth')

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [651, 506]

# set active view
SetActiveView(renderView1)
# set active source
SetActiveSource(warpByScalar1)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [651, 506]

# Properties modified on renderView1
renderView1.OrientationAxesVisibility = 0

# get the material library
materialLibrary1 = GetMaterialLibrary()

## Properties modified on renderView1.AxesGrid
#renderView1.AxesGrid.XTitle = ''
#renderView1.AxesGrid.YTitle = ''
#renderView1.AxesGrid.ZTitle = 'h, m'


# Properties modified on renderView1.AxesGrid
#renderView1.AxesGrid.FacesToRender = 62
renderView1.AxesGrid.FacesToRender = 36


## PNG OUTPUT FONT SIZE
fontsize = 24
## Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XTitleFontSize = fontsize
renderView1.AxesGrid.YTitleFontSize = fontsize
renderView1.AxesGrid.ZTitleFontSize = fontsize

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XLabelFontSize = fontsize
renderView1.AxesGrid.YLabelFontSize = fontsize
renderView1.AxesGrid.ZLabelFontSize = fontsize
# END PNG OUTPUT FONT SIZE

# Z 2D camera
# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XTitle = ''
renderView1.AxesGrid.YTitle = ''
renderView1.AxesGrid.ZTitle = ''
renderView1.InteractionMode = '2D'

renderView1.AxesGrid.FacesToRender = 36

paraview.simple._DisableFirstRenderCameraReset()

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2154, 642]

# get layout
layout1 = GetLayoutByName("Layout #1")

# set active view
SetActiveView(renderView1)

# renderView1.CameraPosition = [300.0, 300.0, 1849.6408276138122]
# renderView1.CameraFocalPoint = [300.0, 300.0, 203.3050537109375]
# renderView1.CameraParallelScale = 388.5759152433676

# get layout
layout1 = GetLayout()

# get display properties
warpByScalar1Display = GetDisplayProperties(warpByScalar1, view=renderView1)

# set scalar coloring
ColorBy(warpByScalar1Display, ('POINTS', 'depth'))

# get color transfer function/color map for 'depth'
depthLUT = GetColorTransferFunction('depth')

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(depthLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
warpByScalar1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
warpByScalar1Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'depth'
depthPWF = GetOpacityTransferFunction('depth')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
depthLUT.ApplyPreset('Linear Blue (8_31f)', True)

# get active source.
warpByScalar1 = GetActiveSource()

# get display properties
warpByScalar1Display = GetDisplayProperties(warpByScalar1, view=renderView1)

# hide color bar/color legend
warpByScalar1Display.SetScalarBarVisibility(renderView1, False)

renderView1.InteractionMode = '3D'
renderView1.CameraPosition = [300.0, 300.0, 1694.6291424366293]
renderView1.CameraFocalPoint = [300.0, 300.0, 54.98125076293945]
renderView1.CameraParallelScale = 424.3721016273458


#END Z 2D camera

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 0

tFinish = time.time()
timeString = '{:.3f}'.format(tFinish - tStart)
print('Overall time: ' + timeString + ' sec\n')

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
