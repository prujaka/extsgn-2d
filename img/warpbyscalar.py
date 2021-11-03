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
warpByScalar1.ScaleFactor = 200.0


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

# # create a new 'Plot Over Line'
# t1 = time.time()
# print('creating a new Plot Over Line...')
# plotOverLine1 = PlotOverLine(Input=cellDatatoPointData1,
#     Source='High Resolution Line Source')
#
# # init the 'High Resolution Line Source' selected for 'Source'
# plotOverLine1.Source.Point1 = [0.0, 300.0, 0.0]
# plotOverLine1.Source.Point2 = [600.0, 300.0, 0.0]
#
# # find source
# gradientOfUnstructuredDataSet1 = FindSource('GradientOfUnstructuredDataSet1')
#
# # show data in view
# plotOverLine1Display = Show(plotOverLine1, renderView1)
#
# # trace defaults for the display properties.
# plotOverLine1Display.Representation = 'Surface'
# plotOverLine1Display.ColorArrayName = ['POINTS', 'depth']
# plotOverLine1Display.LookupTable = depthLUT
# plotOverLine1Display.OSPRayScaleArray = 'depth'
# plotOverLine1Display.OSPRayScaleFunction = 'PiecewiseFunction'
# plotOverLine1Display.SelectOrientationVectors = 'velocity(m/s)'
# plotOverLine1Display.ScaleFactor = 4.0
# plotOverLine1Display.SelectScaleArray = 'depth'
# plotOverLine1Display.GlyphType = 'Arrow'
# plotOverLine1Display.GlyphTableIndexArray = 'depth'
# plotOverLine1Display.GaussianRadius = 0.2
# plotOverLine1Display.SetScaleArray = ['POINTS', 'depth']
# plotOverLine1Display.ScaleTransferFunction = 'PiecewiseFunction'
# plotOverLine1Display.OpacityArray = ['POINTS', 'depth']
# plotOverLine1Display.OpacityTransferFunction = 'PiecewiseFunction'
# plotOverLine1Display.DataAxesGrid = 'GridAxesRepresentation'
# plotOverLine1Display.SelectionCellLabelFontFile = ''
# plotOverLine1Display.SelectionPointLabelFontFile = ''
# plotOverLine1Display.PolarAxes = 'PolarAxesRepresentation'
#
# # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
# plotOverLine1Display.DataAxesGrid.XTitleFontFile = ''
# plotOverLine1Display.DataAxesGrid.YTitleFontFile = ''
# plotOverLine1Display.DataAxesGrid.ZTitleFontFile = ''
# plotOverLine1Display.DataAxesGrid.XLabelFontFile = ''
# plotOverLine1Display.DataAxesGrid.YLabelFontFile = ''
# plotOverLine1Display.DataAxesGrid.ZLabelFontFile = ''
#
# # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
# plotOverLine1Display.PolarAxes.PolarAxisTitleFontFile = ''
# plotOverLine1Display.PolarAxes.PolarAxisLabelFontFile = ''
# plotOverLine1Display.PolarAxes.LastRadialAxisTextFontFile = ''
# plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
#
# # Create a new 'Line Chart View'
# lineChartView1 = CreateView('XYChartView')
# lineChartView1.ViewSize = [635, 482]
# lineChartView1.ChartTitleFontFile = ''
# lineChartView1.LeftAxisTitleFontFile = ''
# lineChartView1.LeftAxisRangeMaximum = 6.66
# lineChartView1.LeftAxisLabelFontFile = ''
# lineChartView1.BottomAxisTitleFontFile = ''
# lineChartView1.BottomAxisRangeMaximum = 6.66
# lineChartView1.BottomAxisLabelFontFile = ''
# lineChartView1.RightAxisRangeMaximum = 6.66
# lineChartView1.RightAxisLabelFontFile = ''
# lineChartView1.TopAxisTitleFontFile = ''
# lineChartView1.TopAxisRangeMaximum = 6.66
# lineChartView1.TopAxisLabelFontFile = ''
#
# # get layout
# layout1 = GetLayout()
#
# # place view in the layout
# layout1.AssignView(2, lineChartView1)
#
# # show data in view
# plotOverLine1Display_1 = Show(plotOverLine1, lineChartView1)
#
# # trace defaults for the display properties.
# plotOverLine1Display_1.CompositeDataSetIndex = [0]
# plotOverLine1Display_1.UseIndexForXAxis = 0
# plotOverLine1Display_1.XArrayName = 'arc_length'
# plotOverLine1Display_1.SeriesVisibility = ['depth']
# plotOverLine1Display_1.SeriesLabel = ['arc_length', 'arc_length', 'depth', 'depth', 'pressure', 'pressure', 'sound_speed', 'sound_speed', 'velocity(m/s)_X', 'velocity(m/s)_X', 'velocity(m/s)_Y', 'velocity(m/s)_Y', 'velocity(m/s)_Z', 'velocity(m/s)_Z', 'velocity(m/s)_Magnitude', 'velocity(m/s)_Magnitude', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
# plotOverLine1Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'depth', '0.89', '0.1', '0.11', 'pressure', '0.22', '0.49', '0.72', 'sound_speed', '0.3', '0.69', '0.29', 'velocity(m/s)_X', '0.6', '0.31', '0.64', 'velocity(m/s)_Y', '1', '0.5', '0', 'velocity(m/s)_Z', '0.65', '0.34', '0.16', 'velocity(m/s)_Magnitude', '0', '0', '0', 'vtkValidPointMask', '0.89', '0.1', '0.11', 'Points_X', '0.22', '0.49', '0.72', 'Points_Y', '0.3', '0.69', '0.29', 'Points_Z', '0.6', '0.31', '0.64', 'Points_Magnitude', '1', '0.5', '0']
# plotOverLine1Display_1.SeriesPlotCorner = ['arc_length', '0', 'depth', '0', 'pressure', '0', 'sound_speed', '0', 'velocity(m/s)_X', '0', 'velocity(m/s)_Y', '0', 'velocity(m/s)_Z', '0', 'velocity(m/s)_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
# plotOverLine1Display_1.SeriesLabelPrefix = ''
# plotOverLine1Display_1.SeriesLineStyle = ['arc_length', '1', 'depth', '1', 'pressure', '1', 'sound_speed', '1', 'velocity(m/s)_X', '1', 'velocity(m/s)_Y', '1', 'velocity(m/s)_Z', '1', 'velocity(m/s)_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
# plotOverLine1Display_1.SeriesLineThickness = ['arc_length', '2', 'depth', '2', 'pressure', '2', 'sound_speed', '2', 'velocity(m/s)_X', '2', 'velocity(m/s)_Y', '2', 'velocity(m/s)_Z', '2', 'velocity(m/s)_Magnitude', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
# plotOverLine1Display_1.SeriesMarkerStyle = ['arc_length', '0', 'depth', '0', 'pressure', '0', 'sound_speed', '0', 'velocity(m/s)_X', '0', 'velocity(m/s)_Y', '0', 'velocity(m/s)_Z', '0', 'velocity(m/s)_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
#
# # update the view to ensure updated data information
# renderView1.Update()
#
# # update the view to ensure updated data information
# lineChartView1.Update()
# t2 = time.time()
# timeString = '{:.3f}'.format(t2 - t1)
# print('Plot Over Line: done\n')
#
# Hide(plotOverLine1, renderView1)

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

#renderView1.AxesGrid.AxesToLabel = 27

# just x and y in 2d mode
renderView1.AxesGrid.AxesToLabel = 3

renderView1.AxesGrid.ZAxisUseCustomLabels = 1
renderView1.AxesGrid.ZAxisLabels = [0.85, 1.0, 1.15]
renderView1.AxesGrid.DataScale = [1.0, 1.0, 400.0]
renderView1.AxesGrid.DataPosition = [300.0, 300.0, 0.0]

# current camera placement for renderView1
renderView1.CameraPosition = [1489.7949448496136, 849.3328587722575, 1199.8204223409357]
renderView1.CameraFocalPoint = [0.0, 0.0, 203.3050537109375]
renderView1.CameraViewUp = [-0.5380390098292285, -0.27812697283309873, 0.7957131461052283]
renderView1.CameraParallelScale = 426.10305291966154

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
renderView1.AxesGrid.XTitle = 'x, m'
renderView1.AxesGrid.YTitle = 'y, m    '
renderView1.AxesGrid.ZTitle = 'h, m'
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

renderView1.CameraPosition = [300.0, 300.0, 1849.6408276138122]
renderView1.CameraFocalPoint = [300.0, 300.0, 203.3050537109375]
renderView1.CameraParallelScale = 388.5759152433676
#END Z 2D camera

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

tFinish = time.time()
timeString = '{:.3f}'.format(tFinish - tStart)
print('Overall time: ' + timeString + ' sec\n')

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
