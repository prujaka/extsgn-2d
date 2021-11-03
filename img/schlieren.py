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

# create a new 'Gradient Of Unstructured DataSet'
t1 = time.time()
print('creating a new Gradient...')
gradient1 = Gradient(Input=cellDatatoPointData1)
gradient1.SelectInputScalars = ['POINTS', 'depth']

# show data in view
gradient1Display = Show(gradient1, renderView1)

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# update the view to ensure updated data information
renderView1.Update()
t2 = time.time()
timeString = '{:.3f}'.format(t2 - t1)
print('Gradient: done\n')

# create a new 'Calculator'
t1 = time.time()
print('creating a new Calculator...')
calculator1 = Calculator(Input=gradient1)

# Properties modified on calculator1
calculator1.ResultArrayName = 'Schlieren'

# Properties modified on calculator1
calculator1.Function = 'ln(1+2*sqrt((depthGradient_X)^2+(depthGradient_Y)^2+(depthGradient_Z)^2))'

# show data in view
calculator1Display = Show(calculator1, renderView1)

# get color transfer function/color map for 'Schlieren'
schlierenLUT = GetColorTransferFunction('Schlieren')

# get opacity transfer function/opacity map for 'Schlieren'
resultPWF = GetOpacityTransferFunction('Schlieren')

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['POINTS', 'Schlieren']
calculator1Display.LookupTable = schlierenLUT
calculator1Display.OSPRayScaleArray = 'Schlieren'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'velocity(m/s)'
calculator1Display.ScaleFactor = 60.0
calculator1Display.SelectScaleArray = 'Schlieren'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'Schlieren'
calculator1Display.GaussianRadius = 3.0
calculator1Display.SetScaleArray = ['POINTS', 'Schlieren']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'Schlieren']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.SelectionCellLabelFontFile = ''
calculator1Display.SelectionPointLabelFontFile = ''
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.ScalarOpacityFunction = resultPWF
calculator1Display.ScalarOpacityUnitDistance = 39.38518727672853

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
calculator1Display.DataAxesGrid.XTitleFontFile = ''
calculator1Display.DataAxesGrid.YTitleFontFile = ''
calculator1Display.DataAxesGrid.ZTitleFontFile = ''
calculator1Display.DataAxesGrid.XLabelFontFile = ''
calculator1Display.DataAxesGrid.YLabelFontFile = ''
calculator1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
calculator1Display.PolarAxes.PolarAxisTitleFontFile = ''
calculator1Display.PolarAxes.PolarAxisLabelFontFile = ''
calculator1Display.PolarAxes.LastRadialAxisTextFontFile = ''
calculator1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(gradient1, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()
t2 = time.time()
timeString = '{:.3f}'.format(t2 - t1)
print('Calculator: done\n')

# invert the transfer function
schlierenLUT.InvertTransferFunction()

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
schlierenLUT.ApplyPreset('X Ray', True)

#### saving camera placements for all active views

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1279, 482]

# reset view to fit data
renderView1.ResetCamera()

# find source
cellDatatoPointData1 = FindSource('CellDatatoPointData1')

# set active source
SetActiveSource(cellDatatoPointData1)

# get color transfer function/color map for 'depth'
depthLUT = GetColorTransferFunction('depth')

# get opacity transfer function/opacity map for 'depth'
depthPWF = GetOpacityTransferFunction('depth')

## create a new 'Plot Over Line'
#t1 = time.time()
#print('creating a new Plot Over Line...')
#plotOverLine1 = PlotOverLine(Input=cellDatatoPointData1,
#    Source='High Resolution Line Source')

## init the 'High Resolution Line Source' selected for 'Source'
#plotOverLine1.Source.Point1 = [0.0, 300.0, 0.0]
#plotOverLine1.Source.Point2 = [600.0, 300.0, 0.0]

## find source
#gradientOfUnstructuredDataSet1 = FindSource('GradientOfUnstructuredDataSet1')

## show data in view
#plotOverLine1Display = Show(plotOverLine1, renderView1)

## trace defaults for the display properties.
#plotOverLine1Display.Representation = 'Surface'
#plotOverLine1Display.ColorArrayName = ['POINTS', 'depth']
#plotOverLine1Display.LookupTable = depthLUT
#plotOverLine1Display.OSPRayScaleArray = 'depth'
#plotOverLine1Display.OSPRayScaleFunction = 'PiecewiseFunction'
#plotOverLine1Display.SelectOrientationVectors = 'velocity(m/s)'
#plotOverLine1Display.ScaleFactor = 4.0
#plotOverLine1Display.SelectScaleArray = 'depth'
#plotOverLine1Display.GlyphType = 'Arrow'
#plotOverLine1Display.GlyphTableIndexArray = 'depth'
#plotOverLine1Display.GaussianRadius = 0.2
#plotOverLine1Display.SetScaleArray = ['POINTS', 'depth']
#plotOverLine1Display.ScaleTransferFunction = 'PiecewiseFunction'
#plotOverLine1Display.OpacityArray = ['POINTS', 'depth']
#plotOverLine1Display.OpacityTransferFunction = 'PiecewiseFunction'
#plotOverLine1Display.DataAxesGrid = 'GridAxesRepresentation'
#plotOverLine1Display.SelectionCellLabelFontFile = ''
#plotOverLine1Display.SelectionPointLabelFontFile = ''
#plotOverLine1Display.PolarAxes = 'PolarAxesRepresentation'

## init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
#plotOverLine1Display.DataAxesGrid.XTitleFontFile = ''
#plotOverLine1Display.DataAxesGrid.YTitleFontFile = ''
#plotOverLine1Display.DataAxesGrid.ZTitleFontFile = ''
#plotOverLine1Display.DataAxesGrid.XLabelFontFile = ''
#plotOverLine1Display.DataAxesGrid.YLabelFontFile = ''
#plotOverLine1Display.DataAxesGrid.ZLabelFontFile = ''

## init the 'PolarAxesRepresentation' selected for 'PolarAxes'
#plotOverLine1Display.PolarAxes.PolarAxisTitleFontFile = ''
#plotOverLine1Display.PolarAxes.PolarAxisLabelFontFile = ''
#plotOverLine1Display.PolarAxes.LastRadialAxisTextFontFile = ''
#plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

## Create a new 'Line Chart View'
#lineChartView1 = CreateView('XYChartView')
#lineChartView1.ViewSize = [635, 482]
#lineChartView1.ChartTitleFontFile = ''
#lineChartView1.LeftAxisTitleFontFile = ''
#lineChartView1.LeftAxisRangeMaximum = 6.66
#lineChartView1.LeftAxisLabelFontFile = ''
#lineChartView1.BottomAxisTitleFontFile = ''
#lineChartView1.BottomAxisRangeMaximum = 6.66
#lineChartView1.BottomAxisLabelFontFile = ''
#lineChartView1.RightAxisRangeMaximum = 6.66
#lineChartView1.RightAxisLabelFontFile = ''
#lineChartView1.TopAxisTitleFontFile = ''
#lineChartView1.TopAxisRangeMaximum = 6.66
#lineChartView1.TopAxisLabelFontFile = ''

## Properties modified on lineChartView1
#lineChartView1.LeftAxisUseCustomRange = 1
## Properties modified on lineChartView1
#lineChartView1.LeftAxisRangeMinimum = 0.8
## Properties modified on lineChartView1
#lineChartView1.LeftAxisRangeMaximum = 1.81
## Properties modified on lineChartView1
##lineChartView1.BottomAxisUseCustomRange = 1

## get layout
#layout1 = GetLayout()

## place view in the layout
#layout1.AssignView(2, lineChartView1)

## show data in view
#plotOverLine1Display_1 = Show(plotOverLine1, lineChartView1)

## trace defaults for the display properties.
#plotOverLine1Display_1.CompositeDataSetIndex = [0]
#plotOverLine1Display_1.UseIndexForXAxis = 0
#plotOverLine1Display_1.XArrayName = 'arc_length'
#plotOverLine1Display_1.SeriesVisibility = ['depth']
#plotOverLine1Display_1.SeriesLabel = ['arc_length', 'arc_length', 'depth', 'depth', 'pressure', 'pressure', 'sound_speed', 'sound_speed', 'velocity(m/s)_X', 'velocity(m/s)_X', 'velocity(m/s)_Y', 'velocity(m/s)_Y', 'velocity(m/s)_Z', 'velocity(m/s)_Z', 'velocity(m/s)_Magnitude', 'velocity(m/s)_Magnitude', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
#plotOverLine1Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'depth', '0.89', '0.1', '0.11', 'pressure', '0.22', '0.49', '0.72', 'sound_speed', '0.3', '0.69', '0.29', 'velocity(m/s)_X', '0.6', '0.31', '0.64', 'velocity(m/s)_Y', '1', '0.5', '0', 'velocity(m/s)_Z', '0.65', '0.34', '0.16', 'velocity(m/s)_Magnitude', '0', '0', '0', 'vtkValidPointMask', '0.89', '0.1', '0.11', 'Points_X', '0.22', '0.49', '0.72', 'Points_Y', '0.3', '0.69', '0.29', 'Points_Z', '0.6', '0.31', '0.64', 'Points_Magnitude', '1', '0.5', '0']
#plotOverLine1Display_1.SeriesPlotCorner = ['arc_length', '0', 'depth', '0', 'pressure', '0', 'sound_speed', '0', 'velocity(m/s)_X', '0', 'velocity(m/s)_Y', '0', 'velocity(m/s)_Z', '0', 'velocity(m/s)_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
#plotOverLine1Display_1.SeriesLabelPrefix = ''
#plotOverLine1Display_1.SeriesLineStyle = ['arc_length', '1', 'depth', '1', 'pressure', '1', 'sound_speed', '1', 'velocity(m/s)_X', '1', 'velocity(m/s)_Y', '1', 'velocity(m/s)_Z', '1', 'velocity(m/s)_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
#plotOverLine1Display_1.SeriesLineThickness = ['arc_length', '2', 'depth', '2', 'pressure', '2', 'sound_speed', '2', 'velocity(m/s)_X', '2', 'velocity(m/s)_Y', '2', 'velocity(m/s)_Z', '2', 'velocity(m/s)_Magnitude', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
#plotOverLine1Display_1.SeriesMarkerStyle = ['arc_length', '0', 'depth', '0', 'pressure', '0', 'sound_speed', '0', 'velocity(m/s)_X', '0', 'velocity(m/s)_Y', '0', 'velocity(m/s)_Z', '0', 'velocity(m/s)_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

## find source
#mergeBlocks1 = FindSource('MergeBlocks1')

## find source
#calculator1 = FindSource('Calculator1')

## find source
#groupDatasets1 = FindSource('GroupDatasets1')

## update the view to ensure updated data information
#renderView1.Update()

## update the view to ensure updated data information
#lineChartView1.Update()
#t2 = time.time()
#timeString = '{:.3f}'.format(t2 - t1)
#print('Plot Over Line: done\n')

## find view
#renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
## uncomment following to set a specific view size
## renderView1.ViewSize = [613, 492]

#Hide(plotOverLine1, renderView1)

## set active view
#SetActiveView(renderView1)

## get color transfer function/color map for 'depth'
#depthLUT = GetColorTransferFunction('depth')

## get opacity transfer function/opacity map for 'depth'
#depthPWF = GetOpacityTransferFunction('depth')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [613, 492]

# FANCY OUTPUT WITH AXES
renderView1.OrientationAxesVisibility = 0

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XTitle = '\n x, m'
renderView1.AxesGrid.YTitle = 'y, m         '
renderView1.AxesGrid.XTitleFontSize = 50
renderView1.AxesGrid.YTitleFontSize = 50
renderView1.AxesGrid.ZTitleFontSize = 50
renderView1.AxesGrid.AxesToLabel = 3
renderView1.AxesGrid.XLabelFontSize = 50
renderView1.AxesGrid.YLabelFontSize = 50
renderView1.AxesGrid.ZLabelFontSize = 50

#renderView1.AxesGrid.XAxisUseCustomLabels = 1
#renderView1.AxesGrid.XAxisLabels = [0, 150, 300, 450, 600]
renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisLabels = [-300, -200, -100, 0, 100, 200, 300]

renderView1.AxesGrid.DataPosition = [300, 300, 0]
renderView1.AxesGrid.DataBoundsScaleFactor = 1.0


# get color transfer function/color map for 'Schlieren'
schlierenLUT = GetColorTransferFunction('Schlieren')

# get color legend/bar for schlierenLUT in view renderView1
schlierenLUTColorBar = GetScalarBar(schlierenLUT, renderView1)

# Properties modified on schlierenLUTColorBar
schlierenLUTColorBar.WindowLocation = 'AnyLocation'
schlierenLUTColorBar.Position = [0.85, 0.15]
schlierenLUTColorBar.Title = """f(\\nabla h)"""
schlierenLUTColorBar.HorizontalTitle = 1
schlierenLUTColorBar.AddRangeLabels = 0
schlierenLUTColorBar.TitleFontSize = 12
schlierenLUTColorBar.LabelFontSize = 12
schlierenLUTColorBar.ScalarBarThickness = 12
schlierenLUTColorBar.ScalarBarLength = 0.4

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# get display properties
calculator1Display = GetDisplayProperties(calculator1, view=renderView1)

# hide color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, False)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [250.0, 280.0, 1300.0]
renderView1.CameraFocalPoint = [250.0, 280.0, 5e-05]
renderView1.CameraParallelScale = 600

#Hide(schlierenLUTColorBar, renderView1)

# get display properties
calculator1Display = GetDisplayProperties(calculator1, view=renderView1)

# rescale color and/or opacity maps used to exactly fit the current data range
calculator1Display.RescaleTransferFunctionToDataRange(False, True)

# set active source
SetActiveSource(calculator1)

tFinish = time.time()
timeString = '{:.3f}'.format(tFinish - tStart)
print('Overall time: ' + timeString + ' sec\n')

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
