# Alireza Rashti
# January 2021

# A VisIt script to visualize patch, scalar, and vector in Elliptica.
# How to invoke the script:
# invoke:
# $ visit_plot.py -h



# import packages:
from __future__ import division
import re
import os
import argparse
from argparse import RawTextHelpFormatter

## data file suffix
FILE_SUFFIX = "silo"

## quick options:
# turn on random mesh coloring
opt_random_mesh_color = 1
opt_annotation        = 1

##{ main
def main():

  # get args
  directory_path,file_names,scalar_name,vector_name = pars_arguments()
  
  # deleting all plots
  DeleteAllPlots()
  
  plot(directory_path,file_names)
  
##} main

##{ plot: given object, data and the section of grid plot the object
def plot(directory_path,file_names):
  
  # collecting all of the files name needed to load
  files_name = collect_files(directory_path,file_names)
  
  n = len(files_name)
  
  for i in range(n):
    print('Openning data file: ' + files_name[i])
  
  SetWindowLayout(1)
  
  SetActiveWindow(1)
  for i in range(0,n):
    OpenDatabase(files_name[i], 0)
    mesh_name = re.sub(r'_xyz_.*\.silo','',files_name[i])
    AddPlot('Mesh',mesh_name, 1, 1)
    DrawPlots()
    # turn off annotation
    if (opt_annotation == 0):
      AnnotationAtts = AnnotationAttributes()
      AnnotationAtts.axes2D.visible = 0
      AnnotationAtts.axes2D.autoSetTicks = 1
      AnnotationAtts.axes2D.autoSetScaling = 1
      AnnotationAtts.axes2D.lineWidth = 0
      AnnotationAtts.axes2D.tickLocation = AnnotationAtts.axes2D.Outside  # Inside, Outside, Both
      AnnotationAtts.axes2D.tickAxes = AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
      AnnotationAtts.axes2D.xAxis.title.visible = 1
      AnnotationAtts.axes2D.xAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
      AnnotationAtts.axes2D.xAxis.title.font.scale = 1
      AnnotationAtts.axes2D.xAxis.title.font.useForegroundColor = 1
      AnnotationAtts.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes2D.xAxis.title.font.bold = 1
      AnnotationAtts.axes2D.xAxis.title.font.italic = 1
      AnnotationAtts.axes2D.xAxis.title.userTitle = 0
      AnnotationAtts.axes2D.xAxis.title.userUnits = 0
      AnnotationAtts.axes2D.xAxis.title.title = "X-Axis"
      AnnotationAtts.axes2D.xAxis.title.units = ""
      AnnotationAtts.axes2D.xAxis.label.visible = 1
      AnnotationAtts.axes2D.xAxis.label.font.font = AnnotationAtts.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
      AnnotationAtts.axes2D.xAxis.label.font.scale = 1
      AnnotationAtts.axes2D.xAxis.label.font.useForegroundColor = 1
      AnnotationAtts.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes2D.xAxis.label.font.bold = 1
      AnnotationAtts.axes2D.xAxis.label.font.italic = 1
      AnnotationAtts.axes2D.xAxis.label.scaling = 0
      AnnotationAtts.axes2D.xAxis.tickMarks.visible = 1
      AnnotationAtts.axes2D.xAxis.tickMarks.majorMinimum = 0
      AnnotationAtts.axes2D.xAxis.tickMarks.majorMaximum = 1
      AnnotationAtts.axes2D.xAxis.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axes2D.xAxis.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axes2D.xAxis.grid = 0
      AnnotationAtts.axes2D.yAxis.title.visible = 1
      AnnotationAtts.axes2D.yAxis.title.font.font = AnnotationAtts.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
      AnnotationAtts.axes2D.yAxis.title.font.scale = 1
      AnnotationAtts.axes2D.yAxis.title.font.useForegroundColor = 1
      AnnotationAtts.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes2D.yAxis.title.font.bold = 1
      AnnotationAtts.axes2D.yAxis.title.font.italic = 1
      AnnotationAtts.axes2D.yAxis.title.userTitle = 0
      AnnotationAtts.axes2D.yAxis.title.userUnits = 0
      AnnotationAtts.axes2D.yAxis.title.title = "Y-Axis"
      AnnotationAtts.axes2D.yAxis.title.units = ""
      AnnotationAtts.axes2D.yAxis.label.visible = 1
      AnnotationAtts.axes2D.yAxis.label.font.font = AnnotationAtts.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
      AnnotationAtts.axes2D.yAxis.label.font.scale = 1
      AnnotationAtts.axes2D.yAxis.label.font.useForegroundColor = 1
      AnnotationAtts.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes2D.yAxis.label.font.bold = 1
      AnnotationAtts.axes2D.yAxis.label.font.italic = 1
      AnnotationAtts.axes2D.yAxis.label.scaling = 0
      AnnotationAtts.axes2D.yAxis.tickMarks.visible = 1
      AnnotationAtts.axes2D.yAxis.tickMarks.majorMinimum = 0
      AnnotationAtts.axes2D.yAxis.tickMarks.majorMaximum = 1
      AnnotationAtts.axes2D.yAxis.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axes2D.yAxis.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axes2D.yAxis.grid = 0
      AnnotationAtts.axes3D.visible = 0
      AnnotationAtts.axes3D.autoSetTicks = 1
      AnnotationAtts.axes3D.autoSetScaling = 1
      AnnotationAtts.axes3D.lineWidth = 0
      AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
      AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.ClosestTriad  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
      AnnotationAtts.axes3D.triadFlag = 0
      AnnotationAtts.axes3D.bboxFlag = 0
      AnnotationAtts.axes3D.xAxis.title.visible = 1
      AnnotationAtts.axes3D.xAxis.title.font.font = AnnotationAtts.axes3D.xAxis.title.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.xAxis.title.font.scale = 1
      AnnotationAtts.axes3D.xAxis.title.font.useForegroundColor = 1
      AnnotationAtts.axes3D.xAxis.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.xAxis.title.font.bold = 0
      AnnotationAtts.axes3D.xAxis.title.font.italic = 0
      AnnotationAtts.axes3D.xAxis.title.userTitle = 0
      AnnotationAtts.axes3D.xAxis.title.userUnits = 0
      AnnotationAtts.axes3D.xAxis.title.title = "X-Axis"
      AnnotationAtts.axes3D.xAxis.title.units = ""
      AnnotationAtts.axes3D.xAxis.label.visible = 1
      AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.xAxis.label.font.scale = 1
      AnnotationAtts.axes3D.xAxis.label.font.useForegroundColor = 1
      AnnotationAtts.axes3D.xAxis.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.xAxis.label.font.bold = 0
      AnnotationAtts.axes3D.xAxis.label.font.italic = 0
      AnnotationAtts.axes3D.xAxis.label.scaling = 0
      AnnotationAtts.axes3D.xAxis.tickMarks.visible = 1
      AnnotationAtts.axes3D.xAxis.tickMarks.majorMinimum = 0
      AnnotationAtts.axes3D.xAxis.tickMarks.majorMaximum = 1
      AnnotationAtts.axes3D.xAxis.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axes3D.xAxis.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axes3D.xAxis.grid = 0
      AnnotationAtts.axes3D.yAxis.title.visible = 1
      AnnotationAtts.axes3D.yAxis.title.font.font = AnnotationAtts.axes3D.yAxis.title.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.yAxis.title.font.scale = 1
      AnnotationAtts.axes3D.yAxis.title.font.useForegroundColor = 1
      AnnotationAtts.axes3D.yAxis.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.yAxis.title.font.bold = 0
      AnnotationAtts.axes3D.yAxis.title.font.italic = 0
      AnnotationAtts.axes3D.yAxis.title.userTitle = 0
      AnnotationAtts.axes3D.yAxis.title.userUnits = 0
      AnnotationAtts.axes3D.yAxis.title.title = "Y-Axis"
      AnnotationAtts.axes3D.yAxis.title.units = ""
      AnnotationAtts.axes3D.yAxis.label.visible = 1
      AnnotationAtts.axes3D.yAxis.label.font.font = AnnotationAtts.axes3D.yAxis.label.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.yAxis.label.font.scale = 1
      AnnotationAtts.axes3D.yAxis.label.font.useForegroundColor = 1
      AnnotationAtts.axes3D.yAxis.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.yAxis.label.font.bold = 0
      AnnotationAtts.axes3D.yAxis.label.font.italic = 0
      AnnotationAtts.axes3D.yAxis.label.scaling = 0
      AnnotationAtts.axes3D.yAxis.tickMarks.visible = 1
      AnnotationAtts.axes3D.yAxis.tickMarks.majorMinimum = 0
      AnnotationAtts.axes3D.yAxis.tickMarks.majorMaximum = 1
      AnnotationAtts.axes3D.yAxis.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axes3D.yAxis.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axes3D.yAxis.grid = 0
      AnnotationAtts.axes3D.zAxis.title.visible = 1
      AnnotationAtts.axes3D.zAxis.title.font.font = AnnotationAtts.axes3D.zAxis.title.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.zAxis.title.font.scale = 1
      AnnotationAtts.axes3D.zAxis.title.font.useForegroundColor = 1
      AnnotationAtts.axes3D.zAxis.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.zAxis.title.font.bold = 0
      AnnotationAtts.axes3D.zAxis.title.font.italic = 0
      AnnotationAtts.axes3D.zAxis.title.userTitle = 0
      AnnotationAtts.axes3D.zAxis.title.userUnits = 0
      AnnotationAtts.axes3D.zAxis.title.title = "Z-Axis"
      AnnotationAtts.axes3D.zAxis.title.units = ""
      AnnotationAtts.axes3D.zAxis.label.visible = 1
      AnnotationAtts.axes3D.zAxis.label.font.font = AnnotationAtts.axes3D.zAxis.label.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.zAxis.label.font.scale = 1
      AnnotationAtts.axes3D.zAxis.label.font.useForegroundColor = 1
      AnnotationAtts.axes3D.zAxis.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.zAxis.label.font.bold = 0
      AnnotationAtts.axes3D.zAxis.label.font.italic = 0
      AnnotationAtts.axes3D.zAxis.label.scaling = 0
      AnnotationAtts.axes3D.zAxis.tickMarks.visible = 1
      AnnotationAtts.axes3D.zAxis.tickMarks.majorMinimum = 0
      AnnotationAtts.axes3D.zAxis.tickMarks.majorMaximum = 1
      AnnotationAtts.axes3D.zAxis.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axes3D.zAxis.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axes3D.zAxis.grid = 0
      AnnotationAtts.axes3D.setBBoxLocation = 0
      AnnotationAtts.axes3D.bboxLocation = (0, 1, 0, 1, 0, 1)
      AnnotationAtts.axes3D.triadColor = (0, 0, 0)
      AnnotationAtts.axes3D.triadLineWidth = 0
      AnnotationAtts.axes3D.triadFont = 0
      AnnotationAtts.axes3D.triadBold = 1
      AnnotationAtts.axes3D.triadItalic = 1
      AnnotationAtts.axes3D.triadSetManually = 0
      AnnotationAtts.userInfoFlag = 0
      AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
      AnnotationAtts.userInfoFont.scale = 1
      AnnotationAtts.userInfoFont.useForegroundColor = 1
      AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
      AnnotationAtts.userInfoFont.bold = 0
      AnnotationAtts.userInfoFont.italic = 0
      AnnotationAtts.databaseInfoFlag = 0
      AnnotationAtts.timeInfoFlag = 1
      AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
      AnnotationAtts.databaseInfoFont.scale = 1
      AnnotationAtts.databaseInfoFont.useForegroundColor = 1
      AnnotationAtts.databaseInfoFont.color = (0, 0, 0, 255)
      AnnotationAtts.databaseInfoFont.bold = 0
      AnnotationAtts.databaseInfoFont.italic = 0
      AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
      AnnotationAtts.databaseInfoTimeScale = 1
      AnnotationAtts.databaseInfoTimeOffset = 0
      AnnotationAtts.legendInfoFlag = 0
      AnnotationAtts.backgroundColor = (255, 255, 255, 255)
      AnnotationAtts.foregroundColor = (0, 0, 0, 255)
      AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
      AnnotationAtts.gradientColor1 = (0, 0, 255, 255)
      AnnotationAtts.gradientColor2 = (0, 0, 0, 255)
      AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
      AnnotationAtts.backgroundImage = ""
      AnnotationAtts.imageRepeatX = 1
      AnnotationAtts.imageRepeatY = 1
      AnnotationAtts.axesArray.visible = 1
      AnnotationAtts.axesArray.ticksVisible = 1
      AnnotationAtts.axesArray.autoSetTicks = 1
      AnnotationAtts.axesArray.autoSetScaling = 1
      AnnotationAtts.axesArray.lineWidth = 0
      AnnotationAtts.axesArray.axes.title.visible = 1
      AnnotationAtts.axesArray.axes.title.font.font = AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axesArray.axes.title.font.scale = 1
      AnnotationAtts.axesArray.axes.title.font.useForegroundColor = 1
      AnnotationAtts.axesArray.axes.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axesArray.axes.title.font.bold = 0
      AnnotationAtts.axesArray.axes.title.font.italic = 0
      AnnotationAtts.axesArray.axes.title.userTitle = 0
      AnnotationAtts.axesArray.axes.title.userUnits = 0
      AnnotationAtts.axesArray.axes.title.title = ""
      AnnotationAtts.axesArray.axes.title.units = ""
      AnnotationAtts.axesArray.axes.label.visible = 1
      AnnotationAtts.axesArray.axes.label.font.font = AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axesArray.axes.label.font.scale = 1
      AnnotationAtts.axesArray.axes.label.font.useForegroundColor = 1
      AnnotationAtts.axesArray.axes.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axesArray.axes.label.font.bold = 0
      AnnotationAtts.axesArray.axes.label.font.italic = 0
      AnnotationAtts.axesArray.axes.label.scaling = 0
      AnnotationAtts.axesArray.axes.tickMarks.visible = 1
      AnnotationAtts.axesArray.axes.tickMarks.majorMinimum = 0
      AnnotationAtts.axesArray.axes.tickMarks.majorMaximum = 1
      AnnotationAtts.axesArray.axes.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axesArray.axes.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axesArray.axes.grid = 0
      SetAnnotationAttributes(AnnotationAtts)
      # Logging for SetAnnotationObjectOptions is not implemented yet.
      AnnotationAtts = AnnotationAttributes()
      AnnotationAtts.axes2D.visible = 0
      AnnotationAtts.axes2D.autoSetTicks = 1
      AnnotationAtts.axes2D.autoSetScaling = 1
      AnnotationAtts.axes2D.lineWidth = 0
      AnnotationAtts.axes2D.tickLocation = AnnotationAtts.axes2D.Outside  # Inside, Outside, Both
      AnnotationAtts.axes2D.tickAxes = AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
      AnnotationAtts.axes2D.xAxis.title.visible = 1
      AnnotationAtts.axes2D.xAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
      AnnotationAtts.axes2D.xAxis.title.font.scale = 1
      AnnotationAtts.axes2D.xAxis.title.font.useForegroundColor = 1
      AnnotationAtts.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes2D.xAxis.title.font.bold = 1
      AnnotationAtts.axes2D.xAxis.title.font.italic = 1
      AnnotationAtts.axes2D.xAxis.title.userTitle = 0
      AnnotationAtts.axes2D.xAxis.title.userUnits = 0
      AnnotationAtts.axes2D.xAxis.title.title = "X-Axis"
      AnnotationAtts.axes2D.xAxis.title.units = ""
      AnnotationAtts.axes2D.xAxis.label.visible = 1
      AnnotationAtts.axes2D.xAxis.label.font.font = AnnotationAtts.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
      AnnotationAtts.axes2D.xAxis.label.font.scale = 1
      AnnotationAtts.axes2D.xAxis.label.font.useForegroundColor = 1
      AnnotationAtts.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes2D.xAxis.label.font.bold = 1
      AnnotationAtts.axes2D.xAxis.label.font.italic = 1
      AnnotationAtts.axes2D.xAxis.label.scaling = 0
      AnnotationAtts.axes2D.xAxis.tickMarks.visible = 1
      AnnotationAtts.axes2D.xAxis.tickMarks.majorMinimum = 0
      AnnotationAtts.axes2D.xAxis.tickMarks.majorMaximum = 1
      AnnotationAtts.axes2D.xAxis.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axes2D.xAxis.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axes2D.xAxis.grid = 0
      AnnotationAtts.axes2D.yAxis.title.visible = 1
      AnnotationAtts.axes2D.yAxis.title.font.font = AnnotationAtts.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
      AnnotationAtts.axes2D.yAxis.title.font.scale = 1
      AnnotationAtts.axes2D.yAxis.title.font.useForegroundColor = 1
      AnnotationAtts.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes2D.yAxis.title.font.bold = 1
      AnnotationAtts.axes2D.yAxis.title.font.italic = 1
      AnnotationAtts.axes2D.yAxis.title.userTitle = 0
      AnnotationAtts.axes2D.yAxis.title.userUnits = 0
      AnnotationAtts.axes2D.yAxis.title.title = "Y-Axis"
      AnnotationAtts.axes2D.yAxis.title.units = ""
      AnnotationAtts.axes2D.yAxis.label.visible = 1
      AnnotationAtts.axes2D.yAxis.label.font.font = AnnotationAtts.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
      AnnotationAtts.axes2D.yAxis.label.font.scale = 1
      AnnotationAtts.axes2D.yAxis.label.font.useForegroundColor = 1
      AnnotationAtts.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes2D.yAxis.label.font.bold = 1
      AnnotationAtts.axes2D.yAxis.label.font.italic = 1
      AnnotationAtts.axes2D.yAxis.label.scaling = 0
      AnnotationAtts.axes2D.yAxis.tickMarks.visible = 1
      AnnotationAtts.axes2D.yAxis.tickMarks.majorMinimum = 0
      AnnotationAtts.axes2D.yAxis.tickMarks.majorMaximum = 1
      AnnotationAtts.axes2D.yAxis.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axes2D.yAxis.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axes2D.yAxis.grid = 0
      AnnotationAtts.axes3D.visible = 0
      AnnotationAtts.axes3D.autoSetTicks = 1
      AnnotationAtts.axes3D.autoSetScaling = 1
      AnnotationAtts.axes3D.lineWidth = 0
      AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
      AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.ClosestTriad  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
      AnnotationAtts.axes3D.triadFlag = 1
      AnnotationAtts.axes3D.bboxFlag = 0
      AnnotationAtts.axes3D.xAxis.title.visible = 1
      AnnotationAtts.axes3D.xAxis.title.font.font = AnnotationAtts.axes3D.xAxis.title.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.xAxis.title.font.scale = 1
      AnnotationAtts.axes3D.xAxis.title.font.useForegroundColor = 1
      AnnotationAtts.axes3D.xAxis.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.xAxis.title.font.bold = 0
      AnnotationAtts.axes3D.xAxis.title.font.italic = 0
      AnnotationAtts.axes3D.xAxis.title.userTitle = 0
      AnnotationAtts.axes3D.xAxis.title.userUnits = 0
      AnnotationAtts.axes3D.xAxis.title.title = "X-Axis"
      AnnotationAtts.axes3D.xAxis.title.units = ""
      AnnotationAtts.axes3D.xAxis.label.visible = 1
      AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.xAxis.label.font.scale = 1
      AnnotationAtts.axes3D.xAxis.label.font.useForegroundColor = 1
      AnnotationAtts.axes3D.xAxis.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.xAxis.label.font.bold = 0
      AnnotationAtts.axes3D.xAxis.label.font.italic = 0
      AnnotationAtts.axes3D.xAxis.label.scaling = 0
      AnnotationAtts.axes3D.xAxis.tickMarks.visible = 1
      AnnotationAtts.axes3D.xAxis.tickMarks.majorMinimum = 0
      AnnotationAtts.axes3D.xAxis.tickMarks.majorMaximum = 1
      AnnotationAtts.axes3D.xAxis.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axes3D.xAxis.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axes3D.xAxis.grid = 0
      AnnotationAtts.axes3D.yAxis.title.visible = 1
      AnnotationAtts.axes3D.yAxis.title.font.font = AnnotationAtts.axes3D.yAxis.title.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.yAxis.title.font.scale = 1
      AnnotationAtts.axes3D.yAxis.title.font.useForegroundColor = 1
      AnnotationAtts.axes3D.yAxis.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.yAxis.title.font.bold = 0
      AnnotationAtts.axes3D.yAxis.title.font.italic = 0
      AnnotationAtts.axes3D.yAxis.title.userTitle = 0
      AnnotationAtts.axes3D.yAxis.title.userUnits = 0
      AnnotationAtts.axes3D.yAxis.title.title = "Y-Axis"
      AnnotationAtts.axes3D.yAxis.title.units = ""
      AnnotationAtts.axes3D.yAxis.label.visible = 1
      AnnotationAtts.axes3D.yAxis.label.font.font = AnnotationAtts.axes3D.yAxis.label.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.yAxis.label.font.scale = 1
      AnnotationAtts.axes3D.yAxis.label.font.useForegroundColor = 1
      AnnotationAtts.axes3D.yAxis.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.yAxis.label.font.bold = 0
      AnnotationAtts.axes3D.yAxis.label.font.italic = 0
      AnnotationAtts.axes3D.yAxis.label.scaling = 0
      AnnotationAtts.axes3D.yAxis.tickMarks.visible = 1
      AnnotationAtts.axes3D.yAxis.tickMarks.majorMinimum = 0
      AnnotationAtts.axes3D.yAxis.tickMarks.majorMaximum = 1
      AnnotationAtts.axes3D.yAxis.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axes3D.yAxis.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axes3D.yAxis.grid = 0
      AnnotationAtts.axes3D.zAxis.title.visible = 1
      AnnotationAtts.axes3D.zAxis.title.font.font = AnnotationAtts.axes3D.zAxis.title.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.zAxis.title.font.scale = 1
      AnnotationAtts.axes3D.zAxis.title.font.useForegroundColor = 1
      AnnotationAtts.axes3D.zAxis.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.zAxis.title.font.bold = 0
      AnnotationAtts.axes3D.zAxis.title.font.italic = 0
      AnnotationAtts.axes3D.zAxis.title.userTitle = 0
      AnnotationAtts.axes3D.zAxis.title.userUnits = 0
      AnnotationAtts.axes3D.zAxis.title.title = "Z-Axis"
      AnnotationAtts.axes3D.zAxis.title.units = ""
      AnnotationAtts.axes3D.zAxis.label.visible = 1
      AnnotationAtts.axes3D.zAxis.label.font.font = AnnotationAtts.axes3D.zAxis.label.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axes3D.zAxis.label.font.scale = 1
      AnnotationAtts.axes3D.zAxis.label.font.useForegroundColor = 1
      AnnotationAtts.axes3D.zAxis.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axes3D.zAxis.label.font.bold = 0
      AnnotationAtts.axes3D.zAxis.label.font.italic = 0
      AnnotationAtts.axes3D.zAxis.label.scaling = 0
      AnnotationAtts.axes3D.zAxis.tickMarks.visible = 1
      AnnotationAtts.axes3D.zAxis.tickMarks.majorMinimum = 0
      AnnotationAtts.axes3D.zAxis.tickMarks.majorMaximum = 1
      AnnotationAtts.axes3D.zAxis.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axes3D.zAxis.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axes3D.zAxis.grid = 0
      AnnotationAtts.axes3D.setBBoxLocation = 0
      AnnotationAtts.axes3D.bboxLocation = (0, 1, 0, 1, 0, 1)
      AnnotationAtts.axes3D.triadColor = (0, 0, 0)
      AnnotationAtts.axes3D.triadLineWidth = 0
      AnnotationAtts.axes3D.triadFont = 0
      AnnotationAtts.axes3D.triadBold = 1
      AnnotationAtts.axes3D.triadItalic = 1
      AnnotationAtts.axes3D.triadSetManually = 0
      AnnotationAtts.userInfoFlag = 0
      AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
      AnnotationAtts.userInfoFont.scale = 1
      AnnotationAtts.userInfoFont.useForegroundColor = 1
      AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
      AnnotationAtts.userInfoFont.bold = 0
      AnnotationAtts.userInfoFont.italic = 0
      AnnotationAtts.databaseInfoFlag = 0
      AnnotationAtts.timeInfoFlag = 1
      AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
      AnnotationAtts.databaseInfoFont.scale = 1
      AnnotationAtts.databaseInfoFont.useForegroundColor = 1
      AnnotationAtts.databaseInfoFont.color = (0, 0, 0, 255)
      AnnotationAtts.databaseInfoFont.bold = 0
      AnnotationAtts.databaseInfoFont.italic = 0
      AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
      AnnotationAtts.databaseInfoTimeScale = 1
      AnnotationAtts.databaseInfoTimeOffset = 0
      AnnotationAtts.legendInfoFlag = 0
      AnnotationAtts.backgroundColor = (255, 255, 255, 255)
      AnnotationAtts.foregroundColor = (0, 0, 0, 255)
      AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
      AnnotationAtts.gradientColor1 = (0, 0, 255, 255)
      AnnotationAtts.gradientColor2 = (0, 0, 0, 255)
      AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
      AnnotationAtts.backgroundImage = ""
      AnnotationAtts.imageRepeatX = 1
      AnnotationAtts.imageRepeatY = 1
      AnnotationAtts.axesArray.visible = 1
      AnnotationAtts.axesArray.ticksVisible = 1
      AnnotationAtts.axesArray.autoSetTicks = 1
      AnnotationAtts.axesArray.autoSetScaling = 1
      AnnotationAtts.axesArray.lineWidth = 0
      AnnotationAtts.axesArray.axes.title.visible = 1
      AnnotationAtts.axesArray.axes.title.font.font = AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axesArray.axes.title.font.scale = 1
      AnnotationAtts.axesArray.axes.title.font.useForegroundColor = 1
      AnnotationAtts.axesArray.axes.title.font.color = (0, 0, 0, 255)
      AnnotationAtts.axesArray.axes.title.font.bold = 0
      AnnotationAtts.axesArray.axes.title.font.italic = 0
      AnnotationAtts.axesArray.axes.title.userTitle = 0
      AnnotationAtts.axesArray.axes.title.userUnits = 0
      AnnotationAtts.axesArray.axes.title.title = ""
      AnnotationAtts.axesArray.axes.title.units = ""
      AnnotationAtts.axesArray.axes.label.visible = 1
      AnnotationAtts.axesArray.axes.label.font.font = AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
      AnnotationAtts.axesArray.axes.label.font.scale = 1
      AnnotationAtts.axesArray.axes.label.font.useForegroundColor = 1
      AnnotationAtts.axesArray.axes.label.font.color = (0, 0, 0, 255)
      AnnotationAtts.axesArray.axes.label.font.bold = 0
      AnnotationAtts.axesArray.axes.label.font.italic = 0
      AnnotationAtts.axesArray.axes.label.scaling = 0
      AnnotationAtts.axesArray.axes.tickMarks.visible = 1
      AnnotationAtts.axesArray.axes.tickMarks.majorMinimum = 0
      AnnotationAtts.axesArray.axes.tickMarks.majorMaximum = 1
      AnnotationAtts.axesArray.axes.tickMarks.minorSpacing = 0.02
      AnnotationAtts.axesArray.axes.tickMarks.majorSpacing = 0.2
      AnnotationAtts.axesArray.axes.grid = 0
      SetAnnotationAttributes(AnnotationAtts)
    # set random color for coordinates lines
    if (opt_random_mesh_color == 1):
      MeshAtts = MeshAttributes()
      MeshAtts.legendFlag = 1
      MeshAtts.lineWidth = 0
      MeshAtts.meshColor = (0, 0, 0, 255)
      MeshAtts.meshColorSource = MeshAtts.MeshRandom  # Foreground, MeshCustom, MeshRandom
      MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom, OpaqueRandom
      MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
      MeshAtts.pointSize = 0.05
      MeshAtts.opaqueColor = (255, 255, 255, 255)
      MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
      MeshAtts.pointSizeVarEnabled = 0
      MeshAtts.pointSizeVar = "default"
      MeshAtts.pointType = MeshAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
      MeshAtts.showInternal = 0
      MeshAtts.pointSizePixels = 2
      MeshAtts.opacity = 1
      SetPlotOptions(MeshAtts)
      MeshAtts = MeshAttributes()
      MeshAtts.legendFlag = 1
      MeshAtts.lineWidth = 0
      MeshAtts.meshColor = (0, 0, 0, 255)
      MeshAtts.meshColorSource = MeshAtts.MeshRandom  # Foreground, MeshCustom, MeshRandom
      MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom, OpaqueRandom
      MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
      MeshAtts.pointSize = 0.05
      MeshAtts.opaqueColor = (255, 255, 255, 255)
      MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
      MeshAtts.pointSizeVarEnabled = 0
      MeshAtts.pointSizeVar = "default"
      MeshAtts.pointType = MeshAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
      MeshAtts.showInternal = 0
      MeshAtts.pointSizePixels = 2
      MeshAtts.opacity = 1
      SetPlotOptions(MeshAtts)
  # saving the session
  SaveSession("{}/visit_plot_mesh{}.session".format(data_path,os.getpid()))

  
##} plot

##{ collect_files to be loaded in visit
def collect_files(data_path,regex):
  allfiles   = os.listdir(data_path)
  files_name = set()
  
  # get all files that match the regex,
  # note files which refer to a different time(cycle) grouped together.
  for f in allfiles:
    if re.search(r'{}\.{}'.format(regex,FILE_SUFFIX),f):
      stem = re.sub(r'\d+\.({})'.format(FILE_SUFFIX),'*.\\1',f)
      files_name.add(stem)
  
  if len(files_name) == 0:
    raise Exception("No data file yielded.")
    
  return files_name
##} collect_files

##{ pars_arguments
def pars_arguments():
  notes = """
  ===================================================================
  A VisIt script to visualize patch, scalar, and vector in Elliptica.
  ===================================================================
  
  1. Usage patch example:
  =======================
  visit -cli -s visit_plot.py -f "(left_central_box.*|left_NS_Up)_X.Y.Z.*"
  
  2. Usage scalar example:
  ========================
  visit -cli -s visit_plot.py -s "enthalpy" -f "left_NS_.*_X.Y.Z.*"
  
  3. Usage patch example:
  =======================
  visit -cli -s visit_plot.py -v "beta_U" -f "left_NS_.*_X.Y.Z.*"
  
  4. Usage no window example:
  ===========================
  visit -nowin -cli -s visit_plot.py -v "beta_U" -f "left_NS_.*_X.Y.Z.*"
          """
  
  parser = argparse.ArgumentParser(description=notes,
                                   formatter_class=RawTextHelpFormatter)
  
  parser.add_argument('-d', action = 'store',default = os.getcwd(), 
                      dest="directory_path", type=str, 
                      help = 'Default is .')
  
  parser.add_argument('-f', action = 'store',default = '', 
                      dest="file_names", type=str, 
                      help = 'Regex format file names.')
  
  parser.add_argument('-s', action = 'store',default = '', 
                      dest="scalar_name", type=str, 
                      help = 'The name of scalar field.')
  
  parser.add_argument('-v', action = 'store',default = '', 
                      dest="vector_name", type=str,
                      help = 'The anem of vector field.')
  
  args = parser.parse_args()
  
  return (str(args.directory_path), str(args.file_names),
          str(args.scalar_name),    str(args.vector_name))
  
##} pars_arguments





if __name__ == '__main__' : main()


