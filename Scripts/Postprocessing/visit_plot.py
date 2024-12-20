# Alireza Rashti
# January 2021

# A VisIt script to visualize patch, scalar, and vector in Elliptica.
# How to invoke the script:
# invoke:
# $ visit -cli -s visit_plot.py -h



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

## global vars:
glob_cycle_flg = 1

##{ main
def main():

  # get args
  directory_path,file_names,scalar_name,vector_name = pars_arguments()

  # deleting all plots
  DeleteAllPlots()
  
  # plot
  plot(directory_path,file_names,scalar_name,vector_name)
  
  exit(0)
  
##} main

##{ plot: given object, data and the section of grid plot the object
def plot(directory_path,file_names,scalar_name,vector_name):
  
  # collecting all of the files name needed to load
  file_set = set_of_files(directory_path,file_names)
  
  for s in file_set:
    print('Openning data file: ' + s)
  
  SetWindowLayout(1)
  SetActiveWindow(1)
  
  # use to activate or deactivate the relevant options
  scalar_flg = 0
  mesh_flg   = 0
  
  for s in file_set:
    # if this dir has cycle files
    if glob_cycle_flg == 1:
      dbfile = directory_path + '/' + s + ' database'
      OpenDatabase(dbfile, 0)
    else:
      dbfile = directory_path + '/' + s
      OpenDatabase(dbfile, 0)
    
    # scalar
    if (scalar_name != ''):
      AddPlot("Pseudocolor", scalar_name, 1, 1)
      scalar_flg = 1
    # patch
    else:
      mesh_name = re.sub(r'_xyz_.*\.silo','',s)
      AddPlot('Mesh',mesh_name, 1, 1)
      mesh_flg = 1
    
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
    if (mesh_flg == 1 and opt_random_mesh_color == 1):
      MeshAtts = MeshAttributes()
      MeshAtts.legendFlag = 1
      MeshAtts.lineWidth = 0
      MeshAtts.meshColor = (0, 0, 0, 255)
      MeshAtts.meshColorSource = MeshAtts.MeshRandom  # Foreground, MeshCustom, MeshRandom
      MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom, OpaqueRandom
      MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
      MeshAtts.pointSize = 0.05
      MeshAtts.opaqueColor = (255, 255, 255, 255)
      #MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
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
      #MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
      MeshAtts.pointSizeVarEnabled = 0
      MeshAtts.pointSizeVar = "default"
      MeshAtts.pointType = MeshAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
      MeshAtts.showInternal = 0
      MeshAtts.pointSizePixels = 2
      MeshAtts.opacity = 1
      SetPlotOptions(MeshAtts)
  # saving the session
  SaveSession("{}/visit_plot_saved_{}.session"
             .format(directory_path,os.getpid()))

  
##} plot

##{ set_of_files to be loaded in visit
def set_of_files(data_path,regex):
  global glob_cycle_flg
  allfiles  = os.listdir(data_path)
  file_set  = set()
  
  # check if this dir has time(cycle) for each file so set glob_cycle_flg = 1.
  # it makes difference for OpenDatabase command.
  # NOTE: it assumes all files are either with cycle or all are not.
  glob_cycle_flg = 0
  for f1 in allfiles:
    
    if glob_cycle_flg == 1:
      break
    
    # if it's not a sile file
    if not re.search(r'{}\.{}'.format(regex,FILE_SUFFIX),f1):
      continue
    
    for f2 in allfiles:
      # if it's not a sile file
      if not re.search(r'{}\.{}'.format(regex,FILE_SUFFIX),f2):
        continue
      if (f1 == f2):
        continue
      
      # compare stems
      f1_stem = re.sub(r'\d+\.({})'.format(FILE_SUFFIX),'',f1)
      if re.search(r'{}'.format(f1_stem),f2):
        glob_cycle_flg = 1
        break
  
  # get all files that match the regex,
  # note files which refer to a different time(cycle) grouped together.
  for f in allfiles:
    if re.search(r'{}\.{}'.format(regex,FILE_SUFFIX),f):
      if glob_cycle_flg == 1:
        stem = re.sub(r'\d+\.({})'.format(FILE_SUFFIX),'*.\\1',f)
        file_set.add(stem)
      else:
        file_set.add(f)
  
  if len(file_set) == 0:
    raise Exception("No data file yielded.")
    
  return file_set
##} set_of_files

##{ pars_arguments
def pars_arguments():
  notes = """
  ===================================================================
  A VisIt script to visualize patch, scalar, and vector in Elliptica.
  ===================================================================
  
  1. Usage patch example:
  =======================
  visit -cli -s visit_plot.py -F "(left_central_box.*|left_NS_Up)_X.Y.Z.*"
  
  2. Usage scalar example:
  ========================
  visit -cli -s visit_plot.py -S "enthalpy" -F "left_NS_.*_X.Y.Z.*"
  
  3. Usage vector example (not ready yet!):
  ========================================
  visit -cli -s visit_plot.py -V "beta_U" -F "left_NS_.*_X.Y.Z.*"
  
  4. Usage no window example:
  ===========================
  visit -nowin -cli -s visit_plot.py -S "phi" -F "left_NS_.*_X.Y.Z.*"
  
  NOTE: flags are case sensitive!
          """
  
  parser = argparse.ArgumentParser(description=notes,
                                   formatter_class=RawTextHelpFormatter)
  
  parser.add_argument('-D', action = 'store',default = os.getcwd(), 
                      dest="directory_path", type=str, 
                      help = 'Default is .')
  
  parser.add_argument('-F', action = 'store',default = '', 
                      dest="file_names", type=str, 
                      help = 'Regex format file names.')
  
  parser.add_argument('-S', action = 'store',default = '', 
                      dest="scalar_name", type=str, 
                      help = 'The name of scalar field.')
  
  parser.add_argument('-V', action = 'store',default = '', 
                      dest="vector_name", type=str,
                      help = 'The anem of vector field.')
  
  args = parser.parse_args()
  
  # remove / 
  args.directory_path = re.sub('/$','',args.directory_path)
  return (str(args.directory_path), str(args.file_names),
          str(args.scalar_name),    str(args.vector_name))
  
##} pars_arguments





if __name__ == '__main__' : main()


