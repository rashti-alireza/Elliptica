# A VisIt script to visualize patch/mesh in Elliptica:

# How to invoke the script:
#
# 1. with windows:
# $ visit -cli -s visit_plot_mesh.py "regex name of file name" --path .
# ex: 
# $ visit -cli -s visit_plot_mesh.py \
#              "(left_central_box_left|left_NS_Up)_X.Y.Z.*"
# 
# 2. without windows: 
# $ visit -nowin -cli -s visit_plot_mesh.py "regex name of file name" --path .


# import packages:
from __future__ import division
import re
import os
import argparse

##{ main
def main():
  # getting field name and type
  regex,data_path = pars_arguments()
  
  # deleting all plots
  DeleteAllPlots()
  
  plot(regex, data_path)
  
##} main

##{ plot: given object, data and the section of grid plot the object
def plot(regex,data_path):
  
  # collecting all of the files name needed to load
  files_name = collect_files(regex,data_path)
  
  n = len(files_name)
  
  print('Plotting: mesh')
  for i in range(n):
    print('Openning data file: ' + files_name[i])
  
  SetWindowLayout(1)
  
  # half of the patches
  SetActiveWindow(1)
  for i in range(0,n):
    OpenDatabase(files_name[i], 0)
    mesh_name = re.sub(r'_xyz_.*\.silo','',files_name[i])
    AddPlot('Mesh',mesh_name, 1, 1)
    DrawPlots()
  
##} plot

##{ collect_files to be loaded in visit
def collect_files(regex,data_path):
  allfiles   = os.listdir(data_path)
  files_name = []
  
  for f in allfiles:
    if re.search(r'{}'.format(regex),f):
      files_name.append(f)
  
  if len(files_name) == 0:
    raise Exception("No data file yielded.")
  
  return files_name
##} collect_files

##{ pars_arguments
def pars_arguments():
  notes = '"A VisIt script to visualize patch/mesh in Elliptica"'
  parser = argparse.ArgumentParser(description=notes)
  parser.add_argument('regex'   , action = 'store',type=str, help = 'Python format regex for file names')
  parser.add_argument('--path =', action = 'store',default = os.getcwd(),dest="data_path", type=str, help = 'Default is .')
  args = parser.parse_args()
  
  return str(args.regex), str(args.data_path)
  
##} pars_arguments





if __name__ == '__main__' : main()


