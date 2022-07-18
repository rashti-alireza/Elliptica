# A VisIt script to visualize quantities in binary black hole neutron star:

# How to invoke the script:
#
# 1. with windows:
# $ visit -cli -s visit_plot_bbn.py field_name field_type grid_section data_directory
# 
# 2. without windows: 
# $ visit -nowin -cli -s visit_plot_bbn.py field_name field_type grid_section data_directory

# import packages:
from __future__ import division
import re
import os
import argparse

##{ main
def main():
  # getting field name and type
  object_name, object_type, section, data_path = pars_arguments()
  
  # in case only the stem of vector specified, complete the name
  object_name = IsVectorName(object_name,object_type)
  
  # deleting all plots
  DeleteAllPlots()
  
  plot(object_name, object_type, data_path,section)
  
  # NS
  #plot(object_name, object_type, data_path,"NS")
  
  # NS surrounding
  #plot(object_name, object_type, data_path,"NS_surrounding")
  
  # BH surrounding
  #plot(object_name, object_type, data_path,"BH_surrounding")
  
  # outermost0
  #plot(object_name, object_type, data_path,"outermost0")
  
  # outermost1
  #plot(object_name, object_type, data_path,"outermost0")
  
  # for all outermost
  #plot(object_name, object_type, data_path,"outermost")
  
  
##} main

##{ plot: given object, data and the section of grid plot the object
def plot(object_name, object_type, data_path,section):
  
  # collecting all of the files name needed to load
  files_name = collect_files(data_path,section)
  
  n = len(files_name)
  
  print('Plotting: ' + object_name)
  for i in range(n):
    print('Openning data file: ' + files_name[i])
  
  SetWindowLayout(4)
  
  # half of the patches
  SetActiveWindow(1)
  for i in range(0,n//2,1):
    OpenDatabase(files_name[i], 0)
    AddPlot(object_type,object_name, 1, 1)
    DrawPlots()
  
  # half of the patches
  SetActiveWindow(2)
  for i in range(n//2,n,1):
    OpenDatabase(files_name[i], 0)
    AddPlot(object_type,object_name, 1, 1)
    DrawPlots()
  
  # half of the patches
  SetActiveWindow(3)
  for i in range(0,n,2):
    OpenDatabase(files_name[i], 0)
    AddPlot(object_type,object_name, 1, 1)
    DrawPlots()
    
  # all patches
  SetActiveWindow(4)
  for i in range(n):
    OpenDatabase(files_name[i], 0)
    AddPlot(object_type,object_name, 1, 1)
    DrawPlots()
##} plot

##{ collect_files to be loaded in visit
def collect_files(data_path,section):
  allfiles   = os.listdir(data_path)
  files_name = []
  
  if ('NS' == section):
    for f in allfiles:
      if re.search(r'left_NS_(up|down|left|right|back|front)_\w*\.silo$',f) or \
         re.search(r'left_centeral_box_\w*\.silo$',f):
        files_name.append(f)
  
  elif ('NS_surrounding' == section):
    for f in allfiles:
      if re.search(r'left_NS_surrounding_(up|down|left|right|back|front)_\w*\.silo$',f):
        files_name.append(f)
        
  elif ('BH_surrounding' == section):
    for f in allfiles:
      if re.search(r'right_BH_surrounding_(up|down|left|right|back|front)_\w*\.silo$',f):
        files_name.append(f)
        
  elif ('outermost0' == section):
    for f in allfiles:
      if re.search(r'outermost0_(up|down|left|right|back|front)_\w*\.silo$',f):
        files_name.append(f)
        
  elif ('outermost1' == section):
    for f in allfiles:
      if re.search(r'outermost0_(up|down|left|right|back|front)_\w*\.silo$',f):
        files_name.append(f)
  
  elif ('outermost' == section):
    for f in allfiles:
      if re.search(r'outermost[0-9]+_(up|down|left|right|back|front)_\w*\.silo$',f):
        files_name.append(f)
  
  else:
    raise Exception("No such section.")
  
  if len(files_name) == 0:
    raise Exception("No data file yielded.")
  
  return files_name
##} collect_files

##{ IsVectorName, in case only the stem of vector specified, complete the name
def IsVectorName(object_name,object_type):

  if re.search(r'Vector',object_type):
    name = 'Vector_'+object_name+'0_'+object_name+'1_'+object_name+'2'
    object_name = name;
    
  return object_name
##} IsVectorName 

##{ pars_arguments
def pars_arguments():
  notes = '"A VisIt script to visualize quantities in binary black hole neutron star"'
  parser = argparse.ArgumentParser(description=notes)
  parser.add_argument('object_name' , type=str, help = 'The name of the object to be visualized.')
  parser.add_argument('object_type' , type=str, help = 'The type of the object (scalar or vector). '\
                                                       'For vector type only specify the name stem, like Beta_U for vector Beta^{i}.')
  parser.add_argument('grid_section', type=str, help = 'The section of grid you want to plot (NS,NS_surrounding,BH_surrounding,outermost?)')
  parser.add_argument('--path ='   , default = os.getcwd(), action = 'store', dest="data_path", type=str, help = 'Default is the current directory unless otherwise specified.')
  args = parser.parse_args()
  
  object_type = ''
  
  if re.search(r'^(?i)scalar',str(args.object_type)):
    object_type = 'Pseudocolor'
    
  elif re.search(r'^(?i)vector',str(args.object_type)):
    object_type = 'Vector'
    
  #elif re.search(r'^(?i)mesh',str(args.object_type)):
  #  object_type = 'Mesh'
  
  else:
    raise Exception("No correct object type was specified.")
    
  return str(args.object_name), object_type , str(args.grid_section), str(args.data_path)
  
##} pars_arguments



if __name__ == '__main__' : main()


