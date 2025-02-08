from cifFileParser import Molecule
import os
import numpy as np
from heatmap import InteractivePlot
import pandas as pd
from collections import Counter
class Main:
  def __init__(self, allfiles: list,location:str):
    self.files=allfiles
    self.folder=location
  
  def magnitude(vector):
    mag=0
    for i in vector:
      mag+=(i**2)
    
    return mag**0.5

  def parseCrownEthers(self):
    #BFS
    for file in self.files:
      molecule=Molecule(self.folder+"/"+file)
      if( not molecule.validFile):
        continue
      structure=molecule.getStructure()
      outfile = open(f"CEResults/{file[0:6]}-Structure.txt","w")
      outString=""
      for key in structure:
        outString+=f"{key} - {structure[key]}\n"
      outfile.write(outString)
      outfile.close()

      cycleResults=molecule.detectCrownEthers()
      outfile = open(f"CEResults/{file[0:6]}-Crowns.txt","w")
      outString=""
      for key in cycleResults:
        outString+=f"{key} - {cycleResults[key]}\n"
      outfile.write(outString)
      outfile.close()
      



    



