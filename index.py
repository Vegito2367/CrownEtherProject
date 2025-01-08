from cifFileParser import Molecule
import os
import numpy as np
from render import Render,ExportUnit
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
    #BFS/DFS
    for file in self.files:
      molecule=Molecule(self.folder+"/"+file)
      print(molecule.structure)



    



