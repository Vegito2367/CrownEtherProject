import numpy as np
import matplotlib.pyplot as plt
from atom import Atom
import os
'''
Need all possible covalent Radii to account for all new atoms, need to put the graph into this class because I need to access the getAtoms
in a radius function.
'''
class CIFParser:
  @staticmethod
  def printNicely(myList):
    output=""
    for j in myList:
      print(j)

  

  def __init__(self,filePath):
    self.covalentRadii={#Did not put nitrogen covalent radius  
  "Cu":1.32,
  "Ag":1.45,
  "Au":1.36,
  "Pt":1.36,
  "Pd":1.39,
  "Hg":1.32,
  "Fe":1.52,
  "Ru":1.46
}
    
    myfile=open(filePath,'r')
    self.fileName=myfile.name
    self.validFile=True
    alllines=myfile.read().split("\n")
    try:
      startIndex=alllines.index("_atom_site_fract_z")+1
    except:
      self.validFile=False
      return
    self.cellvalues={}
    for i in range(startIndex): #Retrieveing cell properties
      if("_cell_length_a" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_length_a"]=(float(val))

      if("_cell_length_b" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_length_b"]=(float(val))

      if("_cell_length_c" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_length_c"]=(float(val))

      if("_cell_angle_alpha" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_angle_alpha"]=(float(val))

      if("_cell_angle_beta" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_angle_beta"]=(float(val))

      if("_cell_angle_gamma" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_angle_gamma"]=(float(val))
    
    textofInterest=alllines[startIndex:alllines.index("#END")]

    
    astarnum=self.alphaStarNumerator()
    astardenom=self.alphaStarDenominator()

    self.cellvalues["cell_astar"]=np.arccos(astarnum/astardenom)

    ConversionMatrix=[
      [self.a(),self.b() * self.cos(self.gamma()),self.c() * self.cos(self.beta())],
                      
      [0, self.b() * self.sin(self.gamma()), -1 * self.c() * self.sin(self.beta()) * self.cos(self.cellvalues["cell_astar"])],

      [0,0,self.c() * self.sin(self.beta()) * self.sin(self.cellvalues["cell_astar"])]
      
    ]
    self.ConversionMatrix=np.array(ConversionMatrix)

    self.Atoms=[]
    for j in textofInterest:
      self.Atoms.append(Atom(j.split(" "),self.ConversionMatrix,self.covalentRadii))

  def getElementAtoms(self,symbol):
    output=[]
    for j in range(len(self.Atoms)):
      if(self.Atoms[j].symbol == symbol):
        output.append(self.Atoms[j])

    return output
  
  def containsAtom(self,symbol):
    atoms=self.getElementAtoms(symbol)
    return len(atoms)>0
  def getAtomsInARadius(self,targetAtom,radius):
    output=[]
    for atom in self.Atoms:
      dist=atom.getDistance(targetAtom)
      if(dist<=radius):
        output.append([atom,dist])
    return output
  
  def getParticularAtom(self,identifier):
    for atom in self.Atoms:
      if(atom.identifier==identifier):
        return atom
    return "Atom Not Found"
  
  def cos(self,x):
    return np.cos(x)
  
  def sin(self,x):
    return np.sin(x)
  
  
  def alphaStarNumerator(self):
    return (self.cos(self.beta()) * self.cos(self.gamma())) - self.cos(self.alpha())
  
  def alphaStarDenominator(self):
    return self.sin(self.beta()) * self.sin(self.gamma())
  

  def b(self):
    return self.cellvalues["cell_length_b"]
  
  def a(self):
    return self.cellvalues["cell_length_a"]
  
  def c(self):
    return self.cellvalues["cell_length_c"]
  
  def gamma(self):
    return np.radians(self.cellvalues["cell_angle_gamma"])
  
  def beta(self):
    return np.radians(self.cellvalues["cell_angle_beta"])
  
  def alpha(self):
    return np.radians(self.cellvalues["cell_angle_alpha"])
  


  
