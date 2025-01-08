import numpy as np
import matplotlib.pyplot as plt
from atom import Atom
import os
from collections import deque
'''
Need all possible covalent Radii to account for all new atoms, need to put the graph into this class because I need to access the getAtoms
in a radius function.
'''
class Molecule:
  @staticmethod
  def printNicely(myList):
    for j in myList:
      print(j)

  

  def __init__(self,filePath):
    self.covalentRadii={#Did not put nitrogen covalent radius 
    "C":[0.75,0.67,0.60],
    "O":[0.63,0.57,0.53]
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

    self.makeGraph()

  def makeGraph(self):
    self.structure={}
    root=None
    for j in self.Atoms:
      if(j.symbol == "O"):
        root=j
        self.structure[j]=[]
        break
    #First lets make structure, then think of bfs to recognize cycles.
    bfsQ=deque()
    bfsQ.append(root)
    visited=set()
    parent=None
    while bfsQ:
      current=bfsQ.popleft()
      visited.add(current)
      surround=self.getAtomsInARadius(current,2)
      self.structure[current]=[]
      for atom,dist in surround:
        if atom not in visited:
          if(current.isBounded(atom,dist)):
            self.addBond(current,atom)
            bfsQ.append(atom)
        else:
          pass


  def getStructure(self):
    return self.structure


    

    
  

  
  def addBond(self,existingAtom,newAtom):
    self.structure[existingAtom].append(newAtom)
    if(newAtom not in self.structure):
      self.structure[newAtom]=[]
    self.structure[newAtom].append(existingAtom)
  
  def returnAtoms(self,symbol):
    output=[]
    keyList=list(self.structure.keys())
    for key in keyList:
      if(key.symbol==symbol):
        output.append(key)
    return output
  
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
  


  
