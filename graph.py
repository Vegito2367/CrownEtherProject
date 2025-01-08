#Each compound will have it own graph, encompassing the entire chem structure
#Each atom will be a node with parent attribute.
class Molecule:
  def __init__(self,root):
    self.structure={}
    self.structure[root]=[]
  
    
  
  
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