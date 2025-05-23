# Each compound will have it own graph, encompassing the entire chem structure
# Each atom will be a node with parent attribute.
from atom import Atom


# class Molecule:
#     def __init__(self, root, parser):
#         self.structure = {}
#         self.structure[root] = []
#         self.parser = parser

#     def addBond(self, existingAtom, newAtom):
#         self.structure[existingAtom].append(newAtom)
#         if newAtom not in self.structure:
#             self.structure[newAtom] = []
#         self.structure[newAtom].append(existingAtom)

#     def returnAtoms(self, symbol):
#         output = []
#         keyList = list(self.structure.keys())
#         for key in keyList:
#             if key.symbol == symbol:
#                 output.append(key)
#         return output


# class Node:
#   def __init__(self,atom:Atom,parent:object=None):
#     self.atom=atom
#     if parent is None:
#       self.parent=None
#     elif isinstance(parent,Atom):
#       self.parent=parent

#   def __hash__(self):
#     return hash(self.atom)

#   def __str__(self):
#     pI=""
#     if self.parent is None:
#       pI="None"
#     else:
#       pI=self.parent.atom
#     return f"{self.atom.identifier} - {self.parent.atom}"

#   def __repr__(self):
#     return f"{self.atom.identifier} - {self.parent.atom}"


class CrownEther:
    def __init__(self, crownAtoms:list[Atom]):
        pass
