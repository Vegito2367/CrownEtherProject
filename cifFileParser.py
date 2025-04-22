import numpy as np
from atom import Atom
from collections import deque

from numpy.typing import NDArray
from typing import Deque


"""
Need all possible covalent Radii to account for all new atoms, need to put the graph into this class because I need to access the getAtoms
in a radius function.
"""

type Structure = dict[Atom, list[Atom]]  # A dictionary of atoms and their neighbours
type DFSStack = list[Atom]  # A stack of atoms for DFS
type AtomDist = tuple[
    Atom, float
]  # A list of atoms and their distances from the target atom


class Molecule:
    @staticmethod
    def printNicely(myList: list[str]):
        for j in myList:
            print(j)

    def __init__(self, filePath: str):
        self.covalentRadii: dict[str, float] = {  # Did not put nitrogen covalent radius
            "H": 0.32,
            "C": 0.75,
            "N": 0.71,
            "O": 0.63,
            "P": 1.11,
            "K": 1.96,
        }
        self.structure: Structure = {}
        self.oxygenCarbons: Structure = {}

        myfile = open(filePath, "r")
        self.fileName = myfile.name
        self.validFile = True
        alllines = myfile.read().split("\n")
        try:
            startIndex = alllines.index("_atom_site_fract_z") + 1
        except:
            print("Invalid File", filePath)
            self.validFile = False
            return
        self.cellvalues: dict[str, float] = {}
        for i in range(startIndex):  # Retrieveing cell properties
            if "_cell_length_a" in alllines[i]:
                val = alllines[i].split(" ")[1]
                if "(" in val:
                    val = val[: val.index("(")]
                self.cellvalues["cell_length_a"] = float(val)

            if "_cell_length_b" in alllines[i]:
                val = alllines[i].split(" ")[1]
                if "(" in val:
                    val = val[: val.index("(")]
                self.cellvalues["cell_length_b"] = float(val)

            if "_cell_length_c" in alllines[i]:
                val = alllines[i].split(" ")[1]
                if "(" in val:
                    val = val[: val.index("(")]
                self.cellvalues["cell_length_c"] = float(val)

            if "_cell_angle_alpha" in alllines[i]:
                val = alllines[i].split(" ")[1]
                if "(" in val:
                    val = val[: val.index("(")]
                self.cellvalues["cell_angle_alpha"] = float(val)

            if "_cell_angle_beta" in alllines[i]:
                val = alllines[i].split(" ")[1]
                if "(" in val:
                    val = val[: val.index("(")]
                self.cellvalues["cell_angle_beta"] = float(val)

            if "_cell_angle_gamma" in alllines[i]:
                val = alllines[i].split(" ")[1]
                if "(" in val:
                    val = val[: val.index("(")]
                self.cellvalues["cell_angle_gamma"] = float(val)

        textofInterest = alllines[startIndex : alllines.index("#END")]

        astarnum = self.alphaStarNumerator()
        astardenom = self.alphaStarDenominator()

        self.cellvalues["cell_astar"] = np.arccos(astarnum / astardenom)

        ConversionMatrix: list[list[float]] = [
            [
                self.a(),
                self.b() * self.cos(self.gamma()),
                self.c() * self.cos(self.beta()),
            ],
            [
                0,
                self.b() * self.sin(self.gamma()),
                -1
                * self.c()
                * self.sin(self.beta())
                * self.cos(self.cellvalues["cell_astar"]),
            ],
            [
                0,
                0,
                self.c()
                * self.sin(self.beta())
                * self.sin(self.cellvalues["cell_astar"]),
            ],
        ]
        self.ConversionMatrix: NDArray[np.float64] = np.array(ConversionMatrix)

        self.Atoms: list[Atom] = []
        for j in textofInterest:
            self.Atoms.append(
                Atom(j.split(" "), self.ConversionMatrix, self.covalentRadii)
            )

        self.makeGraph()

    def makeGraph(self) -> None:
        """
        Generates two graphs of the molecule using a modified BFS approach
        1. A graph of all the atoms in the molecule
        2. A graph of all the oxygen carbon bonds in the molecule
        """
        searchRadius = 5
        alloxygens: list[Atom] = self.getElementAtoms("O")
        visitedOxygens: set[Atom] = set()

        # First lets make structure, then think of bfs to recognize cycles.
        for oxygen in alloxygens:
            if oxygen in visitedOxygens:
                continue
            bfsQ: Deque[Atom] = deque()
            bfsQ.append(oxygen)
            visited: set[Atom] = set()
            while bfsQ:
                current = bfsQ.popleft()
                if current not in visited:
                    visited.add(current)
                    if current.symbol == "O":
                        visitedOxygens.add(current)
                    surround = self.getAtomsInARadius(current, searchRadius)
                    if current not in self.structure:
                        self.structure[current] = []
                    if (
                        current.symbol == "O" or current.symbol == "C"
                    ) and current not in self.oxygenCarbons:
                        self.oxygenCarbons[current] = []
                    for atom, dist in surround:
                        if atom not in visited:
                            if current.isBounded(atom, dist):
                                if atom not in self.structure:
                                    self.structure[atom] = []
                                if (atom.symbol == "O" or atom.symbol == "C") and (
                                    current.symbol == "O" or current.symbol == "C"
                                ):
                                    if atom not in self.oxygenCarbons:
                                        self.oxygenCarbons[atom] = []
                                    self.oxygenCarbons[current].append(atom)
                                    self.oxygenCarbons[atom].append(current)

                                self.structure[current].append(atom)
                                self.structure[atom].append(current)
                                bfsQ.append(atom)

    def detectCrownEthers(self):
        """
        Attempts to detect crowns with the segmented DFS approach
        Only look for C-O bonds to remove unnecessary clutter
        """
        pass
        cycles: dict[Atom, list[list[Atom]]] = {}
        allOxygens = self.getElementAtoms("O")
        visitedOxygens: set[Atom] = set()
        for oxygen in allOxygens:
            cycles[oxygen] = []
            if oxygen in visitedOxygens:
                continue
            stack: Deque[DFSStack] = deque()
            stack.append([oxygen, oxygen])
            visited: set[Atom] = set()
            prev = None
            while stack:
                current, parent = stack[-1]
                if prev != current:
                    prev = current
                else:
                    break
                if current not in visited:
                    visited.add(current)
                    if current.symbol == "O":
                        visitedOxygens.add(current)
                    adjacent = self.oxygenCarbons[current]
                    for atom in adjacent:
                        if atom == parent:
                            continue
                        stack.append([atom, current])
                else:
                    current, _ = stack.pop()
                    cycle = [current]
                    while stack and stack[-1][0] != current:
                        cycle.append(stack.pop()[0])

                    cycles[oxygen].append(cycle)

        return cycles

    def getOxygenCarbons(self):
        return self.oxygenCarbons

    def getStructure(self):
        return self.structure

    def returnAtoms(self, symbol: str) -> list[Atom]:
        """
        Returns all the atoms of the element 'symbol' in the molecule but by traversing the graph
        """
        output: list[Atom] = []
        keyList = list(self.structure.keys())
        for key in keyList:
            if key.symbol == symbol:
                output.append(key)
        return output

    def getElementAtoms(self, symbol: str) -> list[Atom]:
        """
        Returns all the atoms of the element 'symbol' in the molecule
        """
        output: list[Atom] = []
        for j in range(len(self.Atoms)):
            if self.Atoms[j].symbol == symbol:
                output.append(self.Atoms[j])

        return output

    def containsAtom(self, symbol: str) -> bool:
        """
        Returns True if the molecule contains an atom of the element 'symbol'
        """
        atoms = self.getElementAtoms(symbol)
        return len(atoms) > 0

    def getAtomsInARadius(self, targetAtom: Atom, radius: float) -> list[AtomDist]:
        """
        Returns all the atoms in the molecule that are within a certain radius of the target atom
        """
        output: list[AtomDist] = []
        for atom in self.Atoms:
            dist = atom.getDistance(targetAtom)
            if dist <= radius and dist != 0:
                output.append((atom, dist))

        return output

    def getParticularAtom(self, identifier: str) -> Atom | str:
        """
        Returns the atom with the identifier 'identifier' in the molecule
        """
        for atom in self.Atoms:
            if atom.identifier == identifier:
                return atom
        return "Atom Not Found"

    def cos(self, x: float) -> float:
        return np.cos(x)

    def sin(self, x: float) -> float:
        return np.sin(x)

    def alphaStarNumerator(self):
        return (self.cos(self.beta()) * self.cos(self.gamma())) - self.cos(self.alpha())

    def alphaStarDenominator(self):
        return self.sin(self.beta()) * self.sin(self.gamma())

    def b(self) -> float:
        return self.cellvalues["cell_length_b"]

    def a(self) -> float:
        return self.cellvalues["cell_length_a"]

    def c(self) -> float:
        return self.cellvalues["cell_length_c"]

    def gamma(self) -> float:
        return np.radians(self.cellvalues["cell_angle_gamma"])

    def beta(self) -> float:
        return np.radians(self.cellvalues["cell_angle_beta"])

    def alpha(self) -> float:
        return np.radians(self.cellvalues["cell_angle_alpha"])
