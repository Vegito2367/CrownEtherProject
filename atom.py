import numpy as np

from numpy.typing import NDArray


class Atom:
    def __init__(
        self,
        inputList: list[str],
        conversionMatrix: NDArray[np.float64],
        radii: dict[str, float],
    ):

        self.covalentRadius = 0.0
        elemname = inputList.pop(0)
        self.identifier = elemname
        ind1, ind2 = 0, 0
        for i in range(len(elemname)):
            if str.isdecimal(elemname[i]):
                ind1 = i
                break
        for i in range(len(elemname) - 1, -1, -1):
            if str.isdecimal(elemname[i]):
                ind2 = i
                break

        if ind1 == 0 and ind2 == 0:
            self.atomLetter = ""
            self.atomNum = 1
            self.element = elemname
        else:
            self.atomLetter = elemname[ind2 + 1 :]
            self.element = elemname[:ind1]
            self.atomNum = elemname[ind1 : ind2 + 1]

        self.symbol = inputList.pop(0)
        vectorList: list[str] = [inputList.pop(0), inputList.pop(0), inputList.pop(0)]

        position_list: list[float] = []

        for v in range(3):
            elem = vectorList[v]
            if "(" in elem:
                position_list.append(float(elem[: elem.index("(")]))
            else:
                position_list.append(float(elem))

        position_array: NDArray[np.float64] = np.array(position_list)
        self.positionVector = np.matmul(conversionMatrix, position_array)
        self.remainingNumbers = inputList

        if self.symbol in radii:
            self.covalentRadius: float = radii[self.symbol]

    def __str__(self):
        return self.identifier

    def __repr__(self):
        return self.identifier

    def getDistance(self, other: "Atom") -> float:
        distanceVector: list[float] = [
            round(self.positionVector[i] - other.positionVector[i], 6) for i in range(3)
        ]
        output = 0
        for j in distanceVector:
            output += j**2

        return round(output**0.5, 3)

    def __eq__(self, other: object) -> bool:
        if other is None:
            return False
        if not isinstance(other, Atom):
            return False
        for i in range(3):
            if self.positionVector[i] != other.positionVector[i]:
                return False

        return True

    def __hash__(self):
        return hash(tuple(self.positionVector))

    def isBounded(self, other: "Atom", dist: float = 0, fudge: float = 0.5):
        try:
            if dist == 0:
                dist = self.getDistance(other)

            out = dist <= (self.covalentRadius + other.covalentRadius + fudge)
            return out
        except TypeError:
            print(self, other)
            raise TypeError


# class BondType(Enum):
#     SINGLE = 1
#     DOUBLE = 2
#     TRIPLE = 3


# class Bond:
#     def __init__(self, type):
#         self.bondType = type
