from cifFileParser import Molecule


class Main:
    # Main class to parse the CIF files and detect crown ethers

    def __init__(self, allfiles: list, location: str):
        self.files = allfiles
        self.folder = location

    def magnitude(self, vector: list[float]):
        mag = 0
        for i in vector:
            mag += i**2

        return mag**0.5

    def parseCrownEthers(self):
        # BFS
        for file in self.files:
            molecule = Molecule(self.folder + "/" + file)
            if not molecule.validFile:
                continue
            # structure=molecule.getStructure()
            # outfile = open(f"CEResults/{file[0:6]}-Structure.txt","w")
            # outString=""
            # for key in structure:
            #   outString+=f"{key} - {structure[key]}\n"
            # outfile.write(outString)
            # outfile.close()

            oxygenStructure = molecule.detectCrownEthers()
            # cycleResults=molecule.detectCrownEthers()
            outfile = open(f"CEResults/{file[0:6]}-CrownLoop.txt", "w")
            outString = ""
            for key in oxygenStructure:
                outString += f"{key} - {oxygenStructure[key]}\n"
            outfile.write(outString)
            outfile.close()
