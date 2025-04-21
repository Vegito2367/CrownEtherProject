# Main entry point of program. For now i will make it terminal based input but the eventual goal is to have GUI.
# Enter location of dataset
# Run Analysis - Recognize crown ethers. Lets try isolating crown ethers at first.
# Analysis idea: Open file --> create atoms for each symbol in the file --> start with arbitrary root atom and initiate radius BFS
#
#
#
#
#
#
import os
from index import Main

print(
    "Main entry point of program. For now i will make it terminal based input but the eventual goal is to have GUI."
)
location = input("Enter location of dataset relative to ChemResearch folder: ")
allFiles = os.listdir(location)
mainModule = Main(allFiles, location)
mainModule.parseCrownEthers()
print("Analysis complete. Check CEResults folder for results.")
