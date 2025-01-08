from cifFileParser import CIFParser
import os
import numpy as np
from graph import Graph
from render import Render,ExportUnit
from heatmap import InteractivePlot
import pandas as pd
from collections import Counter

def magnitude(vector):
  mag=0
  for i in vector:
    mag+=(i**2)
  
  return mag**0.5

def getAngle(left,center,right):
  v1=[]
  v2=[]
  for i in range(3):
    v1.append(left[i]-center[i])
    v2.append(right[i]-center[i])

  v1=np.array(v1)
  v2=np.array(v2)
  dot=np.dot(v1,v2)
  cosTheta=dot/(magnitude(v1) * magnitude(v2))
  angle=np.arccos(cosTheta)
  return round(np.degrees(angle),6)

def getNormalVector(x,y,z):
  pass
def getTorsionAngle(A,B,C,D):
  #BC is the common edge
  #ABC is first plane and BCD is the second plane
  normABC=np.cross(A.positionVector-B.positionVector,C.positionVector-B.positionVector)
  normBCD=np.cross(B.positionVector-C.positionVector,D.positionVector-C.positionVector)
  return np.degrees(np.arccos(np.dot(normABC,normBCD)/(magnitude(normABC)*magnitude(normBCD))))


######################################################SNSManipulation

def IdentifySNSAngles(parser,lowerLimit,upperLimit,invalidFiles,distanceValues,ExportData):
      sulphurs=[]
      nitrogens=[]

      nitrogens=parser.getElementAtoms("N")
      sulphurs=parser.getElementAtoms("S")
      
      
      for n in nitrogens:
        for s in sulphurs:
          distance=n.getDistance(s)

          if(distance>=lowerLimit and distance<=upperLimit):
            distanceValues[(n,s)]=distance
      
      if(len(distanceValues)==0):
        invalidFiles.append([parser.fileName,"No S-N bonds found"])
        
      # print(f"S-N-S bonds for {file}")
      occurences={}
      for (i,j) in distanceValues:
        if(i in occurences.keys()):
          occurences[i]+=1
        else:
          occurences[i]=1
        
        if(j in occurences.keys()):
          occurences[j]+=1
        else:
          occurences[j]=1
        
      SNSBonds=[]
      for key in occurences.keys():#Cycles over each nitrogen atom in the occurence list
        left,right,center=None,None,None
        if(occurences[key]>=2 and key.symbol=="N"):
            center=key
            bonds=list(distanceValues.keys())
            for n,s in bonds:
              if(n==center):
                if(left is None):
                  left=s
                elif(right is None):
                  right=s

            angle = getAngle(left.positionVector,center.positionVector,right.positionVector)
            g=Graph([left,center,right],angle)
            #Separate the Export stuff from the SNS Bonds Identification
            leftd=distanceValues[(g.center,g.left)]
            rightd=distanceValues[(g.center,g.right)]
            # Extremes(maxProps,minProps,parser.fileName,g.center.symbol,g.left.symbol,leftd)
            # Extremes(maxProps,minProps,parser.fileName,g.center.symbol,g.left.symbol,rightd)
            ExportData.append(ExportUnit(parser.fileName,g.bondAngle,[center,left],leftd,[center,right],rightd))
            SNSBonds.append(g)
      return SNSBonds



def Extremes(maxProps, minProps, file, n, s, distance):
    if(distance>maxProps[0]):
      maxProps[0]=distance
      maxProps[1]=file
      maxProps[2]=f"{n} -- {s}"
    if(distance<minProps[0]):
      minProps[0]=distance
      minProps[1]=file
      minProps[2]=f"{n} -- {s}"
######################################################End of SNSManipulation/Start of Torison Angle


def main():
  pass
  
main()


