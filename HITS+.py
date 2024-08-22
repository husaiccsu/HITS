import argparse
import warnings
import math
import numpy as np
from math import factorial
from datetime import datetime
from collections import OrderedDict
import random
import networkx as nx
PPMatrix=[]
PGMatrix=[]
PDMatrix=[]
PCMatrix=[]
GGMatrix=[]
Hete_ProbMatrix=[]
GoList=[]
Domainlist=[]
Complexlist=[] 
proteinlist=[] 
NList=[]
N2List=[] 
CPList=[]
PCList=[]
PGONum=[]
WeightC=[]
WeightD=[]
WeightP=[] 

α=0.1
β=0.2
ε=0.0001
γ=0.2
GOType='P'
PPIFile='BioGRIDPPI.txt'


def combination(m, n):
    return math.factorial(m) // (math.factorial(n) * math.factorial(m - n))

def HDP_Complex(N,M,n,k):
    return combination(M,k)*combination(N-M,n-k)/combination(N,n)

 
def getscore_complex(N,M,n,k):
    scores=0
    counter = 1
    while counter<=k:
          scores=scores+HDP_Complex(N,M,n,counter)
          counter=counter+1
    return scores

 
def Weight_Complex():
   ComplexScores=[0 for _ in range(len(Complexlist))]
   N=len(proteinlist)
   M=len(GoList)
   MaxScore=0
   MinScore=1  
   for i in range(0,len(Complexlist)):
      L=len(CPList[i])    
      r=0
      for k in range(0,len(CPList[i])):
        if (k>=0)and (PGONum[int(CPList[i][k])]>0):
          r=r+1
      a=getscore_complex(N,M,L,r)
      ComplexScores[i]=a
    

   global WeightC    
   WeightC = [[0 for _ in range(N)] for _ in range(N)]
   for i in range(0,N-1):
     for j in range(i+1,N):
       score_a=0
       score_b=0
       score_common=0
       intersection =list(set(PCList[i]) & set(PCList[j]))
       for k in range(0,len(PCList[i])):
          score_a=score_a+ComplexScores[int(PCList[i][k])]
       for k in range(0,len(PCList[j])):
          score_b=score_b+ComplexScores[int(PCList[j][k])]
       for k in range(0,len(intersection)):
          score_common=score_common+ComplexScores[int(intersection[k])]           
       
       if (score_a*score_b*score_common>0):
         value=(score_common*score_common)/(score_a*score_b)
         WeightC[i][j]=value
         WeightC[j][i]=value
         if WeightC[i][j]>MaxScore:
             MaxScore=WeightC[i][j]
         if WeightC[i][j]<MinScore:
             MinScore=WeightC[i][j]  
   for i in range(0,N-1):
     for j in range(i+1,N):
       WeightC[i][j]=(WeightC[i][j]-MinScore)/(MaxScore-MinScore)
       WeightC[j][i]=WeightC[i][j]
    
     
def Weight_Protein():
   listlen=len(proteinlist)
   global WeightP
   WeightP = [[0 for _ in range(listlen)] for _ in range(listlen)]
   NeighborList=[[0 for j in range(0)] for i in range(listlen)]
   PDList=[[0 for j in range(0)] for i in range(listlen)]
   for i in range(0,listlen):  
     NeighborList[i]=list(set(NList[i]) & set(N2List[i]))
     for j in range(0,len(Domainlist)): 
       if PDMatrix[i][j]==1:
          PDList[i].append(j)
          continue
       for k in range(0,len(NeighborList[i])):
          KPos=NeighborList[i][k]
          if PDMatrix[KPos][j]==1:
            PDList[i].append(j) 
            break
      
   MaxScore=0
   MinScore=1      
   for i in range(0,listlen-1):
     print(i+1,'/',len(proteinlist)-1,end='\r')
     for j in range(i+1,listlen): 
       intersectionset=list(set(NeighborList[i]) & set(NeighborList[j]))
       union = list(set(NeighborList[i]) | set(NeighborList[j]))
       scores=0
       if (len(NeighborList[i])>0) and (len(NeighborList[j])>0):
          scores=len(intersectionset)*1.0/math.sqrt((len(NeighborList[i]))*(len(NeighborList[j])))     
       WeightP[i][j]=scores
       WeightP[j][i]=scores            
       if scores>MaxScore:
          MaxScore=scores
       if scores<MinScore:
          MinScore=scores   

def Weight_Domain():
  listlen=len(proteinlist)
  global WeightD
  WeightD = [[0 for _ in range(listlen)] for _ in range(listlen)]
  N_CPList=[[0 for j in range(0)] for i in range(listlen)]
  PDList=[[0 for j in range(0)] for i in range(listlen)]
  NDList=[[0 for j in range(0)] for i in range(listlen)]
  CDList=[[0 for j in range(0)] for i in range(listlen)]
  for i in range(0,listlen): 
      N_CPList[i].append(i)
      if len(PCList[i])==0:
        for k in range(len(NList[i])):
          N_CPList[i].append(NList[i][k])
      else:
        for j in range(0,len(PCList[i])):
          IPos=int(PCList[i][j])
          for k in range(0,len(CPList[IPos])):
            if CPList[IPos][k] not in N_CPList[i]:
              N_CPList[i].append(CPList[IPos][k])
      
  for i in range(0,listlen): 
      for j in range(0,len(Domainlist)):
       if PDMatrix[i][j]==1:
          PDList[i].append(j)
       for k in range(len(NList[i])):
          if PDMatrix[int(NList[i][k])][j]==1:
            NDList[i].append(j) 
            break
       for k in range(len(N_CPList[i])):
          if PDMatrix[int(N_CPList[i][k])][j]==1:
            CDList[i].append(j) 
            break   

            
  M=len(Domainlist) 
  MaxScore=0
  MinScore=1  
  for i in range(0,listlen-1):
    #print(i+1,'/',len(proteinlist),end='\r')
    if PGONum[i]==0:
      continue
    for j in range(i+1,listlen): 
      if PGONum[j]==0:
       continue 
      a1=len(PDList[i])
      b1=len(PDList[j])
      intersection1 =list(set(PDList[i]) & set(PDList[j]))
      Common1=len(intersection1)
      value1=0
      if  (a1>0) and (b1>0):
        value1=(Common1*Common1)/(a1*b1)

      scores=value1           
      WeightD[i][j]=scores
      WeightD[j][i]=scores
      if scores>MaxScore:
         MaxScore=scores
      if scores<MinScore:
         MinScore=scores        
  for i in range(0,listlen-1):
   for j in range(i+1,listlen):
       WeightD[i][j]=(WeightD[i][j]-MinScore)/(MaxScore-MinScore)
       WeightD[j][i]=WeightD[i][j] 
    
def LoadPPI(Files):
  TempList=[]
  with open('SGD.txt', 'r') as file:
    for line in file:
      line = line.strip()
      beginstr,endstr,types=line.split('\t')
      if types!=GOType:
         continue      
      if beginstr not in TempList:
         TempList.append(beginstr)
         
  with open(Files, 'r') as file:
    for line in file:
      line = line.strip()
      beginstr,endstr=line.split('\t')
      if (beginstr not in TempList) or (endstr not in TempList):
         continue
      if beginstr not in proteinlist:
         proteinlist.append(beginstr)
      if endstr not in proteinlist:
         proteinlist.append(endstr)        
    
    listlen=len(proteinlist)
    file.seek(0)
    global PPMatrix
    PPMatrix = [[0 for _ in range(listlen)] for _ in range(listlen)]
    global NList
    NList=[[0 for j in range(0)] for i in range(listlen)]
    global N2List
    N2List=[[0 for j in range(0)] for i in range(listlen)]
    for line in file:
        line = line.strip()
        beginstr,endstr=line.split('\t')
        if (beginstr not in TempList) or (endstr not in TempList):
          continue
        Ipos=proteinlist.index(beginstr)
        JPos=proteinlist.index(endstr)
        PPMatrix[Ipos][JPos]=1
        PPMatrix[JPos][Ipos]=1  
        NList[Ipos].append(JPos)
        NList[JPos].append(Ipos)
        
    for i in range(0,listlen):
      for j in range(0,len(NList[i])):
         IPos=int(NList[i][j])
         for k in range(0,len(NList[IPos])):
            IPos2=int(NList[IPos][k])
            if (IPos2!=i) and (IPos2 not in N2List[i]):
              N2List[i].append(IPos2)

  with open('SGD.txt', 'r') as file:
    for line in file:
      line = line.strip()
      beginstr,endstr,types=line.split('\t')
      if types!=GOType:
         continue
      if beginstr not in proteinlist:
         continue 
      global GoList
      if endstr not in GoList:
         GoList.append(endstr)         

    global PGMatrix
    PGMatrix=[[0 for j in range(len(GoList))] for i in range(listlen)]
    global GPList
    GPList=[[0 for j in range(0)] for i in range(len(GoList))]
    global PGONum
    PGONum=[0 for i in range(listlen)] 
    file.seek(0)
    for line in file:
      line = line.strip()
      beginstr,endstr,types=line.split('\t')
      try:
          Ipos=proteinlist.index(beginstr)
      except ValueError:
          Ipos=-1
      if types!=GOType:
         continue    
      if Ipos==-1:
         continue     
      JPos=GoList.index(endstr)
      PGMatrix[Ipos][JPos]=1
      GPList[JPos].append(Ipos) 
      PGONum[Ipos]=PGONum[Ipos]+1 
           
      
def LoadMultiData(Files,Types):
   with open(Files, 'r') as file:
    for line in file:
        line = line.strip()
        beginstr,endstr=line.split('\t')
        if Types=='D':
          global Domainlist          
          if endstr not in Domainlist:
            Domainlist.append(endstr)
         
        if Types=='C':
           global Complexlist 
           if endstr not in Complexlist:
            Complexlist.append(endstr) 
 
    PSize=len(proteinlist)
    if Types=='D':
      listlen=len(Domainlist)
      global PDMatrix
      PDMatrix=[[0 for j in range(listlen)] for i in range(PSize)]
    else:
      listlen=len(Complexlist)
      global PCMatrix
      global CPList
      global PCList
      PCMatrix=[[0 for j in range(listlen)] for i in range(PSize)]
      CPList=[[0 for j in range(0)] for i in range(listlen)]
      PCList=[[0 for j in range(0)] for i in range(PSize)]
    file.seek(0)
    for line in file:
        line = line.strip()
        beginstr,endstr=line.split('\t')
        try:
          Ipos=proteinlist.index(beginstr)
        except ValueError:
          Ipos=-1

        if Ipos>=0:
          if Types=='D': 
            JPos=Domainlist.index(endstr)
            PDMatrix[Ipos][JPos]=1
          else:
            JPos=Complexlist.index(endstr)
            PCMatrix[Ipos][JPos]=1
            CPList[JPos].append(Ipos) 
            PCList[Ipos].append(JPos) 
      

def ConstructHeteNet():
    G=nx.DiGraph()
    global GGMatrix
    LP=len(proteinlist)
    LG=len(GoList)
    GGMatrix = [[0 for _ in range(len(GoList))] for _ in range(len(GoList))]
    with open('GO_Part.txt', 'r') as file:  
     nums=0
     for line in file:
      line = line.strip()
      beginstr,endstr,types=line.split('\t')
      if types!=GOType:
         continue   
      if (beginstr not in GoList) or (endstr not in GoList):
         continue
      IPos=GoList.index(beginstr)   
      JPos=GoList.index(endstr)
      G.add_edge(IPos,JPos,weight=1)
    
    for i in range(0,LG):    
      for j in range(0,LG):
        if i in G.nodes() and j in G.nodes():
          try:
            sp=nx.shortest_path_length(G, i, j)
          except  nx.NetworkXNoPath:
            sp=0
          if sp>0:
             GGMatrix[i][j]=1.0/pow(2,sp)
          else:
             GGMatrix[i][j]=0  
            
    global Hete_ProbMatrix

    HeteMatrix=[[0 for j in range(LP+LG)] for i in range(LP+LG)]  
    Hete_ProbMatrix=[[0 for j in range(LP+LG)] for i in range(LP+LG)] 
    for i in range(len(proteinlist)):
       HeteMatrix[i][0:LP]=WeightP[i][0:LP]
    for i in range(len(proteinlist)): 
       HeteMatrix[i][LP:LP+LG]=PGMatrix[i][0:LG]   
    
    for i in range(LG):
       HeteMatrix[LP+i][LP:LP+i]=GGMatrix[i][0:LG]  
   
    transposed_matrix = [[row[i] for row in PGMatrix] for i in range(len(PGMatrix[0]))]
    for i in range(LP,LP+LG): 
       HeteMatrix[i][0:LP]=transposed_matrix[i-LP][0:LP]
   


    for i in range(LP):
       SumP=0
       for j in range(LP): 
         SumP+=HeteMatrix[i][j]
       SumG=0
       for j in range(LG): 
         SumG+=PGMatrix[i][j]
       
       for j in range(LP):  
         if SumP==0:
           Hete_ProbMatrix[i][j]=0
         else: 
              Hete_ProbMatrix[i][j]=(1-β)*HeteMatrix[i][j]/SumP
       for j in range(LG):  
          if SumG==0:
            Hete_ProbMatrix[i][j+LP]=0
          else:
            Hete_ProbMatrix[i][j+LP]=β*HeteMatrix[i][j+LP]/SumG
    
    for j in range(LG): 
       SumP=0
       for i in range(LP): 
         SumP+=PGMatrix[i][j]
       
       SumG=0
       for i in range(LG): 
         SumG+=HeteMatrix[LP+i][LP+j]
         
       for i in range(LG): 
         if SumG==0:
           Hete_ProbMatrix[LP+i][LP+j]=0
         else: 
           Hete_ProbMatrix[LP+i][LP+j]=(1-β)*HeteMatrix[LP+i][LP+j]/SumG 

       
       for i in range(LP):                  
          if SumP==0:
            Hete_ProbMatrix[i+LG][j]=0
          else:
            Hete_ProbMatrix[i+LG][j]=β*HeteMatrix[i+LG][j]/SumP           

def ProductMatrixVector(InputMatrix,InputVector,Index):
    OutputVector=np.zeros(len(proteinlist)+len(GoList))
    for i in range(len(GoList)):
      InputMatrix[Index][len(proteinlist)+i]=0
      InputMatrix[len(proteinlist)+i][Index]=0
    OutputVector=np.dot(InputMatrix,InputVector)
    return OutputVector


def HITS():

    filename='HITS_'+PPIFile
    wfile = open(filename, "w") 
    HVector=np.zeros(len(proteinlist)+len(GoList))
    AVector=np.zeros(len(proteinlist)+len(GoList))
    H0Vector=np.zeros(len(proteinlist)+len(GoList))
    A0Vector=np.zeros(len(proteinlist)+len(GoList))
    NewHVector=np.zeros(len(proteinlist)+len(GoList))
    NewAVector=np.zeros(len(proteinlist)+len(GoList))
    TempVector=np.zeros(len(proteinlist)+len(GoList))
    neighborscores=np.zeros(len(proteinlist)+len(GoList))
    HeteM=np.asarray(Hete_ProbMatrix)
    NewMatrix=np.zeros((len(proteinlist),len(proteinlist)))


    for i in range(len(proteinlist)):
       print(i+1,'/',len(proteinlist),end='\r')
       if PGONum[i]==0:
         continue
       for j in range(0,len(proteinlist)):
           H0Vector[j]=WeightD[i][j]
           A0Vector[j]=WeightC[i][j]
      
       for j in range(len(GoList)):        
          MaxH=0
          MaxA=0
          for k in range(len(GPList[j])):  
            Ipos=GPList[j][k]
            if WeightD[i][Ipos]>MaxH:
              MaxH=WeightD[i][Ipos]
            if WeightC[i][Ipos]>MaxA: 
              MaxA=WeightC[i][Ipos]
   

          H0Vector[len(proteinlist)+j]=MaxH
          A0Vector[len(proteinlist)+j]=MaxA
          
       count=1
       for j in range(len(GoList)+len(proteinlist)):   
          AVector[j]=A0Vector[j]
          HVector[j]=H0Vector[j]
       while count<20:
         TempVector=ProductMatrixVector(HeteM,HVector,i)
         for k in range(len(proteinlist)+len(GoList)):
           NewAVector[k]=(1-α)*A0Vector[k]+α*TempVector[k]
         
         TempVector=ProductMatrixVector(HeteM,NewAVector,i)
         for k in range(len(proteinlist)+len(GoList)):
           NewHVector[k]=(1-α)*H0Vector[k]+α*TempVector[k]
        

         
         suma=0
         sumh=0
         for k in range(len(proteinlist)):
           suma=suma+abs(NewAVector[k]-AVector[k])
           sumh=sumh+abs(NewHVector[k]-HVector[k])
         if (suma+sumh)>ε:
           for k in range(len(proteinlist)+len(GoList)): 
             HVector[k]=NewHVector[k]
             AVector[k]=NewAVector[k]
         else:
           break 
         count=count+1
       maxvalue=0
       maxindex=-1 
       for j in range(0,len(proteinlist)):
          neighborscores[j]=γ*HVector[j]+(1-γ)*AVector[j]
          NewMatrix[i][j]=neighborscores[j]
          if maxvalue<NewMatrix[i][j]:
             maxvalue=NewMatrix[i][j]
             maxindex=j 

    
    for i in range(len(proteinlist)):
       print(i+1,'/',len(proteinlist),end='\r')
       if PGONum[i]==0:
         continue
       maxvalue=0
       maxindex=-1 
       GoNeighbor=[]  
       for j in range(0,len(proteinlist)):
         if (i!=j) and (PGONum[j]>0) and (NewMatrix[i][j]>ε):
           GoNeighbor.append(j)
           if maxvalue<NewMatrix[i][j]:
             maxvalue=NewMatrix[i][j]
             maxindex=j       
   
       Scores=dict()
       Benchbark=[]
       for k in range(0,len(GoList)):
        GoScore=0
        if PGMatrix[i][k]==1:
          Benchbark.append(k)
        for j in range(0,len(GoNeighbor)):
            JPos=int(GoNeighbor[j])
            if PGMatrix[JPos][k]==1:
              GoScore=GoScore+NewMatrix[i][JPos]
        if GoScore>0:
          Scores[k]=GoScore
       sorted_dictionary = sorted(Scores.items(), key=lambda item: item[1], reverse=True) 
                
       Ipos=PGONum[maxindex]
       if len(sorted_dictionary)<Ipos:
         Ipos=len(sorted_dictionary)
       keys_list=[]
       for key, value in sorted_dictionary:
         keys_list.append(key)
       pridictedlist=keys_list[:Ipos]
       s=proteinlist[i]+'\t'
       for k in range(0,len(pridictedlist)):
         s=s+GoList[pridictedlist[k]]+'\t'
       wfile.write(s+'\n') 
       
      

    
def main():
  
    print("Loading Data:"+PPIFile, datetime.now().strftime('%Y-%m-%d %H:%M:%S')) 
    LoadPPI(PPIFile)    

    LoadMultiData('Domain.txt','D')
    LoadMultiData('complexes.txt','C')
    print('Weighting based on PPI network:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    Weight_Protein()
    print('Weighting based on complex:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    Weight_Complex()
    print('Weighting based on domain:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    Weight_Domain()

    print('Constructing heterogeneous network:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    ConstructHeteNet()
    print('Running HITS:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    HITS()  


    print('Over', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
   
    
    
if __name__ == "__main__":
	main()
    