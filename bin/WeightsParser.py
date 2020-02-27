#!/usr/bin/env python3

"""
WeightsParser.py
This script is a part of piTargetClassifier project.
This script contains the functions to get the weights, either from default values, or from a file.
Yu Sun, ysun43@ur.rochester.edu, 2018-11
"""

#Default values
WeightsDict_Default_Match={'AT': '1', 'TA': '1', 'CG': '1', 'GC': '1', 'GG': '1'}
WeightsDict_Default_Weight={'AT': '2', 'TA': '2', 'CG': '3', 'GC': '3', 'GG': '1'}

#This function appends the key and its corresponding weight if the key hasn't occurred in the dictionary before.
def AddWeightsToDict(Dict,CuurrKey,CurrWeight,Warn):
    if CuurrKey not in Dict:
        Dict[CuurrKey] = CurrWeight
    else:
        if Warn==True:
            print("Please be aware that the weight file has duplications for pairing: ",CuurrKey,". Use the former ones.")

#This function parses the arguments and gets the WeightsDict for further analysis
def WeightFileParser(UseDefault_Match,UseDefault_Weight,UseDefault_UserDefined,WeightFile):
    Dict={}
    if (UseDefault_Match==True) and (UseDefault_Weight==False) and (UseDefault_UserDefined==False):
        Dict = WeightsDict_Default_Match.copy()
    elif (UseDefault_Match==False) and (UseDefault_Weight==True) and (UseDefault_UserDefined==False):
        Dict = WeightsDict_Default_Weight.copy()
    elif (UseDefault_Match==False) and (UseDefault_Weight==False) and (UseDefault_UserDefined==True):
        Weights = open(WeightFile, "r")
        for line in Weights:
            Pairing = line.strip().split()[0]
            PairingWeight = line.strip().split()[1]
            AddWeightsToDict(Dict, Pairing, PairingWeight, True)
            if Pairing[0] != Pairing[1]:
                PairingRev = Pairing[1] + Pairing[0]
                AddWeightsToDict(Dict, PairingRev, PairingWeight, False)
        Weights.close()
    else:
        print("The input parameters are ambiguous...")
    return Dict
