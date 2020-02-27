#!/usr/bin/env python3

"""
MachineLearningDataParser.py
This script is a part of piTargetClassifier project.
This script takes the output data of the alignment (*AlignmentPattern.txt) and shapes the data to a numpy array,
 which can be used for the machine learning classifier.
Labels: 0 for piRNA-Control, 1 for piRNA-Target
Requied packages: numpy
Yu Sun, ysun43@ur.rochester.edu, 2018-11
"""

import numpy as np

def DataParserFromPatternOnlyFile(File):
    #Parse file (single column patterns) to np.array
    fi = open(File, "r")
    ML_List = []
    for i in fi:
        data = i.strip()
        datasplit = [int(d) for d in data]
        ML_List.append(datasplit)
    ML_Array = np.array(ML_List)

    return ML_Array

def DataParserFromPatternFile(File, Group):
    if Group == "Control":
        Type=0
    elif Group == "Target":
        Type=1

    #Parse file to np.array
    fi = open(File, "r")
    ML_List = []
    for i in fi:
        data = i.strip().split()[3]
        datasplit = [int(d) for d in data]
        ML_List.append(datasplit)
    ML_Array = np.array(ML_List)

    # Add labels
    datanum = ML_Array.shape[0]
    Labels=np.repeat(Type, datanum)
    return (ML_Array, Labels)

def DataParserFromPatternArray(Array, Group):
    if Group == "Control":
        Type=0
    elif Group == "Target":
        Type=1

    # Parse Data from np.array
    datanum = Array.shape[0]
    Labels = np.repeat(Type, datanum)
    return (Array, Labels)

def DataParserFromPatternFileRf(File, Group):
    if Group == "Control":
        Type="Control"
    elif Group == "Target":
        Type="Target"

    #Parse file to np.array
    fi = open(File, "r")
    ML_List = []
    for i in fi:
        data = i.strip().split()[3]
        datasplit = [int(d) for d in data]
        ML_List.append(datasplit)
    ML_Array = np.array(ML_List)

    # Add labels
    datanum = ML_Array.shape[0]
    Labels=np.repeat(Type, datanum)
    return (ML_Array, Labels)
