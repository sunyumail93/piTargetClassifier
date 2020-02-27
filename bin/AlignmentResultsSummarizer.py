#!/usr/bin/env python3

"""
AlignmentResultsSummarizer.py
This script is a part of piTargetClassifier project.
This script takes the output from ObjectsMatrixGenerator (an object matrix) and provides a summary.
Requied packages: numpy
Yu Sun, ysun43@ur.rochester.edu, 2018-11
"""

import numpy as np

#This class provides a summary of alignment results.
class AlignmentSummary:
    def __init__(self, AlignBestPatterns, AlignBestPatternsSum, AlignBestPatternsAve, \
                 AlignBestScores, AlignBestScoresFreqDict, AlignBestMatchings, AlignBestMatchingsSum, AlignBestMatchingsAve, AlignBestPosList, \
                 AlignAllPatterns, AlignAllPatternsSum, AlignAllPatternsAve, \
                 AlignAllScores, AlignAllScoresFreqDict, AlignAllMatchings, AlignAllMatchingsSum, AlignAllMatchingsAve, AlignAllPosList):
        self.AlignBestPatterns = AlignBestPatterns
        self.AlignBestPatternsSum = AlignBestPatternsSum
        self.AlignBestPatternsAve=AlignBestPatternsAve
        self.AlignBestScores=AlignBestScores
        self.AlignBestScoresFreqDict=AlignBestScoresFreqDict
        self.AlignBestMatchings=AlignBestMatchings
        self.AlignBestMatchingsSum=AlignBestMatchingsSum
        self.AlignBestMatchingsAve=AlignBestMatchingsAve
        self.AlignBestPosList=AlignBestPosList
        self.AlignAllPatterns=AlignAllPatterns
        self.AlignAllPatternsSum=AlignAllPatternsSum
        self.AlignAllPatternsAve=AlignAllPatternsAve
        self.AlignAllScores=AlignAllScores
        self.AlignAllScoresFreqDict=AlignAllScoresFreqDict
        self.AlignAllMatchings=AlignAllMatchings
        self.AlignAllMatchingsSum=AlignAllMatchingsSum
        self.AlignAllMatchingsAve=AlignAllMatchingsAve
        self.AlignAllPosList=AlignAllPosList

#This function gets the best matching (rather than pattern) from an numpy pattern array.
#It will change all non-zero values to 1, and keeps all zero values.
def GetMatching(Array):
    list = []
    for i in Array:
        rows = []
        for j in i:
            if j != 0:
                z = 1
            else:
                z = 0
            rows.append(z)
        list.append(rows)
    return np.array(list)

#This function provides a summary of piRNA-RNA alignment object matrix, by reporting the sum, ave, scores of all pairings
def SummairzeObjectsMatrix(Matrix):
    AlignBestPatternsList = []
    AlignAllPatternsList = []
    AlignBestPosList=[]
    AlignAllPosList=[]
    for i in range(Matrix.shape[0]):
        for j in range(Matrix.shape[1]):
            CurrObj = Matrix[i, j]
            #print(CurrObj.TargetName, CurrObj.piRNAName)
            #print(CurrObj.Matching1_Pattern, CurrObj.Matching1_Score)
            AlignBestPatternsList.append(CurrObj.Matching1_Pattern)  # Append to list, and finally convert to array
            AlignAllPatternsList.append(CurrObj.Matching1_Pattern)
            AlignAllPatternsList.append(CurrObj.Matching2_Pattern)
            AlignAllPatternsList.append(CurrObj.Matching3_Pattern)
            AlignBestPosList.append(CurrObj.Matching1_BinPos)
            AlignAllPosList.append(CurrObj.Matching1_BinPos)
            AlignAllPosList.append(CurrObj.Matching2_BinPos)
            AlignAllPosList.append(CurrObj.Matching3_BinPos)
    AlignBestPatterns = np.array(AlignBestPatternsList)  # Convert to numpy array
    AlignBestPatternsSum = np.sum(AlignBestPatterns, axis=0).tolist()
    AlignBestPatternsAve = np.mean(AlignBestPatterns, axis=0).tolist()
    AlignBestScores = np.sum(AlignBestPatterns, axis=1).tolist()
    AlignBestScoresFreqDict = {item: AlignBestScores.count(item) for item in set(AlignBestScores)}

    AlignBestMatchings = GetMatching(AlignBestPatterns)
    AlignBestMatchingsSum = np.sum(AlignBestMatchings, axis=0).tolist()
    AlignBestMatchingsAve = np.mean(AlignBestMatchings, axis=0).tolist()

    AlignAllPatterns = np.array(AlignAllPatternsList)  # Convert to numpy array
    AlignAllPatternsSum = np.sum(AlignAllPatternsList, axis=0).tolist()
    AlignAllPatternsAve = np.mean(AlignAllPatternsList, axis=0).tolist()
    AlignAllScores = np.sum(AlignAllPatternsList, axis=1).tolist()
    AlignAllScoresFreqDict = {item: AlignAllScores.count(item) for item in set(AlignAllScores)}

    AlignAllMatchings = GetMatching(AlignAllPatterns)
    AlignAllMatchingsSum = np.sum(AlignAllMatchings, axis=0).tolist()
    AlignAllMatchingsAve = np.mean(AlignAllMatchings, axis=0).tolist()

    Results=AlignmentSummary(AlignBestPatterns,AlignBestPatternsSum,AlignBestPatternsAve,AlignBestScores,AlignBestScoresFreqDict, \
                             AlignBestMatchings, AlignBestMatchingsSum, AlignBestMatchingsAve, AlignBestPosList, \
                             AlignAllPatterns,AlignAllPatternsSum,AlignAllPatternsAve,AlignAllScores,AlignAllScoresFreqDict, \
                             AlignAllMatchings, AlignAllMatchingsSum, AlignAllMatchingsAve, AlignAllPosList)
    return Results
