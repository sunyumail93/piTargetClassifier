#!/usr/bin/env python3

"""
MultiSequencesAligner.py
This script is a part of piTargetClassifier project.
This script takes the functions from SequencesAligner.py and perform all-to-all alignment
Requied packages: numpy
Yu Sun, ysun43@ur.rochester.edu, 2018-11
"""

import sys
import numpy as np
from bin import SequencesAligner

#This function cleans up the output files (before writing anything in them)
def OutputAlignFilesFlusher(OutputPrefix):
    foalignment_ControlBest = open(OutputPrefix + ".Control.BestAlignment.txt", "w")
    foalignment_ControlAll = open(OutputPrefix + ".Control.AllAlignment.txt", "w")
    fopattern_ControlBest = open(OutputPrefix + ".Control.BestAlignmentPattern.txt", "w")
    fopattern_ControlAll = open(OutputPrefix + ".Control.AllAlignmentPattern.txt", "w")
    foalignment_ControlBest.close()
    foalignment_ControlAll.close()
    fopattern_ControlBest.close()
    fopattern_ControlAll.close()

    foalignment_TargetBest = open(OutputPrefix + ".Targets.BestAlignment.txt", "w")
    foalignment_TargetAll = open(OutputPrefix + ".Targets.AllAlignment.txt", "w")
    fopattern_TargetBest = open(OutputPrefix + ".Targets.BestAlignmentPattern.txt", "w")
    fopattern_TargetAll = open(OutputPrefix + ".Targets.AllAlignmentPattern.txt", "w")
    foalignment_TargetBest.close()
    foalignment_TargetAll.close()
    fopattern_TargetBest.close()
    fopattern_TargetAll.close()

#Use finished percentage to create a progress bar for piRNA-RNA alignment.
def ProgressBar(CurrRNAId, RNAlen):
    Finished = CurrRNAId / RNAlen * 100
    barid = round(Finished / 2 + 0.5)
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %d%%" % ('=' * barid, 2 * barid))
    sys.stdout.flush()

#This function generates a numpy array of objects that stores matching results. Each row is an RNA, and each column is a piRNA.
#piRNA sequences have been reversed in this function.
def ObjectsMatrixGenerator(RNAData, piRNAData, RNAlen, piRNAlen, Weights, OutputFilePrefix, Verbose):
    piRNARNA_AlignObjects = np.zeros((int(RNAlen / 2), int(piRNAlen / 2)), dtype=object)
    Counter=0
    for rnaIndex in range(0, RNAlen, 2):
        Header = RNAData[rnaIndex].strip()[1:]
        Sequence = RNAData[rnaIndex + 1].strip().upper()
        Counter=Counter+1
        if Verbose == True:
            print("Processing RNA: "+Header)
        elif Verbose == False:
            if (Counter==RNAlen/2):
                ProgressBar(Counter, RNAlen/2)
                print("\nDone.")
            else:
                ProgressBar(Counter, RNAlen / 2)

        for piIndex in range(0, piRNAlen, 2):
            piHeader = piRNAData[piIndex].strip()[1:]
            #Updated on 2018-11-26, use reversed (but not complement) sequence!
            piSequence = piRNAData[piIndex + 1].strip().upper()[::-1]
            #print(piHeader + piSequence)

            # Do alignment: piRNA-RNA
            AlignTg_CurrpiRNA = SequencesAligner.RNAClass(piHeader, piSequence)
            AlignTg_CurrTarget = SequencesAligner.RNAClass(Header, Sequence)
            AlignTg_CurrMatrix = SequencesAligner.PairingMatrixGenerator(AlignTg_CurrpiRNA, AlignTg_CurrTarget, Weights)
            AlignTg_CurrResult = SequencesAligner.PairingResultsWraper(AlignTg_CurrpiRNA, AlignTg_CurrTarget, AlignTg_CurrMatrix)
            piRNARNA_AlignObjects[int(rnaIndex / 2), int(piIndex / 2)] = AlignTg_CurrResult

            #Pass the pairing info to the alignment visualizer and save them to output files, for the best and all pairings
            #Alignment.txt contains all the alignments and position, score info; while the Pattern.txt only contains critical
            # info for the classifier: RNA name, piRNA name, topID and scoring pattern
            SequencesAligner.VisualizePairing(AlignTg_CurrResult, 1, OutputFilePrefix+".BestAlignment.txt", OutputFilePrefix+".BestAlignmentPattern.txt")
            SequencesAligner.VisualizePairing(AlignTg_CurrResult, 1, OutputFilePrefix+".AllAlignment.txt", OutputFilePrefix+".AllAlignmentPattern.txt")
            SequencesAligner.VisualizePairing(AlignTg_CurrResult, 2, OutputFilePrefix+".AllAlignment.txt", OutputFilePrefix+".AllAlignmentPattern.txt")
            SequencesAligner.VisualizePairing(AlignTg_CurrResult, 3, OutputFilePrefix+".AllAlignment.txt", OutputFilePrefix+".AllAlignmentPattern.txt")
    return piRNARNA_AlignObjects
