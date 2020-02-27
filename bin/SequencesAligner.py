#!/usr/bin/env python3

"""
SequencesAligner.py
This script is a part of piTargetClassifier project.
This script contains the functions to compute piRNA-RNA parining, and returns objects that store matching information.
Requied packages: numpy
Yu Sun, ysun43@ur.rochester.edu, 2018-11
"""

"""
Scoring matrix for Dynamic Programming:
piRNA Length=m
Target Length=N
Matrix dimension: m rows and N+2(m-1) columns
"""

import numpy as np

#This object stores RNA sequence information
class RNAClass:
    def __init__(self, Name, Sequence):
        self.Name = Name
        self.Sequence = Sequence
        self.Length=len(Sequence)

#This object stores pairing results for 1 piRNA and 1 mRNA
class Pairing:
    def __init__(self, TargetName, TargetSeq, TargetLength, piRNAName, piRNASeq, piRNALength, \
                 Matching1_Start, Matching1_End, Matching1_Pattern, Matching1_Score, Top1_BinPos, \
                 Matching2_Start, Matching2_End, Matching2_Pattern, Matching2_Score, Top2_BinPos, \
                 Matching3_Start, Matching3_End, Matching3_Pattern, Matching3_Score, Top3_BinPos):
        self.TargetName = TargetName
        self.TargetSeq = TargetSeq
        self.TargetLength = TargetLength
        self.piRNAName = piRNAName
        self.piRNASeq = piRNASeq
        self.piRNALength = piRNALength

        self.Matching1_Start = Matching1_Start
        self.Matching1_End = Matching1_End
        self.Matching1_Pattern = Matching1_Pattern
        self.Matching1_Score = Matching1_Score
        self.Matching1_BinPos=Top1_BinPos

        self.Matching2_Start = Matching2_Start
        self.Matching2_End = Matching2_End
        self.Matching2_Pattern = Matching2_Pattern
        self.Matching2_Score = Matching2_Score
        self.Matching2_BinPos=Top2_BinPos

        self.Matching3_Start = Matching3_Start
        self.Matching3_End = Matching3_End
        self.Matching3_Pattern = Matching3_Pattern
        self.Matching3_Score = Matching3_Score
        self.Matching3_BinPos=Top3_BinPos

#Get the visualization of alignments from the object Pairing, for the top 1, 2, 3 alignments.
#For partial matching, it will add 'E' (which means empty) in the target sequence.
#In the pattern output file, the pattern is in piRNA 5'->3' direction!
def VisualizePairing(PairingObject,TopID, OutputFileAln, OutputFilePattern):
    fo=open(OutputFileAln,"a")
    fop=open(OutputFilePattern,"a")
    if TopID==1:
        Start=PairingObject.Matching1_Start
        End=PairingObject.Matching1_End
        Pattern=PairingObject.Matching1_Pattern
        Score=PairingObject.Matching1_Score
        BinPos=PairingObject.Matching1_BinPos
    elif TopID==2:
        Start = PairingObject.Matching2_Start
        End = PairingObject.Matching2_End
        Pattern = PairingObject.Matching2_Pattern
        Score = PairingObject.Matching2_Score
        BinPos = PairingObject.Matching2_BinPos
    elif TopID==3:
        Start = PairingObject.Matching3_Start
        End = PairingObject.Matching3_End
        Pattern = PairingObject.Matching3_Pattern
        Score = PairingObject.Matching3_Score
        BinPos = PairingObject.Matching3_BinPos
    else:
        print("TopID out of bound...")
    fo.write("Target RNA: " + PairingObject.TargetName + "; "+ str(Start + 1) + " - " + str(End + 1)+"\n")
    OriSeqTargetSeq=PairingObject.TargetSeq
    EmptyE="E"*(PairingObject.piRNALength-1)
    TargetSeqEmptyE=EmptyE+OriSeqTargetSeq+EmptyE
    fo.write(TargetSeqEmptyE[(Start+(PairingObject.piRNALength-1)):(End + 1 + PairingObject.piRNALength-1)]+"\n")
    MatchingString = ""
    ScoringString = ""
    for pos in range(PairingObject.piRNALength):
        if Pattern[pos] > 0:
            MatchingString = MatchingString + "|"
            ScoringString = ScoringString + str(Pattern[pos])
        else:
            MatchingString = MatchingString + " "
            ScoringString = ScoringString + "0"
    fo.write(MatchingString+"\n")
    fo.write(PairingObject.piRNASeq+"\n")
    fo.write(ScoringString+"\n")
    fop.write(PairingObject.TargetName+"\t"+PairingObject.piRNAName+"\t"+"m"+str(TopID)+"\t"+ScoringString[::-1] + "\t"+ str(BinPos) +"\n")
    fo.write("piRNA: "+ PairingObject.piRNAName + "; " + str(PairingObject.piRNALength) +" - 1 ; m"+str(TopID)+ ", Score: " + str(Score)+"\n")
    fo.write("\n")
    fo.close()
    fop.close()

#This function searches the paired nucleotide (CurrPairing)
# and looks up the weight in the dictionary (WeightsDict) of that pairing.
def GetPairingWeight(WeightsDict, CurrPairing):
    if CurrPairing in WeightsDict:
        w=WeightsDict[CurrPairing]
    else:
        w=0
    return w

#This function takes the Alignment matrix and the end position (at the last row, 0 based) as input,
# and traces back the alignment pattern.
#CurrTrace is the aggragated values, while CurrPattern is the raw matching score list.
def PairingPatternTracer(Matrix,EndPos,piRNALength):
    CurrTrace = [999] * piRNALength
    CurrPattern = [999] * piRNALength
    for i in range(piRNALength):
        CurrTrace[i]=Matrix[i][EndPos-(piRNALength-1)+i]
    for j in range(piRNALength):
        if j==0:
            CurrPattern[j] = CurrTrace[j]
        else:
            CurrPattern[j] = CurrTrace[j]-CurrTrace[j-1]
    return CurrPattern

#CurrPairing: Target nt goes first, and piRNA nt goes second.
#This function generates the scoring matrix. The two inputs are RNAClass objects.
def PairingMatrixGenerator(CurrpiRNA,CurrTarget,Weight):
    AlignmentMatrix = [999] * CurrpiRNA.Length
    AlignmentMatrix_Length = CurrTarget.Length + 2 * (CurrpiRNA.Length - 1)
    for i in range(CurrpiRNA.Length):
        AlignmentMatrix[i] = [999] * AlignmentMatrix_Length
    for j in range(CurrpiRNA.Length):
        for i in range(j, CurrTarget.Length + j + (CurrpiRNA.Length - 1)):
            if (i < CurrpiRNA.Length - 1):
                AlignmentMatrix[j][i] = 0  # Initialize the non-matching parts
            elif j == 0:
                CurrPairing = CurrTarget.Sequence[i - (CurrpiRNA.Length - 1)] + CurrpiRNA.Sequence[j]  #1st nt of piRNA
                CurrWeight = GetPairingWeight(Weight, CurrPairing)
                AlignmentMatrix[j][i] = int(CurrWeight)
                # print(CurrPairing,CurrWeight)
            elif (j > 0) and (i < CurrTarget.Length + CurrpiRNA.Length - 1):
                CurrPairing = CurrTarget.Sequence[i - (CurrpiRNA.Length - 1)] + CurrpiRNA.Sequence[j]
                CurrWeight = GetPairingWeight(Weight, CurrPairing)
                AlignmentMatrix[j][i] = int(AlignmentMatrix[j - 1][i - 1]) + int(CurrWeight)
                # print(CurrPairing)
            else:
                CurrWeight = 0
                AlignmentMatrix[j][i] = int(AlignmentMatrix[j - 1][i - 1]) + int(CurrWeight)
                # print(j,i)
    return AlignmentMatrix

#This function takes the StartPos of the piRNA-RNA targeting event (based on target RNA nt), bins the position on a
# [0, RNALength) scale, and returns a integer
def GetBinnedPosOnRNA(StartPos, RNALength):
    if StartPos<0:
        Frac=0
    elif StartPos == RNALength:
        Frac=RNALength-1
    else:
        Frac=int(StartPos / RNALength * 100)
    if Frac < 100:
        return Frac
    else:
        return (Frac -1)

#This function is a wraper of the pairing result. It takes the alignment matrix and original RNA objects,
# picks the top 3 alignments, and summairzes them into a Pairing object.
def PairingResultsWraper(CurrpiRNA,CurrTarget,AlignmentMatrix):
    ScoringVector = AlignmentMatrix[CurrpiRNA.Length - 1]
    ScoringVector_Cut = ScoringVector[CurrpiRNA.Length - 1:]
    ScoringVector_Cut_Order = np.argsort(ScoringVector_Cut)
    # Top 1, 2, 3 matching. The id may be out-of bound, which indicates partially matching. 0 based pos
    Top1_EndPos = ScoringVector_Cut_Order[-1]
    Top2_EndPos = ScoringVector_Cut_Order[-2]
    Top3_EndPos = ScoringVector_Cut_Order[-3]
    Top1_StartPos = ScoringVector_Cut_Order[-1] - (CurrpiRNA.Length - 1)
    Top2_StartPos = ScoringVector_Cut_Order[-2] - (CurrpiRNA.Length - 1)
    Top3_StartPos = ScoringVector_Cut_Order[-3] - (CurrpiRNA.Length - 1)
    #Get the binned pos on target RNA:
    Top1_BinPos = GetBinnedPosOnRNA(Top1_StartPos, CurrTarget.Length)
    Top2_BinPos = GetBinnedPosOnRNA(Top2_StartPos, CurrTarget.Length)
    Top3_BinPos = GetBinnedPosOnRNA(Top3_StartPos, CurrTarget.Length)
    # Get other information: Top Matching Score
    Top1_Score = ScoringVector_Cut[Top1_EndPos]
    Top2_Score = ScoringVector_Cut[Top2_EndPos]
    Top3_Score = ScoringVector_Cut[Top3_EndPos]
    # Get pattern
    Top1_Pattern = PairingPatternTracer(AlignmentMatrix, Top1_EndPos + (CurrpiRNA.Length - 1), CurrpiRNA.Length)
    Top2_Pattern = PairingPatternTracer(AlignmentMatrix, Top2_EndPos + (CurrpiRNA.Length - 1), CurrpiRNA.Length)
    Top3_Pattern = PairingPatternTracer(AlignmentMatrix, Top3_EndPos + (CurrpiRNA.Length - 1), CurrpiRNA.Length)
    PairingResult = Pairing(CurrTarget.Name, CurrTarget.Sequence, CurrTarget.Length, \
                            CurrpiRNA.Name, CurrpiRNA.Sequence, CurrpiRNA.Length, \
                            Top1_StartPos, Top1_EndPos, Top1_Pattern, Top1_Score, Top1_BinPos, \
                            Top2_StartPos, Top2_EndPos, Top2_Pattern, Top2_Score, Top2_BinPos, \
                            Top3_StartPos, Top3_EndPos, Top3_Pattern, Top3_Score, Top3_BinPos)
    return PairingResult
