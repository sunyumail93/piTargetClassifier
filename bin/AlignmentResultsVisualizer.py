#!/usr/bin/env python3

"""
AlignmentResultsVisualizer.py
This script is a part of piTargetClassifier project.
This script takes the summarized output, and performs paired visualization between piRNA-Target, piRNA-Control.
After running VisualizePairedSummary(), you will get 10 paired plots.
Requied packages: numpy, matplotlib
Yu Sun, ysun43@ur.rochester.edu, 2018-11
"""

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')  #This is required for Linux system to plot figures (MacOS doesn't require)
from bin import MultiSequencesAligner

#Plot paired histograms
def PairedPlotHist(Left, Right, Types, Label, OutputPDFName):
    plt.figure(figsize=(16, 9))
    x1=plt.subplot(1, 2, 1)
    plt.hist(Left)
    plt.xlabel("Scores")
    plt.ylabel("piRNA : Control RNA, " + Label)
    plt.title("piRNA: Control RNA")

    plt.subplot(1, 2, 2, sharey=x1)
    plt.hist(Right)
    plt.xlabel("Scores")
    plt.ylabel("piRNA : Target RNA, " + Label)
    plt.suptitle(Types + ", " + Label)
    plt.title("piRNA: Target RNA")
    plt.savefig(OutputPDFName)

#Plot paired barplots
def PairedPlot(Left, Right, Types, Label, OutputPDFName):
    plt.figure(figsize=(16, 9))
    #plt.title(Types+","+Label)
    x1=plt.subplot(1, 2, 1)
    plt.bar(np.arange(len(Left)), Left)
    plt.xticks(np.arange(len(Left)),range(1, len(Left) + 1))
    plt.xlabel("Positions")
    plt.ylabel("piRNA : Control RNA, "+Label)
    plt.title("piRNA: Control RNA")

    plt.subplot(1, 2, 2, sharey=x1)
    plt.bar(np.arange(len(Right)), Right)
    plt.xticks(np.arange(len(Right)),range(1, len(Right) + 1))
    plt.xlabel("Positions")
    plt.ylabel("piRNA : Target RNA, " + Label)
    plt.title("piRNA: Target RNA")
    plt.suptitle(Types + ", " + Label)
    plt.savefig(OutputPDFName)

def PairedPlotPosHist(Left, Right, Types, Label, OutputPDFName):
    plt.figure(figsize=(16, 9))
    x1=plt.subplot(2, 1, 1)
    plt.hist(Left, bins=100)
    plt.xlabel("Positions on mRNA (binned length 1 - 100)")
    plt.ylabel("Counts, " + Label)
    plt.title("piRNA: Control RNA")

    plt.subplot(2, 1, 2, sharex=x1)
    plt.hist(Right, bins=100)
    plt.xlabel("Positions on mRNA (binned length 1 - 100)")
    plt.ylabel("Counts, " + Label)
    plt.suptitle(Types + ", " + Label)
    plt.title("piRNA: Target RNA")
    plt.savefig(OutputPDFName)

#Wrap all plots into this function
def VisualizePairedSummary(Summary_piRNA_Control, Summary_piRNA_Targets, Prefix):
    print("Generating paired comparison figures...")
    #Plot1: AlignBest Weight Sum Barplot
    PairedPlot(Summary_piRNA_Control.AlignBestPatternsSum[::-1], Summary_piRNA_Targets.AlignBestPatternsSum[::-1], "Best Patterns", \
               "Sum of Weights", Prefix+"PairedPlot.BestPattern.WeightSum.pdf")
    MultiSequencesAligner.ProgressBar(1,12)
    #Plot2: AlignBest Weight Ave Barplot
    PairedPlot(Summary_piRNA_Control.AlignBestPatternsAve[::-1], Summary_piRNA_Targets.AlignBestPatternsAve[::-1], "Best Patterns", \
               "Average of Weights", Prefix+"PairedPlot.BestPattern.WeightAve.pdf")
    MultiSequencesAligner.ProgressBar(2, 12)
    #Plot3: AlignBest Matching Sum Barplot
    PairedPlot(Summary_piRNA_Control.AlignBestMatchingsSum[::-1], Summary_piRNA_Targets.AlignBestMatchingsSum[::-1], "Best Matchings", \
               "Sum of Matchings", Prefix+"PairedPlot.BestMatching.WeightSum.pdf")
    MultiSequencesAligner.ProgressBar(3, 12)
    #Plot4: AlignBest Matching Ave Barplot
    PairedPlot(Summary_piRNA_Control.AlignBestMatchingsAve[::-1], Summary_piRNA_Targets.AlignBestMatchingsAve[::-1], "Best Matchings", \
               "Average of Matchings", Prefix+"PairedPlot.BestMatching.WeightAve.pdf")
    MultiSequencesAligner.ProgressBar(4, 12)
    #Plot5: AlignBest Score Histogram
    PairedPlotHist(Summary_piRNA_Control.AlignBestScores, Summary_piRNA_Targets.AlignBestScores, "Best Patterns", \
                   "Histogram of Scores", Prefix+"PairedPlot.BestPattern.HistoScores.pdf")
    MultiSequencesAligner.ProgressBar(5, 12)
    #Plot6: AlignBest Targed Pos Histogram
    PairedPlotPosHist(Summary_piRNA_Control.AlignBestPosList, Summary_piRNA_Targets.AlignBestPosList, "Best Positions", \
                      "Histogram of Positions", Prefix + "PairedPlot.BestPattern.HistoPos.pdf")
    MultiSequencesAligner.ProgressBar(6, 12)

    #Plot7: AlignAll Weight Sum Barplot
    PairedPlot(Summary_piRNA_Control.AlignAllPatternsSum[::-1], Summary_piRNA_Targets.AlignAllPatternsSum[::-1], "All Patterns", \
               "Sum of Weights", Prefix+"PairedPlot.AllPattern.WeightSum.pdf")
    MultiSequencesAligner.ProgressBar(7, 12)
    #Plot8: AlignAll Weight Ave Barplot
    PairedPlot(Summary_piRNA_Control.AlignAllPatternsAve[::-1], Summary_piRNA_Targets.AlignAllPatternsAve[::-1], "All Patterns", \
               "Average of Weights", Prefix+"PairedPlot.AllPattern.WeightAve.pdf")
    MultiSequencesAligner.ProgressBar(8, 12)
    #Plot9: AlignAll Matching Sum Barplot
    PairedPlot(Summary_piRNA_Control.AlignAllMatchingsSum[::-1], Summary_piRNA_Targets.AlignAllMatchingsSum[::-1], "All Matchings", \
               "Sum of Matchings", Prefix+"PairedPlot.AllMatching.WeightSum.pdf")
    MultiSequencesAligner.ProgressBar(9, 12)
    #Plot10: AlignAll Matching Ave Barplot
    PairedPlot(Summary_piRNA_Control.AlignAllMatchingsAve[::-1], Summary_piRNA_Targets.AlignAllMatchingsAve[::-1], "All Matchings", \
               "Average of Matchings", Prefix+"PairedPlot.AllMatching.WeightAve.pdf")
    MultiSequencesAligner.ProgressBar(10, 12)
    #Plot11: AlignAll Score Histogram
    PairedPlotHist(Summary_piRNA_Control.AlignAllScores, Summary_piRNA_Targets.AlignAllScores, "All Patterns", \
                   "Histogram of Scores", Prefix+"PairedPlot.AllPattern.HistoScores.pdf")
    MultiSequencesAligner.ProgressBar(11, 12)
    #Plot12: AlignAll Score Histogram
    PairedPlotPosHist(Summary_piRNA_Control.AlignBestPosList, Summary_piRNA_Targets.AlignBestPosList, "All Positions", \
                      "Histogram of Positions", Prefix + "PairedPlot.AllPattern.HistoPos.pdf")
    MultiSequencesAligner.ProgressBar(12, 12)
    print("\nDone.")
