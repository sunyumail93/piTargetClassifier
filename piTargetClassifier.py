#!/usr/bin/env python3
"""
piTargetClassifier:
A classifier for biological and non-biological relevant piRNA-mRNA pairing, based on an experimental data.
Requied packages: numpy, matplotlib, sklearn
Output:
Prefix.*txt
  In the *Alignment.txt files, piRNAs are reversed, while in *Pattern.txt files, patterns are in piRNA 5´-> 3´ direction.
Prefix.*pdf
  Paired comparison of scoring, patterns and matchings.
Yu Sun, ysun43@ur.rochester.edu, 2018-11
"""

import sys
import time
import datetime
import argparse
import numpy as np
import matplotlib.pyplot as plt
from bin import ProgramUsagePrinter
from bin import SequencesAligner
from bin import WeightsParser
from bin import MultiSequencesAligner
from bin import AlignmentResultsSummarizer
from bin import AlignmentResultsVisualizer
from bin import MachineLearningDataParser
from bin import MachineLearningModels

#Starting time stamp
start = time.time()
if len(sys.argv) == 1:
    print("***********************************************************************")
    print("*                Welcome to use piTargetClassifier                    *")
    print("*                        Y. Sun, 2018-11                              *")
    print("***********************************************************************")
    ProgramUsagePrinter.PrintUsage()
else:
    print("***********************************************************************")
    print("*                Welcome to use piTargetClassifier                    *")
    print("*                        Y. Sun, 2018-11                              *")
    print("***********************************************************************")

    parser = argparse.ArgumentParser()
    parser.add_argument("--Demo", help="Use files in ./demo. No other inputs needed.",action="store_true")
    parser.add_argument("--Align", help="Use Align mode", action="store_true")
    parser.add_argument("-v", type=int, dest="Verbose", help="(Optional) Verbosity, default 1",default=1)
    parser.add_argument("-p","--pi", type=str, dest="PiRNA_File", help="piRNA FASTA file")   #Set an argument with file inputs using dest option.
    parser.add_argument("-c", "--mc", type=str, dest="RNA_Control_File", help="Control mRNA FASTA file")
    parser.add_argument("-t", "--mt", type=str, dest="RNA_Target_File", help="Target mRNA FASTA file")
    parser.add_argument("-w", type=str, dest="Weight", help="Weight: 'match' or 'hy', or a file", default="hy")
    parser.add_argument("-o", type=str, dest="Out_Prefix", help="Prefix of output files")
    parser.add_argument("--import", type=str, dest="Import_Prefix", help="Prefix of pre-aligned data")

    parser.add_argument("--Learn", help="Use Learn mode", action="store_true")
    parser.add_argument("--TestFrac", type=float, dest="TestFrac", help="(Optional) Test data fraction, default 0.05", default=0.05)
    parser.add_argument("-l", "--logi", action="store_true", help="Use logistic regression", default=False)
    parser.add_argument("--CVNum", type=int, dest="CVNum", help="(Optional) Cross Validation N, default 5", default=5)
    parser.add_argument("-r", "--rf", action="store_true", help="Use random forest classifier", default=False)
    parser.add_argument("--TNum", type=int, dest="TNum", help="(Optional) Tree Number N, default 100", default=100)
    parser.add_argument("--Depth", type=int, dest="Depth", help="(Optional) Tree Depth N, default 8", default=8)
    parser.add_argument("-s", "--svm", action="store_true", help="Use SVM classifier", default=False)
    parser.add_argument("--Cpen", type=int, dest="Cpen", help="(Optional) SVM penalty C, default 1", default=1)
    parser.add_argument("-a", "--all", action="store_true", help="Use all above three methods", default=False)
    parser.add_argument("-m", "--mode", type=str, dest="Mode", help="(Optional) All/Best/AllBest modes, default All", default="All")

    parser.add_argument("--Predict", help="Use Predict mode", action="store_true")
    parser.add_argument("--prealign", type=str, dest="Pre_Align", help="pre-aligned patterns (single column).")
    parser.add_argument("--preout", type=str, dest="Pre_Out", help="Prefix of the output file.", default="DefaultPre")

    #args is a list of arguments. If default values not set, it will be None or False
    args = parser.parse_args()

    #Set default values
    PiRNAFile="None"
    RNAControlFile="None"
    RNATargetingFile="None"
    OutputPrefix="None"
    WeightFile="None"
    UseDefault_Match = False
    UseDefault_Weight = False
    UseDefault_UserDefined = False
    V_pi=False
    V_con=False
    V_ML=False
    if args.Verbose == 2:
        V_pi=True
    elif args.Verbose == 3:
        V_pi=True
        V_con=True
        V_ML=True
        print("Using highest verbosity. A full list of parameters:")
        print(args)
    else:
        pass

    #Setting input files
    if (args.Demo == True) or (args.Align == True):
        # Get inputs for Demo or Align (de novo) modes
        if (args.Demo == True) or ((args.Align == True) and (args.Import_Prefix == None)):
            if args.Demo == True:
                print("0. Starting Demo mode using input files from ./demo/")
                print("This Demo data will take ~1 min to finish.")
                PiRNAFile = "./demo/demo.piRNA.fa"
                RNAControlFile = "./demo/demo.Control.fa"
                RNATargetingFile = "./demo/demo.Target.fa"
                UseDefault_Match = False
                UseDefault_Weight = False
                UseDefault_UserDefined = True
                WeightFile = "./demo/weight.hy.txt"
                OutputPrefix = "DemoRun"
                V_pi = False
                V_con = False
                V_ML = False
                args.TestFrac=0.33
                args.all=True
                args.Pre_Align="./demo/demo.UnknownPatterns.txt"
                args.Pre_Out="DemoPre"
                print("This test run uses the following command:")
                print("python3 piTargetClassifier.py \\")
                print("      --Align -v 1 -p ./demo/demo.piRNA.fa ")
                print("          -c ./demo/demo.Control.fa -t ./demo/demo.Target.fa\\")
                print("          -o DemoRun --TestFrac 0.33 -w ./demo/weight.hy.txt \\")
                print("      --Learn -a -m AllBest \\")
                print("      --Predict --prealign ./demo/demo.UnknownPatterns.txt \\")
                print("          --preout DemoPre\n")
            elif (args.Align == True) and (args.Import_Prefix == None):
                PiRNAFile = args.PiRNA_File
                RNAControlFile = args.RNA_Control_File
                RNATargetingFile = args.RNA_Target_File
                if PiRNAFile == None or RNAControlFile == None or RNATargetingFile == None:
                    print("Missing input files...")
                    exit()
                WeightValue=args.Weight
                if WeightValue == "hy":
                    UseDefault_Match = False
                    UseDefault_Weight = True
                    UseDefault_UserDefined = False
                elif WeightValue == "match":
                    UseDefault_Match = False
                    UseDefault_Weight = True
                    UseDefault_UserDefined = False
                else:
                    UseDefault_Match = False
                    UseDefault_Weight = False
                    UseDefault_UserDefined = True
                    WeightFile = args.Weight
                OutputPrefix=args.Out_Prefix

            # Starting Align mode
            print("1. Starting piTargetClassifier Align mode with input FASTA files...")

            # Read inputs and get statistics
            PiRNAFi = open(PiRNAFile, "r")
            RNAControlFi = open(RNAControlFile, "r")
            RNATargetingFi = open(RNATargetingFile, "r")
            PiRNA = PiRNAFi.readlines()
            RNAControl = RNAControlFi.readlines()
            RNATargeting = RNATargetingFi.readlines()
            PiRNA_Len = len(PiRNA)
            RNAControl_Len = len(RNAControl)
            RNATargeting_Len = len(RNATargeting)
            print(">>> Input statistics:")
            print("Control RNAs: " + str(int(RNAControl_Len / 2)))
            print("Target RNAs: " + str(int(RNATargeting_Len / 2)))
            print("piRNAs: " + str(int(PiRNA_Len / 2)))

            # Get Weights
            print(">>> Getting weights and input files")
            WeightsDict = WeightsParser.WeightFileParser(UseDefault_Match, UseDefault_Weight, UseDefault_UserDefined, WeightFile)
            print(">>> Using weight: ")
            print(WeightsDict)

            # Report other parameters:
            print(">>> Other parameters")
            print("Verbosity: "+str(args.Verbose))
            print("Output Prefix: "+str(OutputPrefix))

            # Perform multi-alignment. The returned result is a numpy Array of Objects
            MultiSequencesAligner.OutputAlignFilesFlusher(OutputPrefix)
            # print(V_pi, V_con)
            print(">>> Aligning piRNAs with Target RNAs")
            AlignMatrix_piRNA_Targets = MultiSequencesAligner.ObjectsMatrixGenerator(RNATargeting, PiRNA,
                                                                                     RNATargeting_Len, PiRNA_Len,
                                                                                     WeightsDict,
                                                                                     OutputPrefix + ".Targets", V_pi)
            print(">>> Aligning piRNAs with Control RNAs")
            AlignMatrix_piRNA_Control = MultiSequencesAligner.ObjectsMatrixGenerator(RNAControl, PiRNA, RNAControl_Len,
                                                                                     PiRNA_Len, WeightsDict,
                                                                                     OutputPrefix + ".Control", V_con)

            # Get summary. The result is an AlignmentSummary object, which contains patterns/matchings, sum and average values.
            Summary_piRNA_Targets = AlignmentResultsSummarizer.SummairzeObjectsMatrix(AlignMatrix_piRNA_Targets)
            Summary_piRNA_Control = AlignmentResultsSummarizer.SummairzeObjectsMatrix(AlignMatrix_piRNA_Control)

            AlignmentResultsVisualizer.VisualizePairedSummary(Summary_piRNA_Control, Summary_piRNA_Targets, OutputPrefix)
            PiRNAFi.close()
            RNAControlFi.close()
            RNATargetingFi.close()
            print("***********************************************************************")

        elif (args.Align == True) and (args.Import_Prefix != None):
            # This part won't run into real alignment steps.
            print("1. Starting piTargetClassifier Align mode with input FASTA files...")
            OutputPrefix = args.Import_Prefix
            print(">>> Getting pre-aligned result Prefix: " + str(OutputPrefix))
            print("***********************************************************************")

    else:
        print("At least --Demo or --Align mode needs to be set.")
        print("To check manual, please use -h, --help")
        print("Exitting...")
        print("*                      Programme finished.                         *")
        print("********************************************************************")
        exit()

    if (args.Demo == True) or (args.Learn == True):
        print("\n2. Starting piTargetClassifier Learn mode...")
        print(">>> Learn Mode parameters:")
        print("Test dataset fraction: "+str(args.TestFrac))

        #All
        (MLControl_Data_All, MLControl_Label_All) = MachineLearningDataParser.DataParserFromPatternFile(OutputPrefix + ".Control.AllAlignmentPattern.txt", "Control")
        (MLTarget_Data_All, MLTarget_Label_All) = MachineLearningDataParser.DataParserFromPatternFile(OutputPrefix + ".Targets.AllAlignmentPattern.txt", "Target")
        #Best
        (MLControl_Data_Best, MLControl_Label_Best) = MachineLearningDataParser.DataParserFromPatternFile(OutputPrefix + ".Control.BestAlignmentPattern.txt", "Control")
        (MLTarget_Data_Best, MLTarget_Label_Best) = MachineLearningDataParser.DataParserFromPatternFile(OutputPrefix + ".Targets.BestAlignmentPattern.txt", "Target")

        if (args.all == True) or (args.logi ==True):
            if args.Mode == "All" or args.Mode == "AllBest":
                #All
                print(">>> Logistic: Model Fitting using All Alignments")
                print("Cross Validation N: "+str(args.CVNum))
                OutputSummaryFile_LogiA = OutputPrefix + ".ML.All.Logistic.CV" + str(args.CVNum) + ".txt"
                OutputSummaryFig_LogiA = OutputPrefix + ".ML.All.Logistic.CV" + str(args.CVNum) + ".pdf"
                X_train_logia, X_test_logia, y_train_logia, y_test_logia, logi_a = MachineLearningModels.FitClassifierLogistic(MLControl_Data_All,
                        MLControl_Label_All,MLTarget_Data_All,MLTarget_Label_All,args.CVNum,args.TestFrac,OutputSummaryFile_LogiA)
                MachineLearningModels.CompareMLResultsAndROCCurves(X_train_logia, X_test_logia, y_train_logia, y_test_logia, logi_a, "Logistic",OutputSummaryFile_LogiA, OutputSummaryFig_LogiA, V_ML)
                if args.Mode == "All":
                    print("***********************************************************************")

            if args.Mode == "Best" or args.Mode == "AllBest":
                #Best
                print(">>> Logistic: Model Fitting using Best Alignments")
                print("Cross Validation N: " + str(args.CVNum))
                OutputSummaryFile_LogiB = OutputPrefix + ".ML.Best.Logistic.CV" + str(args.CVNum) + ".txt"
                OutputSummaryFig_LogiB = OutputPrefix + ".ML.Best.Logistic.CV" + str(args.CVNum) + ".pdf"
                X_train_logib, X_test_logib, y_train_logib, y_test_logib, logi_b = MachineLearningModels.FitClassifierLogistic(MLControl_Data_Best,
                        MLControl_Label_Best, MLTarget_Data_Best, MLTarget_Label_Best, args.CVNum, args.TestFrac, OutputSummaryFile_LogiB)
                MachineLearningModels.CompareMLResultsAndROCCurves(X_train_logib, X_test_logib, y_train_logib, y_test_logib,logi_b, "Logistic", OutputSummaryFile_LogiB,OutputSummaryFig_LogiB, V_ML)
                print("***********************************************************************")

        if args.all == True or args.rf ==True:
            if args.Mode == "All" or args.Mode == "AllBest":
                # All
                print(">>> Random Forest: Model Fitting using All Alignments")
                OutputSummaryFile_RFA = OutputPrefix + ".ML.All.RandomForest.T" + str(args.TNum) + "D" + str(args.Depth) + ".txt"
                OutputSummaryFig_RFA = OutputPrefix + ".ML.All.RandomForest.T" + str(args.TNum) + "D" + str(args.Depth) + ".pdf"
                X_train_rfa, X_test_rfa, y_train_rfa, y_test_rfa, rf_a = MachineLearningModels.FitClassifierRandomForest(MLControl_Data_All,
                        MLControl_Label_All, MLTarget_Data_All, MLTarget_Label_All, args.TNum, args.Depth, args.TestFrac,OutputSummaryFile_RFA)
                MachineLearningModels.CompareMLResultsAndROCCurves(X_train_rfa, X_test_rfa, y_train_rfa, y_test_rfa, rf_a,"Random Forest", OutputSummaryFile_RFA,OutputSummaryFig_RFA, V_ML)
                if args.Mode == "All":
                    print("***********************************************************************")

            if args.Mode == "Best" or args.Mode == "AllBest":
                #Best
                print(">>> Random Forest: Model Fitting using Best Alignments")
                OutputSummaryFile_RFB = OutputPrefix + ".ML.Best.RandomForest.T" + str(args.TNum) + "D" + str(args.Depth) + ".txt"
                OutputSummaryFig_RFB = OutputPrefix + ".ML.Best.RandomForest.T" + str(args.TNum) + "D" + str(args.Depth) + ".pdf"
                X_train_rfb, X_test_rfb, y_train_rfb, y_test_rfb, rf_b = MachineLearningModels.FitClassifierRandomForest(MLControl_Data_Best,
                        MLControl_Label_Best, MLTarget_Data_Best, MLTarget_Label_Best, args.TNum, args.Depth, args.TestFrac, OutputSummaryFile_RFB)
                MachineLearningModels.CompareMLResultsAndROCCurves(X_train_rfb, X_test_rfb, y_train_rfb, y_test_rfb, rf_b,"Random Forest", OutputSummaryFile_RFB,OutputSummaryFig_RFB, V_ML)
                print("***********************************************************************")

        if args.all == True or args.svm ==True:
            if args.Mode == "All" or args.Mode == "AllBest":
                #All
                print(">>> SVM: Model Fitting using All Alignments")
                OutputSummaryFile_SVMA = OutputPrefix + ".ML.All.SVM.C" + str(args.Cpen) + ".txt"
                OutputSummaryFig_SVMA = OutputPrefix + ".ML.All.SVM.C" + str(args.Cpen) + ".pdf"
                X_train_svma, X_test_svma, y_train_svma, y_test_svma, svcfig_a = MachineLearningModels.FitClassifierSVM(MLControl_Data_All,
                        MLControl_Label_All,MLTarget_Data_All,MLTarget_Label_All, args.Cpen,args.TestFrac,OutputSummaryFile_SVMA)
                MachineLearningModels.CompareMLResultsAndROCCurves(X_train_svma, X_test_svma, y_train_svma, y_test_svma, svcfig_a, "SVM",OutputSummaryFile_SVMA, OutputSummaryFig_SVMA, V_ML)
                if args.Mode == "All":
                    print("***********************************************************************")

            if args.Mode == "Best" or args.Mode == "AllBest":
                #Best
                print(">>> SVM: Model Fitting using Best Alignments")
                OutputSummaryFile_SVMB = OutputPrefix + ".ML.Best.SVM.C" + str(args.Cpen) + ".txt"
                OutputSummaryFig_SVMB = OutputPrefix + ".ML.Best.SVM.C" + str(args.Cpen) + ".pdf"
                X_train_svmb, X_test_svmb, y_train_svmb, y_test_svmb, svcfig_b = MachineLearningModels.FitClassifierSVM(MLControl_Data_Best,
                        MLControl_Label_Best, MLTarget_Data_Best, MLTarget_Label_Best, args.Cpen, args.TestFrac, OutputSummaryFile_SVMB)
                MachineLearningModels.CompareMLResultsAndROCCurves(X_train_svmb, X_test_svmb, y_train_svmb, y_test_svmb, svcfig_b, "SVM", OutputSummaryFile_SVMB,OutputSummaryFig_SVMB, V_ML)
                print("***********************************************************************")

    if (args.Demo == True) or (args.Predict == True):
        print("3. Starting piTargetClassifier Predict mode...")
        Predict_File=args.Pre_Align
        print("Use the learnt patterns from "+str(args.Mode)+" mode.")
        #All
        print(">>> Reading the input alignment patterns from: "+str(Predict_File))
        Predict_Data = MachineLearningDataParser.DataParserFromPatternOnlyFile(Predict_File)

        if (args.all == True) or (args.logi ==True):
            if args.Mode == "All" or args.Mode == "AllBest":
                #All
                PredictResult_LogiA = args.Pre_Out + ".Predict.ML.All.Logistic.CV" + str(args.CVNum) + ".txt"
                MachineLearningModels.PredictResultsFromModel(Predict_Data, Predict_File, logi_a, "Logistic", PredictResult_LogiA)
            if args.Mode == "Best" or args.Mode == "AllBest":
                #Best
                PredictResult_LogiB = args.Pre_Out + ".Predict.ML.Best.Logistic.CV" + str(args.CVNum) + ".txt"
                MachineLearningModels.PredictResultsFromModel(Predict_Data, Predict_File, logi_b, "Logistic",PredictResult_LogiB)

        if (args.all == True) or (args.rf ==True):
            if args.Mode == "All" or args.Mode == "AllBest":
                #All
                PredictResult_RFA = args.Pre_Out + ".Predict.ML.All.RandomForest.T" + str(args.TNum) + "D" + str(args.Depth) + ".txt"
                MachineLearningModels.PredictResultsFromModel(Predict_Data, Predict_File, rf_a, "Random Forest",PredictResult_RFA)
            if args.Mode == "Best" or args.Mode == "AllBest":
                #Best
                PredictResult_RFB = args.Pre_Out + ".Predict.ML.Best.RandomForest.T" + str(args.TNum) + "D" + str(args.Depth) + ".txt"
                MachineLearningModels.PredictResultsFromModel(Predict_Data, Predict_File, rf_b, "Random Forest",PredictResult_RFB)

        if (args.all == True) or (args.svm ==True):
            if args.Mode == "All" or args.Mode == "AllBest":
                #All
                PredictResult_SVMA = args.Pre_Out + ".Predict.ML.All.SVM.C" + str(args.Cpen) + ".txt"
                MachineLearningModels.PredictResultsFromModel(Predict_Data, Predict_File, svcfig_a, "SVM",PredictResult_SVMA)
            if args.Mode == "Best" or args.Mode == "AllBest":
                #Best
                PredictResult_SVMB = args.Pre_Out + ".Predict.ML.Best.SVM.C" + str(args.Cpen) + ".txt"
                MachineLearningModels.PredictResultsFromModel(Predict_Data, Predict_File, svcfig_b, "SVM",PredictResult_SVMB)
        print("Output file prefix: "+str(args.Pre_Out))

    # Finishing
    end = time.time()
    elapse = datetime.timedelta(seconds=end - start)
    print("Thanks for useing piTargetClassifier. Time elapsed: " + str(elapse))
    print("*                        Programme finished.                          *")
    print("***********************************************************************")
