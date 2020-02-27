#!/usr/bin/env python3

"""
ProgramUsagePrinter.py
This script is a part of piTargetClassifier project.
This script contains the usage which will be printed out if no inputs passed to this program.
Yu Sun, ysun43@ur.rochester.edu, 2018-11
"""

def PrintUsage():
    print("-h , --help                          Show manual")
    print("Demo mode (--Demo)")
    print("      Demo mode uses the prepared files in ./demo.")
    print("      No other inputs needed. This will take 1 min to run.")

    print("Align mode (--Align)")
    print("           -v        INT             Verbosity, [1, 2, 3]")
    print("      Perform de novo alignments:")
    print("           -p, -pi   FILE            piRNA FASTA file")
    print("           -c, -mc   FILE            Control mRNA FASTA file")
    print("           -t, -mt   FILE            Target mRNA FASTA file")
    print("           -w        Weight/FILE     Weight: 'match' or 'hy'")
    print("                                      or an additional file")
    print("           -o        Prefix          Prefix of output files")
    print("      or Use ready-to-use aligned results:")
    print("           --import  Prefix          Prefix of pre-aligned data")

    print("Learn mode (--Learn)")
    print("           --TestFrac                Test data set fraction")
    print("           -l, --logi                Use logistic regression")
    print("               --CVNum INT           (Optional) Cross Validation N")
    print("           -r, --rf                  Use random forest classifier")
    print("               --TNum  INT           (Optional) Tree Number N")
    print("               --Depth INT           (Optional) Tree Depth N")
    print("           -s, --svm                 Use SVM classifier")
    print("               --Cpen  INT           (Optional) SVM penalty C")
    print("           -a, --all                 Use all above three methods")
    print("           -m, --mode  MODE          (Optional) All/Best/AllBest modes.")

    print("Predict mode (--Predict)")
    print("           --prealign                pre-aligned pattern file")
    print("           --preout                  Prefix of the output file")
    print("")
    print("Demo or Align mode can be run separately.")
    print("Align+Learn, or Align+Learn+Predict modes can be run together.")
    print("For Align mode, please pick one sub-mode: de novo or import.")
    print("    within each sub-mode, all arguments are required")
    print("For Learn mode, one or more modes can be used. Default none.")
    print("For Predict mode, the input contains only a single pattern column.")
    print("    The output file is Predictfile.pre.txt")
    print("Default optional values: CVNum=5, TNum=100, Depth=8, Cpen=1, mode=All")
    print("               TestFrac=0.05")
    print("***********************************************************************")
