#!/usr/bin/env python3

"""
MachineLearningModels.py
This script is a part of piTargetClassifier project.
This script takes the training and testing datasets as input, use assigned tree number and depth to perform Random Forest classification.
The outputs will be a summary report and a paired comparison figure for training and testing datasets.
Requied packages: numpy, sklearn, matplotlib.
Yu Sun, ysun43@ur.rochester.edu, 2018-11
"""

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
plt.switch_backend('agg')   #This is required for Linux system to plot figures (MacOS doesn't require)
import numpy as np
from bin import MachineLearningDataParser

#This is a Logistic classifier. The Cross validation number CVNum can be adjusted.
def FitClassifierLogistic(Control_Data, Control_Label, Target_Data, Target_Label, CVNum, test_frac ,OutputSummaryFileName):
    fo = open(OutputSummaryFileName, "w")
    ML_Data = np.concatenate((Control_Data, Target_Data))
    ML_Label = np.append(Control_Label, Target_Label)
    X_train, X_test, y_train, y_test = train_test_split(ML_Data, ML_Label, test_size=test_frac, random_state=0)
    print("Staring Logistic Regressor")
    print("Training data: " + str(X_train.shape[0]))
    print("Testing data: " + str(X_test.shape[0]))
    print("Fitting Logistic Regressor with +"+str(CVNum)+" Cross Validation")
    fo.write("Summary report of the Logistic (LOG) Regressor with "+str(CVNum)+" Cross Validation\n")
    logi = LogisticRegressionCV(cv=CVNum, solver="liblinear", random_state=42)
    logi.fit(X_train, y_train)
    fo.close()
    return X_train, X_test, y_train, y_test, logi

#This is a Random Forest classifier. Tree depth and tree number can be adjusted.
def FitClassifierRandomForest(Control_Data, Control_Label, Target_Data, Target_Label, TreeNum, Depth, test_frac ,OutputSummaryFileName):
    fo = open(OutputSummaryFileName, "w")
    ML_Data = np.concatenate((Control_Data, Target_Data))
    ML_Label = np.append(Control_Label, Target_Label)
    X_train, X_test, y_train, y_test = train_test_split(ML_Data, ML_Label, test_size=test_frac, random_state=0)
    print("Staring Random Forest Classifier")
    print("Training data: " + str(X_train.shape[0]))
    print("Testing data: " + str(X_test.shape[0]))
    print("Fitting Random Forest Classifier with "+str(TreeNum)+" trees and maximum depth "+str(Depth)+".")
    fo.write("Summary report of the Random Forest (RF) Classifier\n")
    fo.write("Random Forest Classifier with "+str(TreeNum)+" trees and maximum depth "+str(Depth)+".\n")
    rf = RandomForestClassifier(n_estimators=TreeNum, max_depth=Depth, random_state=42)
    rf.fit(X_train, y_train)
    fo.close()
    return X_train, X_test, y_train, y_test, rf

#This is a SVM classifier, using rbf kernal and auto gamma. C is the penalty score which is adjustable.
def FitClassifierSVM(Control_Data, Control_Label, Target_Data, Target_Label, C, test_frac ,OutputSummaryFileName):
    fo = open(OutputSummaryFileName, "w")
    ML_Data = np.concatenate((Control_Data, Target_Data))
    ML_Label = np.append(Control_Label, Target_Label)
    X_train, X_test, y_train, y_test = train_test_split(ML_Data, ML_Label, test_size=test_frac, random_state=0)
    print("Staring Support Vector Machine (SVM) Classifier")
    print("Training data: " + str(X_train.shape[0]))
    print("Testing data: " + str(X_test.shape[0]))
    print("Fitting Support Vector Machine (SVM)...This may take a while...")
    fo.write("Summary report of the Support Vector Machine (SVM) Classifier\n")
    svcfit = SVC(kernel='rbf', gamma='auto', C=C, probability=True, random_state=42)  #Default to use Gaussian kernel (rbf), turn on probability so that the ROC curve can be drawn.
    svcfit.fit(X_train, y_train)
    fo.close()
    return X_train, X_test, y_train, y_test, svcfit

#Convert 0/1 numpy array into Predicted Results
def Convert01ToPredictedResults(Num):
    if Num == 1:
        return "Target"
    elif Num == 0:
        return "Control"

#Predict results
def PredictResultsFromModel(Data, InputFile, Classifier, ClassifierType, OutputSummaryFileName):
    fi=open(InputFile, "r")
    fo = open(OutputSummaryFileName, "w")
    print("Running Prediction using "+str(ClassifierType)+"...")
    Predicted=Classifier.predict(Data)
    #print(Predicted)
    count=0
    for line in fi:
        count=count+1
        Pattern=line.strip()
        fo.write(Pattern+"\t"+Convert01ToPredictedResults(Predicted[count-1])+"\n")
    fi.close()
    fo.close()

#Evaluate classifier
def CompareMLResultsAndROCCurves(X_train, X_test, y_train, y_test, Classifier, ClassifierType, OutputSummaryFileName, OutputPDFName, Verbose):
    print("Evaluating the classifier: "+ClassifierType)
    fo=open(OutputSummaryFileName,"a")
    fo.write("Evaluating the classifier: "+ClassifierType+"\n")
    fo.write("Training data: " + str(X_train.shape[0])+"\n"+"Testing data: " + str(X_test.shape[0])+"\n")
    #Analyze training set:
    datasettype = "training"
    data_x = X_train
    data_y = y_train
    y_pred = Classifier.predict(data_x)
    print("Use training datasets:")
    fo.write("\nUse training datasets:\n")
    print("Accuracy of "+ClassifierType+" Classifier on " + datasettype + " dataset: {:.2f}".format(
        Classifier.score(data_x, data_y) * 100) + " %")
    fo.write("Accuracy of "+ClassifierType+" Classifier on " + datasettype + " dataset: {:.2f}".format(
        Classifier.score(data_x, data_y) * 100) + " %\n")

    if Verbose == True:
        print("Confusion matrix:")
    fo.write("Confusion matrix:\n")
    tablewidth = "{0:20}{1:10}{2:10}"
    tn, fp, fn, tp=confusion_matrix(data_y, y_pred).ravel()
    if Verbose == True:
        print(tablewidth.format("Actual\Predicted","0","1"))
        print(tablewidth.format("   Actual 0   ", str(tn), str(fp)))
        print(tablewidth.format("   Actual 1   ", str(fn), str(tp)))
        print("TN="+str(tp)+"\tFP="+str(fp)+"\tFN="+str(fn)+"\tTP="+str(tp))
    fo.write(tablewidth.format("Actual\Predicted","0","1")+"\n")
    fo.write(tablewidth.format("   Actual 0   ", str(tn), str(fp))+"\n")
    fo.write(tablewidth.format("   Actual 1   ", str(fn), str(tp))+"\n")
    fo.write("TN="+str(tn)+"\tFP="+str(fp)+"\tFN="+str(fn)+"\tTP="+str(tp)+"\n")
    if Verbose == True:
        print("Summary report:")
        print(classification_report(data_y, y_pred))
    fo.write("Summary report:\n")
    fo.write(classification_report(data_y, y_pred)+"\n\n")

    plt.figure(figsize=(16, 9))
    x1 = plt.subplot(1, 2, 1)
    curve_roc_auc = roc_auc_score(data_y, Classifier.predict(data_x))
    fpr, tpr, thresholds = roc_curve(data_y, Classifier.predict_proba(data_x)[:, 1])
    plt.plot(fpr, tpr, label=ClassifierType+" Classifier (area = %0.3f)" % curve_roc_auc)
    plt.plot([0, 1], [0, 1], marker='.', lw=2, ls='--', mfc='g', mec='g', color='r')
    plt.xlim([0.0, 1.02])
    plt.ylim([0.0, 1.02])
    plt.xlabel("False Positive Rate (1-Specificity)")
    plt.ylabel("True Positive Rate (Sensitivity)")
    plt.grid(True, linewidth=1, linestyle='--')
    plt.title("Training dataset, accuracy: "+"{:.2f}".format(Classifier.score(data_x, data_y) * 100) + " %")
    plt.legend(loc="lower right")


    # Analyze testing set:
    datasettype = "testing"
    data_x = X_test
    data_y = y_test
    y_pred = Classifier.predict(data_x)
    print("Use testing datasets:")
    fo.write("Use testing datasets:\n")
    print("Accuracy of "+ClassifierType+" Classifier on " + datasettype + " dataset: {:.2f}".format(
        Classifier.score(data_x, data_y) * 100) + " %\n")
    fo.write("Accuracy of "+ClassifierType+" Classifier on " + datasettype + " dataset: {:.2f}".format(
        Classifier.score(data_x, data_y) * 100) + " %\n")
    if Verbose == True:
        print("Confusion matrix:")
    fo.write("Confusion matrix:\n")
    tablewidth = "{0:20}{1:10}{2:10}"
    tn, fp, fn, tp=confusion_matrix(data_y, y_pred).ravel()
    if Verbose == True:
        print(tablewidth.format("Actual\Predicted","0","1"))
        print(tablewidth.format("   Actual 0   ", str(tn), str(fp)))
        print(tablewidth.format("   Actual 1   ", str(fn), str(tp)))
        print("TN="+str(tp)+"\tFP="+str(fp)+"\tFN="+str(fn)+"\tTP="+str(tp))
    fo.write(tablewidth.format("Actual\Predicted","0","1")+"\n")
    fo.write(tablewidth.format("   Actual 0   ", str(tn), str(fp))+"\n")
    fo.write(tablewidth.format("   Actual 1   ", str(fn), str(tp))+"\n")
    fo.write("TN="+str(tn)+"\tFP="+str(fp)+"\tFN="+str(fn)+"\tTP="+str(tp)+"\n")
    if Verbose == True:
        print("Summary report:")
        print(classification_report(data_y, y_pred))
    fo.write("Summary report:\n")
    fo.write(classification_report(data_y, y_pred)+"\n")

    plt.subplot(1, 2, 2, sharey=x1)
    curve_roc_auc = roc_auc_score(data_y, Classifier.predict(data_x))
    fpr, tpr, thresholds = roc_curve(data_y, Classifier.predict_proba(data_x)[:, 1])
    plt.plot(fpr, tpr, label=ClassifierType+" Classifier (area = %0.3f)" % curve_roc_auc)
    plt.plot([0, 1], [0, 1], marker='.', lw=2, ls='--', mfc='g', mec='g', color='r')
    plt.xlim([0.0, 1.02])
    plt.ylim([0.0, 1.02])
    plt.xlabel("False Positive Rate (1-Specificity)")
    plt.ylabel("True Positive Rate (Sensitivity)")
    plt.grid(True, linewidth=1, linestyle='--')
    plt.title("Testing dataset, accuracy: " + "{:.2f}".format(Classifier.score(data_x, data_y) * 100) + " %")
    plt.legend(loc="lower right")

    plt.suptitle("Comparison of Receiver operating characteristic (ROC) curves on training and testing datasets: "+ClassifierType)
    plt.savefig(OutputPDFName)
    fo.close()
