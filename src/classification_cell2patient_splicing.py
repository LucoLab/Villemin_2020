
'''

:date: October 22, 2020
:platform: MacOS 

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis:     A cell-to-patient machine learning transfer approach to classified patients using cell lines.  

'''
# ignore all future warnings
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)
from pathlib import Path
import glob
import argparse,textwrap
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import os 
import re
from os import listdir
from os.path import isfile, join
import codecs
import numpy as np
import scipy.stats as stats
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import csv
import scikitplot as skplt
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.datasets import load_iris
from sklearn import tree
import graphviz 
from sklearn.preprocessing import Imputer
import seaborn as sns
import pandas as pd
from itertools import chain
from sklearn.manifold import TSNE
pd.set_option('display.max_rows', 30)
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from matplotlib.legend_handler import HandlerLine2D
from joblib import dump, load
from sklearn.model_selection import RepeatedStratifiedKFold,StratifiedKFold,cross_val_score
from sklearn.ensemble import (RandomTreesEmbedding, RandomForestClassifier,GradientBoostingClassifier)
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import roc_curve,auc, roc_auc_score, f1_score, accuracy_score , confusion_matrix ,classification_report
from scipy import interp
from sklearn.feature_selection import RFECV
from datetime import datetime
from sklearn.preprocessing import Binarizer,KBinsDiscretizer 
import sklearn
import matplotlib
import time
import sys
from collections import Counter
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from boruta import BorutaPy
import copy 
import scipy.cluster.hierarchy as sch
from sklearn.decomposition import IncrementalPCA

###########################################################################################################
########################################   Functions   ####################################################
###########################################################################################################


def wordListToFreqDict(wordlist):
    """
    Count number of occurencies of a list of string.

    Args:
        wordlist(array)     : list of strings

    Returns:
        dict(list(zip(wordlist,wordfreq))) : Dictionary with string associated with number of occurencies.

    """
    wordfreq = [wordlist.count(p) for p in wordlist]
    return dict(list(zip(wordlist,wordfreq)))

def all_features_importances (importances,genes,X) :
    
    indices = np.argsort(importances)[::-1]
    dict_best = {}
    
    # GO throught the whole list of genes.
    for f in range(X.shape[1]):
        genes[indices[f]]
        #print("%d. feature %d %s (%f)" % (f + 1, indices[f],genes[indices[f]], importances[indices[f]]))
        #1. feature 54 SLK_chr10_104010815-104010908 (0.102410)
        dict_best[genes[indices[f]]] = importances[indices[f]]
           
    return sorted(dict_best.items(), key=lambda x: x[1], reverse=True)       
           
def plot_features_importances (importances,clf,genes,id_cv,nbfeaturesOfImportance,pathOutput,tag,X) :
    """
    Plot the features of importance for the classification.
    
    http://scikit-learn.org/stable/auto_examples/ensemble/plot_forest_importances.html#sphx-glr-auto-examples-ensemble-plot-forest-importances-py

    Args:
        importances(object)     : Features of importance return by sklean - clf.feature_importances_
        clf(object)             : model used
        genes(list)             : list of features - genes.
        id_cv                    : tag for the iteration of the cross validation.
        nbfeaturesOfImportance(int) : number of features of importance.
        pathOutput(str)          :  The output path to write the files.
        tag(str)                 :  Tag to annotated the output files.
        X(?)                     : Values for all features group by samples.
    
    
    Returns:
        genesToPlot : Genes of importance.
        
        
    Note : Need to improve the plot.It's buggy depending the number of features, size of the picture.

        
    """
    std = np.std([tree.feature_importances_ for tree in clf.estimators_],axis=0)
    
    # Return indice of the biggest value to the smalest 
    indices = np.argsort(importances)[::-1]

    result = open(pathOutput+tag+"_features"+str(id_cv)+".txt","w")

    genesToPlot = []
    dict_best = {}
    
    # GO throught the whole list of genes. Select by the indices the gene.
    for f in range(X.shape[1]):
        genesToPlot.append(genes[indices[f]])
        if (f < nbfeaturesOfImportance) :
           genes[indices[f]]
            #print("%d. feature %d %s (%f)" % (f + 1, indices[f],genes[indices[f]], importances[indices[f]]))
            #1. feature 54 SLK_chr10_104010815-104010908 (0.102410)
           dict_best[genes[indices[f]]] = importances[indices[f]]
            #result.write(genes[indices[f]]+"\t"+importances[indices[f]]+"\n")
   
    '''
    for index, value in sorted(dict_best.items(), key=lambda x: x[1], reverse=True) : 
        result.write(index+"\t"+str(value)+"\n")
    result.close()
    
    # Plot the feature importances of the forest
    plt.figure() ##figsize=(80,60)
    plt.rc('xtick', labelsize=10) 

    plt.title("Feature importances")
    plt.bar(range(X.shape[1]), importances[indices],color="r", yerr=std[indices], align="center")
    plt.xticks(range(nbfeaturesOfImportance), genesToPlot[:nbfeaturesOfImportance],rotation='vertical')
    plt.margins(0.3)
    
    # Tweak spacing to prevent clipping of tick-labels
    plt.subplots_adjust(bottom=0.60)
    #plt.tight_layout
    plt.xlim([-1, nbfeaturesOfImportance+1])
    plt.show()
    plt.savefig(pathOutput+tag+"_features"+str(id_cv)+".png" )
    plt.close()
    '''
    return sorted(dict_best.items(), key=lambda x: x[1], reverse=True)




def import_data_used2_create_model(ccle,label,zscore_choice):
    """
    Import data from dataframe. Here cell lines.
  
    Args:
        ccle (object): Dataframe of the cellLines
        label (str): Label of the group you want to classify.


    Returns:
        X  (arrayOfarray) : Values organized by sample. 1 array correspond 1 to one sample
        Y (array) :  Labels for each array of X, for each sample
        samplesID (array)  : List of the sampleID.(ordered as the initial matrice)
        genesASarray(array) : List of genes of the whole matrice.  (ordered as the initial matrice)
        
    """
    
    #print(ccle)

    labels    =  []
    labels_raw    =  []

    samplesID =  []
    for name in list(ccle.columns.values) :
        
        if name[3].replace("Group2: ","") == label : labels.append(name[3].replace("Group2: ",""))
        else : labels.append("Other")
        samplesID.append(name[1])
        labels_raw.append(name[3].replace("Group2: ",""))
        
    ccle2 = ccle.reset_index() # makes date column part of your data
    genes = ccle2.level_0
    genesASarray = np.asarray(genes)
    
    # Impute Values you have NA
    #imp = SimpleImputer(strategy="mean",verbose=1,axis = 1)

    #Axis 0 will act on all the ROWS in each COLUMN
    #Axis 1 will act on all the COLUMNS in each ROW
    imp = Imputer(strategy="mean",verbose=1,axis = 1)
    imputed_matrice = imp.fit_transform(ccle.values)
    
    if(zscore_choice != "No") :
        imputed_matrice = zscoreMyStuff (imputed_matrice)
    
    ccle_values_transposed = imputed_matrice.T 
    #print("Values Transposed for each label")
    #print(ccle_values_transposed)
    
    X = ccle_values_transposed
    
    # Convert list to ndarray (necessary to remove error after)
    Y = labels
    Y = np.asarray(Y)
  
    
    return X,Y,samplesID,genesASarray,labels_raw


    
def createDirUsingAPath(path,parameters):
    """
    Create a dir with datetime using a pre-defined path.
  
    Args:
        path (str): Pre-defined path.

    Returns:
        path (str): Path where to save all the png , tsv etc....
    
        
    """
    today = datetime.now()

    pathOutput = path +today.strftime("%Y_%m_%d-%H%M")
    pathOutput = pathOutput +"_"+ str(parameters.randomState)+"_"+ str(parameters.nestimators)+"_"+ str(parameters.maxdepth)
    print(pathOutput)
    try:
        os.mkdir(pathOutput)
    except OSError:
        print ("Creation of the directory %s failed" % pathOutput)
    else:
        print ("Successfully created the directory %s " % pathOutput)
        
    return pathOutput

def show_split(type,label,valueForlabelToControl): 
    """
    Show how many samples are used for cross validation. How many are picked up from each group.
  
    Args:
        type(str): Type of label group were you picked the sample to cross validate . Can be TRAIN or TEST.
        labels(list) : labels in training or labels in test
        valueForlabelToControl (str): The label you want to classify.

     """
   
    countG1 = 0
    countG2 = 0
    #https://stackoverflow.com/questions/34842405/parameter-stratify-from-method-train-test-split-scikit-learn
    for i in label:
        if i == valueForlabelToControl: countG1+=1
        else : countG2+=1
    print(type) 
    print(valueForlabelToControl+" : "+str(countG1)+"/"+str(countG2)+" (basalB/other) in  "+str(len(label)) +"  samples selected.")
 

def train_vs_test_AUC(x_train,y_train,x_test,y_test,classes,index,id_cv,pathOutput,tag) : 
    """
    TODO
   #https://medium.com/@mohtedibf/indepth-parameter-tuning-for-decision-tree-6753118a03c3

    Args:
        

    Returns:
       
     """
    max_depths       = np.linspace(1, 10, 10, endpoint=True)
    #print(max_depths)
    train_results    = []
    test_results     = []
    
    for max_depth in max_depths:
        
       dt = tree.DecisionTreeClassifier(max_depth=max_depth,random_state=5)
       
       dt.fit(x_train, y_train)
       
       ### TRAIN
       predicted_probas = dt.predict_proba(x_train)

       false_positive_rate, true_positive_rate, thresholds = roc_curve(y_train, predicted_probas[:,index],pos_label=classes[index])
       
       roc_auc = auc(false_positive_rate, true_positive_rate)
       # Add auc score to previous train results
       train_results.append(roc_auc)
       
       ### TEST
       predicted_probas = dt.predict_proba(x_test)

       false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, predicted_probas[:,index],pos_label=classes[index])
       
       roc_auc = auc(false_positive_rate, true_positive_rate)
       # Add auc score to previous test results
       test_results.append(roc_auc)
       
    line1, = plt.plot(max_depths, train_results, 'b', label = "Train AUC")
    line2, = plt.plot(max_depths, test_results, 'r', label  = "Test AUC")
    
    plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)})
    plt.ylabel('AUC score')
    plt.xlabel('Tree depth')
    plt.show()
    plt.savefig(pathOutput+tag+"_"+str(id_cv)+".png" )
    plt.close()

def test_and_train (clf,X,Y,nb_split,pathOutput,tag,genes,isfitted) :
    """
    TODO
  
    Args:
        clf         :
        X           :
        Y           :
        nb_split    :
        pathOutput  :
        tag         :
        genes       :

    Returns:
       
     """
    tprs = []
    aucs = []
    tprs_jp = []
    fprs_jp = []
    accuracies = []
    
    mean_fpr = np.linspace(0, 1, 100)    
    id_cv=0
    
    #This cross-validation object is a variation of KFold that returns stratified folds.
    #1 The folds are made by preserving the percentage of samples for each class.
    cv_skf = StratifiedKFold(n_splits=nb_split,random_state=1)
    
    index_features_of_interest = []

    
    for train_index, test_index in cv_skf.split(X,Y):
        
        X_train, X_test = X[train_index], X[test_index]
        Y_train, Y_test = Y[train_index], Y[test_index]
        #print ("TEST AND TRAIN :")
        #print("TRAIN:", train_index, "TEST:", test_index)

        classes = np.unique(Y_test)
        print(classes)
        
        #print(clf)
        if (isfitted=="no") :
           clf.fit(X_train, Y_train)
        
        # all
        importances = clf.feature_importances_

        # Je garde toute les features
        model       = SelectFromModel(clf, prefit=True)# selecting features based on importance weights.
        
        arrayToKeep = model.get_support(indices=True)
        
        index_features_of_interest.append(arrayToKeep)
        
        genesASarray               = np.asarray(genes)
        genes2                     = genesASarray[arrayToKeep]
        number_feature_of_interest = len(genes2)
        
        #print("number_feature_of_interest"+str(number_feature_of_interest))
        
        important_features_dict = {}
        for x,i in enumerate(clf.feature_importances_):
            #print(x,i)
            important_features_dict[x]=i
            
        important_features_list = sorted(important_features_dict,key=important_features_dict.get,reverse=True)

        # Show how the CF worked out. Number of samples picked up in train/test.
        show_split("TRAIN",Y_train,classes[0])
        show_split("TEST",Y_test,classes[0])

        predicted_probas = clf.predict_proba(X_test)
        Y_predict        = clf.predict(X_test)
        
        train_vs_test_AUC(X_train,Y_train,X_test,Y_test,classes,0,id_cv,pathOutput,tag) 
        
        print("==> Accuracy : ")
        print(accuracy_score(Y_test, Y_predict))
        accuracies.append(accuracy_score(Y_test, Y_predict))
        #print("==> Confusion matrice : ")
        conf = pd.DataFrame(confusion_matrix(Y_test, Y_predict),columns=["Predicted "+classes[0], "Predicted "+classes[1]],index=[classes[0], classes[1]])
        #print(conf)
       
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(Y_test, predicted_probas[:, 0],pos_label=classes[0])
        
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        
        roc_auc = auc(fpr, tpr)
        print("==> AUC : ")
        print(roc_auc)
        aucs.append(roc_auc)
        
        skplt.estimators.plot_learning_curve(clf, X_train, Y_train)
        plt.show()
        plt.savefig(pathOutput+tag+"_learning_curve_"+str(id_cv)+".png" )
        plt.close()
        
        skplt.metrics.plot_confusion_matrix(Y_test, Y_predict, normalize=True)
        plt.show()
        plt.savefig(pathOutput+tag+"_confusion_norm_"+str(id_cv)+".png" )
        plt.close()
    
        skplt.metrics.plot_confusion_matrix(Y_test, Y_predict)
        plt.show()
        plt.savefig(pathOutput+tag+"_confusion_"+str(id_cv)+".png" )
        plt.close()
        
        skplt.metrics.plot_roc(Y_test, predicted_probas,plot_micro=False,plot_macro=False,classes_to_plot=[classes[0]])
        plt.show()
        plt.savefig(pathOutput+tag+"_roc_"+str(id_cv)+".png" ,)
        plt.close()
    
        skplt.metrics.plot_precision_recall(Y_test, predicted_probas,plot_micro=False,classes_to_plot=[classes[0]])
        plt.show()
        plt.savefig(pathOutput+tag+"_"+str(id_cv)+".png" ,)
        plt.close()
        id_cv += 1
 
    recap_of_common_important_features(genes,index_features_of_interest,pathOutput,tag)
    
    return accuracies,aucs
    
def recap_of_common_important_features(genes,index_features_of_interest,pathOutput,tag) :    
    """
    Print all the genes in a file that seems to play a role in classification.
  
    Args:
        
        genes (array) : Ordered list of the genes from the whole matrix.
        index_features_of_interest(listofinteger) : List of index to retrieve the different event that can be usefull for classification.
        pathOutput(str) : Path for the ouput.
        tag (str) : Tag to describe the output in the filename.

    """   
       
    index_kept=set()
    for i in index_features_of_interest :
        for v in i :
            index_kept.add(v)
   
    list_index_kept = list(index_kept)
    genesASarray = np.asarray(genes)
    genes2        = genesASarray[list_index_kept]
    
    final = open(pathOutput+tag+"_features_final.txt","w")
  
    print( "\nTotal Gene Kept: "  +str(len(genes2)))
    for f in genes2:
        final.write(str(f)+"\n")
    final.close()
    


def import_data_for_prediction(patients,zscore_choice):
    """
    Import data you want to use to make some predictions on it.
  
    Args:
        patients (object): Dataframe of the patients

    Returns:
        X(arrayOfarray) : Values organized by patients. 1 array correspond 1 to one patient
        patientsID (array)  : List of the patientID.(ordered as the initial matrice)
        genesASarray(array) : List of genes of the whole matrice.  (ordered as the initial matrice)
        
    """

    
    patientsID = []
    patientsMolSubType = []
    
    for name in list(patients.columns.values) :
 
        patientsID.append(name[0])
        patientsMolSubType.append(name[3].replace("Group2:",""))

    patients1 = patients.reset_index() # makes date column part of your data
    genes     = patients1['index']
    genesASarray = np.asarray(genes)

    # Impute Values you have NA
    imp             = Imputer(strategy="mean",verbose=1,axis = 1) # axis = 0  along the column
    #imp             = SimpleImputer(strategy="mean",verbose=1,axis = 1) # axis = 0  along the column

    imputed_matrice = imp.fit_transform(patients.values)
    
    if(zscore_choice != "No") :
        imputed_matrice = zscoreMyStuff (imputed_matrice)

    patients_values_transposed = imputed_matrice.T 

    X = patients_values_transposed
    
    return X,patientsID,genesASarray,patientsMolSubType

def zscoreMyStuff (matrice_values) :
    
    matrice_zcored = stats.zscore(matrice_values, axis=1, ddof=1)

    return matrice_zcored


def read_score(fileToPathScore,list_patient_ordered) : 

    list_patient_ordered_renamed = renamingPatients(list_patient_ordered)
    
    list_event_MES  = [None] * len(list_patient_ordered_renamed)
    list_event_EPI  = [None] * len(list_patient_ordered_renamed)
    
    count = -1
    with open(fileToPathScore) as lines:
        for line in lines:
           if count==-1 :
               count+=1
               continue
             
           elements    = line.strip().split(";") 

           m    = re.search('(.*\.*\..*)\.(.*)$', elements[0])
    
           if(m):

            id_renamed = m.group(1)
            if id_renamed in list_patient_ordered_renamed :

               list_event_MES[list_patient_ordered_renamed.index(id_renamed)]       = elements[1]
               list_event_EPI[list_patient_ordered_renamed.index(id_renamed)]  = elements[3]
           
    lines.close()

    M           =(np.asarray(list_event_MES).astype(np.float)).reshape(1, -1)
    transformer = Binarizer(np.percentile(M[0], 75)).fit(M) 
    best_M      = transformer.transform(M)   
    
    E           =(np.asarray(list_event_EPI).astype(np.float)).reshape(1, -1)
    transformer = Binarizer(np.percentile(E[0], 75)).fit(E) 
    best_E      = transformer.transform(E)   

   
    list_score = [list_event_MES,best_M[0].tolist()]
    
 
    return list_score



def read_claudinLow(fileToPath,list_patient_ordered):
    
    list_patient_ordered_renamed = renamingPatients(list_patient_ordered)

    list_event_CL       = [None] * len(list_patient_ordered_renamed)
    list_event_Stromal  = [None] * len(list_patient_ordered_renamed)
    list_event_Immune   = [None] * len(list_patient_ordered_renamed)
    #list_event_ESTIMATE = [None] * len(list_patient_ordered_renamed) ESTIMATE IS SUM
    list_event_intClust = [None] * len(list_patient_ordered_renamed)
    count = -1
    with open(fileToPath) as lines:
        for line in lines:
           if count==-1 :
               count+=1
               continue
           
           elements    = line.strip().split("\t") 
           
           #print(elements)
           id_renamed= elements[0].replace("-",".")
           if id_renamed in list_patient_ordered_renamed :
              
             
               list_event_CL[list_patient_ordered_renamed.index(id_renamed)]       = elements[30]
               
               list_event_Stromal[list_patient_ordered_renamed.index(id_renamed)]  = elements[31]  if elements[31] != "None"  else 0
               list_event_Immune[list_patient_ordered_renamed.index(id_renamed)]   = elements[32] if elements[32] != "None" else 0
               list_event_intClust[list_patient_ordered_renamed.index(id_renamed)] = elements[34]
           
    lines.close()
    
    list_annot_global = [list_event_CL,list_event_Stromal,list_event_Immune,list_event_intClust]
    
    
    return list_annot_global


def read_clinical_end_points (fileToPath,list_patient_ordered) :
    
    list_patient_ordered_renamed = renamingPatients(list_patient_ordered)
    
    list_event_OS  = [None] * len(list_patient_ordered_renamed)
    list_event_DSS = [None] * len(list_patient_ordered_renamed)
    list_event_DFI = [None] * len(list_patient_ordered_renamed)
    list_event_PFI = [None] * len(list_patient_ordered_renamed)
    count = -1
    with open(fileToPath) as lines:
        for line in lines:
           if count==-1 :
               count+=1
               continue
             
           elements    = line.strip().split(";") 
           id_renamed= elements[0].replace("-",".")
           if id_renamed in list_patient_ordered_renamed :
              

               list_event_OS[list_patient_ordered_renamed.index(id_renamed)] = elements[24]
               list_event_DSS[list_patient_ordered_renamed.index(id_renamed)] = elements[26]
               list_event_DFI[list_patient_ordered_renamed.index(id_renamed)] = elements[28]
               list_event_PFI[list_patient_ordered_renamed.index(id_renamed)] = elements[30]
           
    lines.close()
    
    list_event_global = [list_event_OS,list_event_DSS,list_event_DFI,list_event_PFI]
    
    
    return list_event_global
    
def rewriteWithUpdatedHeader(pathTofile,headerlines,tag,genes,list_event_global,id):

    result      = open(os.path.dirname(pathTofile)+"/"+tag+"_"+id+"_BASAL_HEADER_ADDED.tsv","w")
    
    if(tag=="expression"):
        for headerline in headerlines : 
            result.write("\t"+"\t"+("\t".join(str(x) for x in headerline))+"\n")
        for list_event in list_event_global :
            result.write("\t"+"\t"+("\t".join(str(x) for x in list_event))+"\n")
    else :
        for headerline in headerlines : 
            result.write("\t"+("\t".join(str(x) for x in headerline))+"\n")
        for list_event in list_event_global :
            result.write("\t"+("\t".join(str(x) for x in list_event))+"\n")

    lineCount=0
    with open(pathTofile) as lines:
        
        for line in lines:
           
            if(tag=="expression"):
                
                if (lineCount == 0) :

                    result.write(line)
                    
                elements = line.strip().split("\t") 

                if elements[1] in genes : result.write(line)

            if(tag!="expression"):
                result.write(line)
            lineCount+=1
            
    lines.close()
    
    result.close()

    
def renamingPatients(patientHeaderToRename) : 
    
    patientHeaderRenamed = []
    
    for patient in  patientHeaderToRename   :
     
        m    = re.search('Patient:\s+(TCGA.*)_(.*)_(.*)', patient)
    
        if(m):
            patientHeaderRenamed.append(m.group(1).replace("-","."))
            
    return patientHeaderRenamed


def compute_proba_and_add_best_from_fitted_model(N,clf,X,patientsID,threshold,iteration) :
    """
    Compute_proba_and_add_best_from_fitted_model2
  
    Args:
       N,
       clf,
       X,
       patientsID : ID 
       threshold : cut off for the probability 0.6
       iteration : Number of the run

    Returns:
    
        basalB_probas : All probability to be basal B computed for Patients.
        annotationNewHeader : For the heatmap, each run will add the patients annotated as B or O in a list added to annotationNewHeader
        bestPatientsMES : ID  Patients annotated as MES
        bestPatientsEPI: ID  Patients annotated as EPI
        allPatients : For the heatmap, in each run we modified allPatients list to show if patient has been annotated B or N in at least one of the run
        X[index_best_patient] : List of patients that have been annoted in the 10 B , 10 O passing the threshold
        YlistOfNewLabel : The label Modified for patients ,  Patient are Annotated B or O
     
    """   
    
    count               = 0   
    annotationNewHeader = ["-"] * len(X)
    allPatients         = {}
    allPatients_nothing = {}
    index_best_patient  = []
    global probas_selected 
    print("Patients  {} ".format(len(X)))
    #exit()
    # X  =  PATIENTS. FIRST RUN.
    # X  =  CELL LINES + PATIENTS 
    global n
    n = N * iteration
  
    print("Number we try to add in both :")
    print(n)
    global subset_dict_best_features_of_importance

    YlistOfNewLabel     = []
    basalB_probas             = clf.predict_proba(X)[:,0]

    index               = list(range(0 ,len(X))) #< 1 - threshold
    global indicesRunAddedBasalB
    global indicesRunAddedOther
  
    all_patients   = dict(zip(index, basalB_probas.tolist()))
    
    ######### ANNOTATED THE BAD ONE AS MANY AS THERE IS GOOD PATIENTS #####
    listofTuples = sorted(all_patients.items() ,  key=lambda x: x[1] ,reverse=True)

    number_bad = 0
    for elem in listofTuples[-n:] :#[-n:][-n:]
     
        if(1 - elem[1] < threshold  ) : continue
        
        allPatients[patientsID[elem[0]]] = "Inverse_BASALB"
        # Comment these lines if you want to add them in the model
        index_best_patient.append(elem[0])
        YlistOfNewLabel.append('Other')
        ###
        annotationNewHeader[elem[0]]     = "O" #str(iteration+1) # just had for plotting
        global_annot[elem[0]]            = "Inverse_BASALB"   
        number_bad+=1

    
    number_good = 0
    for elem in listofTuples[:n] : #[:n]
       
        if(elem[1] < threshold  ) : continue
        allPatients[patientsID[elem[0]]] = "BASALB"
        index_best_patient.append(elem[0])
        YlistOfNewLabel.append('BASALB')       # YlistOfNewLabel.append('Group2: BASALB')

        annotationNewHeader[elem[0]]     = "B" #str(iteration) # just had for plotting
        global_annot[elem[0]]            = "BASALB"   
        number_good+=1   
    
    # TRICK 
    for elem in listofTuples :    
        if  (patientsID[elem[0]] not in allPatients) :
            allPatients_nothing[patientsID[elem[0]]] = "-"

    print("{} total pass a threshold ".format(len(allPatients)))
    print("{} BASALB  ".format(sum(value == "BASALB" for value in allPatients.values())))
    print("{} BasalA  .".format(sum(value == "Inverse_BASALB" for value in allPatients.values())))

 
    ### CREATE THE LIST OF PATIENTS YOU WANT TO COMPARE ####
    bestPatientsMES = [k for k,v in allPatients.items() if v == 'BASALB'] 
    bestPatientsEPI = [k for k,v in allPatients.items() if v == 'Inverse_BASALB'] 
    bestPatientsEPI_ALL =  [k for k,v in allPatients_nothing.items() if v == '-'] + [k for k,v in allPatients.items() if v == 'Inverse_BASALB']
    
    bestPatientsMES_renamed     = renamingPatients(bestPatientsMES)
    bestPatientsEPI_renamed     = renamingPatients(bestPatientsEPI)
    bestPatientsEPI_ALL_renamed = renamingPatients(bestPatientsEPI_ALL)
    
    ##### PRINT TO FILE PATIENT MES TO USE TO PLOT SURVIVAL ###
    path2MES = pathOutput+str(iteration)+"_Patients_MES.txt"
    final = open(path2MES,"w")
    for patient in bestPatientsMES_renamed:
        final.write(str(patient.replace(".","-"))+"\n")
    final.close()
    
    ##### PRINT TO FILE PATIENT EPI TO USE TO PLOT SURVIVAL ###
    path2EPI = pathOutput+str(iteration)+"_Patients_EPI.txt"
    final = open(path2EPI,"w")
    for patient in bestPatientsEPI_renamed:
        final.write(str(patient.replace(".","-"))+"\n")
    final.close()

    ##### PRINT TO FILE PATIENT EPI TO USE TO PLOT SURVIVAL ###
    path2EPI_ALL = pathOutput+str(iteration)+"_Patients_EPI_ALL.txt"
    final = open(path2EPI_ALL,"w")
    for patient in bestPatientsEPI_ALL_renamed:
        final.write(str(patient.replace(".","-"))+"\n")
    final.close()

    #### PLOT SURVIVAL ####
    for endPoint in ["OS","DSS","DFI","PFI"]:
        R_command = "Rscript "+Rscript_survival+" -e "+endPoint+" -s "+fileToPath+" -a " +path2MES+ " -b "+path2EPI+" -o "+str(iteration)+"_"+endPoint
        #print(R_command)
        R = subprocess.run((R_command),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        
        R_command2 = "Rscript "+Rscript_survival+" -e "+endPoint+" -s "+fileToPath+" -a " +path2MES+ " -b "+path2EPI_ALL+" -o "+str(iteration)+"_"+endPoint+"_ALL"
        R2 = subprocess.run((R_command2),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)


    #time.sleep(10)

    return basalB_probas,annotationNewHeader,bestPatientsMES,bestPatientsEPI,allPatients,X[index_best_patient],YlistOfNewLabel 



def best_features_of_importance (clf,genesPatients,pathOutput,tag,X,id,i) :
        
    importances = clf.feature_importances_

    model                      = SelectFromModel(clf, prefit=True)# selecting features based on importance weights.
    arrayToKeep                = model.get_support(indices=True)
    genesASarray               = np.asarray(genesPatients)
    genes2                     = genesASarray[arrayToKeep]
    number_feature_of_interest = len(genes2)

     
    ###### Plot the best features selected and return only the best in a dict #####
    subset_sortedDict_bestfeatures = plot_features_importances(importances,clf,np.asarray(genesPatients),id,number_feature_of_interest,pathOutput,tag,X)
    
    all_sortedDict_bestfeatures = all_features_importances (importances,np.asarray(genesPatients),X)
    
    global subset_dict_best_features_of_importance
    global global_dict_best_features_of_importance
    global in_best_features_for_run_num
    
    for index,value in subset_sortedDict_bestfeatures :
        
        if index not in subset_dict_best_features_of_importance : 
            subset_dict_best_features_of_importance[index] = []
        subset_dict_best_features_of_importance[index].append(value)
        
        if index not in in_best_features_for_run_num : 
            in_best_features_for_run_num[index] = []
        in_best_features_for_run_num[index].append(i)
        
    for index,value in all_sortedDict_bestfeatures :
        if index not in global_dict_best_features_of_importance : 
            global_dict_best_features_of_importance[index] = []
        global_dict_best_features_of_importance[index].append(value)    
    #return global_dict_best_features_of_importance

def return_index(bigList,list) : 
    
    indexes = []
    for value in list :
     
       indexes.append(int(bigList.index(value)))
       
    return indexes 

###########################################################################################################
########################################   Main   ####################################################
###########################################################################################################

       
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    
     A cell-to-patient machine learning transfer approach uncovers novel basal-like breast  \  
     cancer prognostic markers amongst alternative splice variants.


    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-c","--celllines",action="store",help="Matrice of Cell lines",required=True,type=str,dest='matrice_cellLines')
    parser.add_argument("-p","--patients",action="store",help="Matrice of Patients",required=True,type=str,dest='matrice_patients')
    parser.add_argument("-t","--threshold",action="store",help="Threshold of the best cell lines. 0 to 1.",required=True,type=float,dest='threshold')
    parser.add_argument("-n","--nestimators",action="store",help="Zscore the values",required=True,type=int,default="No",dest='nestimators')

    parser.add_argument("-r","--randomState",action="store",help="Fix RandomState",required=False,type=int,default=2955,dest='randomState')
    parser.add_argument("-rb","--randomStateBoruta",action="store",help="Fix several RandomState for Boruta",required=False,default="Yes",type=str,dest='randomStateBoruta')
    
    parser.add_argument("-z","--zscore",action="store",help="Zscore the values",required=False,type=str,default="No",dest='zscore')
    parser.add_argument("-m","--maxdepth",action="store",help="Zscore the values",default=0,type=int,dest='maxdepth')
 
    parameters = parser.parse_args()
    
    pathname                = os.path.dirname(sys.argv[0])        
    full_path_script_python = os.path.abspath(pathname)
    
    basepath =  full_path_script_python.replace("src", "")
    basepathforData =  basepath+"data/"

    print("Files will be written in :")

    pathOutput = createDirUsingAPath(basepath,parameters)
    pathOutput = pathOutput+"/"
    logger     = open(basepath+"/logger.txt","a")
    
    logger.write(str(pathOutput)+"\n")
    logger.write(str(parameters.threshold)+"\n")
    logger.write(str(parameters.zscore)+"\n")
    logger.write(str(parameters.randomState)+"\n")
    logger.write(str(parameters.nestimators)+"\n")
    logger.write(str(parameters.matrice_cellLines)+"\n")
    logger.write(str(parameters.matrice_patients)+"\n")
    logger.write(str(parameters.maxdepth)+"\n")
   
    maxdepth =None
    if (str(parameters.maxdepth) != "0") : 
        maxdepth = parameters.maxdepth 
    
    print('The scikit-learn version is {}.'.format(sklearn.__version__))
    print('The matplotlib version is {}. \n'.format(matplotlib.__version__))

    fileToPath              =   basepathforData+"TCGA_CDR.csv"
    fileToPathClaudinLow    =   basepathforData+"TCGA-BRCA_ClaudinLow_sorly.txt"
    fileToPathScore         =   basepathforData+"mes_epi_tcga.csv"
    
    Rscript_survival  = full_path_script_python+"/R/"+"Survival_AvsB.R"

    ccle     = pd.read_csv(parameters.matrice_cellLines,sep='\t',header=[0,1,2,3],index_col=[0,1,2])
    
    patients = pd.read_csv(parameters.matrice_patients,sep='\t',header=list(range (0,32,1)),index_col=0)

    XcellLines,YcellLines,samplesCellLinesName,genesCellLines,labelcellLines_raw = import_data_used2_create_model(ccle,"BASALB",parameters.zscore)
    XPatients,patientsID,genesPatients,patientsMolSubType                        = import_data_for_prediction(patients,parameters.zscore)
    
    list_event_global      = read_clinical_end_points(fileToPath,patientsID)
    list_annotation_global = read_claudinLow(fileToPathClaudinLow,patientsID)
    
    list_scoreMesEPi = read_score(fileToPathScore,patientsID)

    mainlistAnnotations = list_event_global + list_annotation_global + list_scoreMesEPi
    
    concat_func = lambda x,y: x + "_" + str(y)

    ######################################################################################################################################################################################
    ######################################################################################################################################################################################
    ######################################################################################################################################################################################
    
    if ( parameters.randomState  == False ) :
      random_state = np.random.randint(low=0, high=10000, size=1)[0]
    else : random_state = parameters.randomState  
    
    clf = ""


    clf = RandomForestClassifier(max_depth=maxdepth,n_estimators=parameters.nestimators,n_jobs=-1,class_weight="balanced",random_state=random_state )#random_state=13#random_state=5,,class_weight="balanced"

    # Global variable
    global_annot =  ["UNKNOW"] * len(XPatients)  #list(range(len(XPatients)))

    # Fit model  on Cell Lines
    print("\nFit model simply on all the cell lines : ")
    clf.fit(XcellLines,YcellLines)
   
    
    global_new_header      = []
    global_bestCandidates  = []
    patients_added_perRun  = []
    global_dict_best_features_of_importance = {} # use as global
    subset_dict_best_features_of_importance = {} # use as global
    in_best_features_for_run_num = {} # use as global
    accuracies_per_run    = {} 
    probas_selected       = [] # use as global
    indicesRunAddedBasalB = [[]]
    indicesRunAddedOther  = [[]]
    all_probas_per_run    = []
    allproba_per_run1  = {} 
    bestPatientsMES_final = 0
    
    patients_added     =[ [0,0]]
    t       = 100000
    N       = 10 # Number of patients. N * RunIter Patients added at each Run.
    run        = 0
    run_stop = 0
    for i in range(1,t) :
        
        print("\n=====> Incremental analysis id # : {} \n".format(i))
      
        all_probas,annotationNewHeader,bestPatientsMES,bestPatientsEPI,allPatients,Xnew,Ynew = compute_proba_and_add_best_from_fitted_model(N,clf,XPatients,patientsID,parameters.threshold,i)
        all_probas_per_run.append(all_probas)
        
        bestPatientsMES_final = bestPatientsMES
        
        indexes = return_index(patientsID,bestPatientsMES)

        run = i 
        
        print("allPatients : ")
        print(len(allPatients))
  
        print("BestCandidates annotated in Patients and added to Cell Lines : {} ".format(len(allPatients)))
        patients_added.append([len(bestPatientsMES),len(bestPatientsEPI)])
        
        if (len(bestPatientsMES) == (patients_added[i-1][0])) :
            print("You are not adding more Basal B.")
            #run_stop = i
        
        if (len(allPatients) == (patients_added[i-1][0]+patients_added[i-1][1])) :
            print("You are not adding more patient in both group.")
            run_stop = i # Run_stop use for diff absolute but it's deprecated
            break
        
        # We add the best to the model
        X2test = np.append(XcellLines,Xnew,axis=0) #CellLines + Patients 
        Y2test = np.append(YcellLines,Ynew,axis=0) #CellLines + Patients 

        clf = RandomForestClassifier(max_depth=maxdepth, min_samples_split = 0.1,n_estimators=parameters.nestimators,n_jobs=-1,class_weight="balanced",random_state=random_state)##random_state=5,
        clf.fit(X2test,Y2test)

        # Call it once be carefull if not this will create buggy results
        best_features_of_importance(clf,genesPatients,pathOutput,"cellLinesAndPatients",X2test,i*10,i)
     
        global_new_header.append(annotationNewHeader)
        
        if(run == run_stop) : break # If you do that clf is not 11 but 12....hummm
    
    logger.write("BasalFinal : "+str(len(bestPatientsMES_final))+"\n")
    
    global_new_header.append(global_annot)
    
    print("\n=====> Feature selection with boruta algorithm \n")

    random_state2 = []
    if (parameters.randomStateBoruta != "No"):
        random_state2 = [2274,931,3891,2845,6538,7524,5051,6298,877,7403]
        #random_state2 = [4358,7315,3137,4079,8288,7711,7987,8850,7823,9659]
    else :     random_state2 = np.random.randint(low=0, high=10000, size=10)
    
    veryBestof  = open(pathOutput+"outputBorutaPy.txt","w")
    rsgenerated = open(pathOutput+"rsBorutaPy.txt","w")

    allgeneskeep = []
    for rs in random_state2 :
        
        clf2 = copy.deepcopy(clf) 
        rsgenerated.write(str(rs)+"\n")
        feat_selector    = BorutaPy(clf2, n_estimators='auto', verbose=0, random_state=rs)
        feat_selector.fit(Xnew,Ynew)
        geneskeep        = genesPatients[feat_selector.support_]
        for i in geneskeep :  allgeneskeep.append(i)
    
    gene_to_count = wordListToFreqDict(allgeneskeep)    
    

    bedSubset = open(pathOutput+"outputBorutaPy.bed","w")
    
    veryBestofs = []
    # Should be in the subbset of best feature
    for gene2keep in sorted(list(set(allgeneskeep))):
        rsgenerated.write(gene2keep+"\t"+str(gene_to_count[gene2keep])+"\n")
        if (gene_to_count[gene2keep] >= 7 ) :
            if (gene2keep in subset_dict_best_features_of_importance ) :
                veryBestof.write(gene2keep+"\n")
                veryBestofs.append(gene2keep)
                id  = "NA_"+str(gene2keep.split("_")[0])+"_NA_NA_NA"
                pos = str(gene2keep.split("_")[1])
                chr = pos.split(":")[0]
                coords = pos.split(":")[1]
                coordStart =  coords.split("-")[0]
                coordEnd   =  coords.split("-")[1]
                bedSubset.write(chr+"    "+coordStart+"    "+coordEnd+"    "+id+"    "+"0"+"    "+"."+"\n")
    veryBestof.close()  
    bedSubset.close()  
    rsgenerated.close()

    indexes = return_index(genesPatients.tolist(),veryBestofs )#np.asarray()
    X_bestfeatureboruta = Xnew[:, indexes]

    #############################################################################################
    #############################################################################################
    #####################           Rscript Plotting                        #####################
    #############################################################################################
    #############################################################################################

    num_of_best_features = len(global_dict_best_features_of_importance) 

    ### WRITE THE PROBABILITY OF EACH RUN ####
    probaFile = open(pathOutput+"/probas.txt","w")
    a = 1
    for run_per_proba in all_probas_per_run :
        for proba in run_per_proba :
            probaFile.write(str(a)+"\t"+str(proba)+"\n")
        a+=1
    probaFile.close()
    
    
    ### WRITE THE Percentage of occurance between the run ####
    occ = open(pathOutput+"/occurrencies.txt","w")
    for gene in in_best_features_for_run_num:
        for value_run in  in_best_features_for_run_num[gene] :
            occ.write(gene+"\t"+str(value_run)+"\n")
    occ.close()  
    
    ### WRITE THE NUMBER OF PATIENTS ADDED ####
    patientsFile = open(pathOutput+"/patients_added.txt","w")
    o = 0
    tot = len(patients_added)
    for add in patients_added:
        if o == 0 : 
            o+=1
            continue
        if o == tot : break 
        patientsFile.write(str(o)+"\t"+"BASALB-LIKE"+"\t"+str(add[0])+"\n")
        patientsFile.write(str(o)+"\t"+"BASALA-LIKE"+"\t"+str(add[1])+"\n")

        o+=1
    patientsFile.close()    
    
    ### WRITE THE FEATURES AND COMPUTE SUM FOR ALL FEATURE THAT YOU WILL WRITE TO ANOTHER FILE LATER ####
    final = open(pathOutput+"/features_final.txt","w")
    final2 = open(pathOutput+"/features_final_normalised.txt","w")

    dict_to_order_by_meanvalue_of_features = {}
    for gene in global_dict_best_features_of_importance:
        dict_to_order_by_meanvalue_of_features[gene] = np.sum(global_dict_best_features_of_importance[gene])
        index = 0 
        mean = np.mean(global_dict_best_features_of_importance[gene])
        if(gene in subset_dict_best_features_of_importance) :
            for feature in global_dict_best_features_of_importance[gene] : 
                
                bestBoruta = "NO"
                if (gene in veryBestofs) : bestBoruta = "YES"
                final.write(str(index+1)+"\t"+gene+"\t"+str(feature)+"\t"+bestBoruta+"\n")

                bestBoruta = "NO"
                if (gene in veryBestofs) : bestBoruta = "YES"
                if (feature==0) :
                    final2.write(str(index+1)+"\t"+gene+"\t"+str(np.log(1))+"\t"+bestBoruta+"\n")
                    continue
                
                if (global_dict_best_features_of_importance[gene][0]==0) : # Don't plot the weirds where the ratio will fail
                    continue
                    
                final2.write(str(index+1)+"\t"+gene+"\t"+str(np.log((feature/global_dict_best_features_of_importance[gene][0])+1))+"\t"+bestBoruta+"\n")

                index+=1
    final.close()  
    final2.close()  

    ### WRITE THE DIFF FROM THE START POINT AND THE STOP ####
    # Be carefull with run_stop
    final3 = open(pathOutput+"/features_final_absolute_diff.txt","w")
    index = 0 
    #print(run_stop)
    for gene in global_dict_best_features_of_importance:  
        if(gene in subset_dict_best_features_of_importance) :
            for feature in global_dict_best_features_of_importance[gene] :
                bestBoruta = "NO"
                if (gene in veryBestofs) : bestBoruta = "YES"
                final3.write(str(index+1)+"\t"+gene+"\t"+str(global_dict_best_features_of_importance[gene][run_stop-2] - global_dict_best_features_of_importance[gene][0])+"\t"+bestBoruta+"\n")
                break
        index+=1
    final3.close()

    ### WRITE THE SUM OF ALL FEATURES ####
    final2 = open(pathOutput+"/bed_features_final.bed","w")
    final7 = open(pathOutput+"/symbol_features_final.txt","w")

    final1 = open(pathOutput+"/sum_features_final.txt","w")
    final6 = open(pathOutput+"/sum_features_final_selectedatLeastOnceInaRun.txt","w")

    p = 1
    genes_set3= []
    for index, value in sorted(dict_to_order_by_meanvalue_of_features.items(), key=lambda x: x[1], reverse=True) :
        bestBoruta = "NO"
        if (index in veryBestofs) : bestBoruta = "YES" 
        final1.write(index+"\t"+str(value)+"\t"+bestBoruta+"\n")
        
        if(index in subset_dict_best_features_of_importance) :
            bestBoruta = "NO"
            if (index in veryBestofs) : bestBoruta = "YES"
            final6.write(index+"\t"+str(value)+"\t"+bestBoruta+"\n")

        if ( p<= num_of_best_features):

            final7.write(str(index.split("_")[0])+"\n")
            genes_set3.append(str(index.split("_")[0])) 
            id  = "NA_"+str(index.split("_")[0])+"_NA_NA_NA"
            pos = str(index.split("_")[1])
            chr = pos.split(":")[0]
            coords = pos.split(":")[1]
            coordStart =  coords.split("-")[0]
            coordEnd   =  coords.split("-")[1]
            final2.write(chr+"    "+coordStart+"    "+coordEnd+"    "+id+"    "+"0"+"    "+"."+"\n")

        p+=1
    final1.close()
    final2.close()
    final6.close()
    final7.close()

    #####################    #####################
    #####################    #####################
    
    # WRITE THE SPLICING OF TCGA PATIENT #####
    #### This is crapy ####
    patients.to_csv(path_or_buf=pathOutput+"_NEW_BASALB.csv",sep="\t")
    # I reopen the file and I do it badely.
    id="TCGA"
    rewriteWithUpdatedHeader(pathOutput+"_NEW_BASALB.csv",global_new_header,"splicing",[],mainlistAnnotations,id)
    # Remove this copy on disk
    remove ="rm "+pathOutput+"_NEW_BASALB.csv"
    removeCommand = subprocess.run((remove),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    print(removeCommand)

    print("random_state Splicing : "+str(random_state))
    rs = open(pathOutput+"/random_state.txt","w")
    rs.write(str(random_state)+"\n")
    rs.close()  
    
    print ("The values are normalised ? : {} \n" .format(parameters.zscore))
    logger.write("+++++++++++"+"\n")
    logger.close()
  
  
