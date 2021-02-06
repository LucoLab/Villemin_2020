from sklearn import clone
import sys
import os
import numpy as np 
import logging
logger = logging.getLogger("Cell2Patients")
from sklearn.preprocessing import Imputer
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 
import warnings

class Cell2Patients():
    def __init__(self, clf, threshold, max_run=None, out_dir="./semi_supervised/", field_separator="\t", increment_rate =10, verbose=True):
        self._clf=clone(clf)
        self._increment_rate=increment_rate
        self._max_run=max_run
        self._out_dir=os.path.realpath(out_dir)
        os.makedirs(out_dir, exist_ok=True)
        self._Xl=[]
        self._Xu=[]
        self._names_l=[]
        self._features_l=[]
        self._names_u=[]
        self._features_u=[]
        self._Y=[]
        self._verbose=verbose
        self._FS=field_separator
        self._leftover=""
        self._thr=threshold
        if self._thr < 0.51 or self._thr > 0.99:
            raise ValueError("Threshold must be between 0.51 and 0.99")  
    
    def _mex(self, txt, end="\n"):
        if self._verbose:
            print(txt, flush=True, end=end)
        else:
            if end == "\n":
                logger.debug("{}{}".format(self._leftover , txt))
                self._leftover=""
            else:
                self._leftover=txt
                              
    
    def import_data(self, input_file, labelled=True):
        if labelled:
            X, features, names, Y=self._Xl , self._features_l, self._names_l , self._Y
        else:
            X, features, names, Y=self._Xu , self._features_u, self._names_u , None
        if not os.path.exists(input_file):
            raise FileNotFoundError("File {} not found!".format(input_file))
        head_n= 2 if labelled else 1
        with open(input_file, "r") as f:
            while head_n != 0:
                line=[ v.strip() for v in f.readline().split(self._FS)]
                if line[0]=="group":
                    for idx in range(1, len(line)):
                        Y.append(line[idx])
                elif line[0]=="name":
                    for idx in range(1, len(line)):
                        names.append(line[idx])
                else:
                    raise ValueError("File {} is not well formatted. Found {} as first field, should be 'name' or 'group'.".format(input_file, line[0]))
                head_n-=1
            line=f.readline()
            while line:
                line=[ v.strip() for v in line.split(self._FS)]
                if len(line) == len(names)+1:
                    features.append(line[0])
                    vals=[]
                    for idx in range(1, len(line)):
                        try:
                            vals.append(float(line[idx]))
                        except ValueError:
                            vals.append("NaN")
                    X.append(vals)
                line=f.readline()
                
    def _initData(self):
        self._class_names=sorted(list(set(self._Y)), reverse=True)
        if len(self._class_names ) != 2:
            raise ValueError("Cell2Patients works only with binary classification.\n{} labels found: {}".format(len(self._class_names), ", ".join(list(self._class_names)) ))
        for i in range(len(self._Y)):
            self._Y[i]=self._class_names.index(self._Y[i])
        if sorted(self._features_l) != sorted(self._features_u):
            raise ValueError("Features in labelled and unlabelled are not the same!")
        self._features=self._features_l.copy()
        for i in range(len(self._features_l)):
            self._features_l[i]=self._features.index(self._features_l[i])
        for i in range(len(self._features_u)):
            self._features_u[i]=self._features.index(self._features_u[i])
        warnings.filterwarnings("ignore")
        self._Y=np.asarray(self._Y)
        imp = Imputer(strategy="mean",verbose=1,axis = 1)
        self._Xl=imp.fit_transform(np.array(self._Xl)[self._features_l,:]).T
        imp = Imputer(strategy="mean",verbose=1,axis = 1)
        self._Xu=imp.fit_transform(np.array(self._Xu)[self._features_u,:]).T
        self._Ypatients=[-1 for _ in self._names_u]
        
    def run(self):
        self._initData()
        self._mex("Starting the tranfer learning of {} cell lines to {} patients".format(
            len(self._names_l), len(self._names_u)))
        self._patients_proba=[]
        self._cell_proba=[]
        running=True
        gen=0
        Xl=self._Xl.copy()
        Y=self._Y.copy()
        self._feature_importances=[]
        transferred=[]
        while running:
            gen+=1
            self._mex("-"*20)
            self._mex("Training the generation number {}".format(gen))
            clf=clone(self._clf)
            clf.fit(Xl, Y)
            self._feature_importances.append(clf.feature_importances_)
            newX, newY , novel = self._evaluate(clf, transferred, gen*self._increment_rate)
            if novel[0] + novel[1] == 0 :
                self._mex("No patients added.")
                running=False
            else:
                self._mex("Adding {} patients, {} novel:".format(len(newY), novel[0]+novel[1]))
                self._mex(" - {} labelled as {} ( {} novel ) ".format(len(newY)-sum(newY), self._class_names[0], novel[0]))
                self._mex(" - {} labelled as {} ( {} novel ) ".format(sum(newY), self._class_names[1], novel[1]))
                Xl=np.concatenate((self._Xl, newX))
                Y=np.concatenate((self._Y, newY))
        self._mex("-"*20)
        self._write_results()
        return self.get_patient_labelled_data()
            
    def _evaluate(self, clf, transferred, limit):
        proba=clf.predict_proba(self._Xu)[:,0]
        self._patients_proba.append(proba.copy())
        self._cell_proba.append(clf.predict_proba(self._Xl)[:,0])
        newX, newY=[],[]
        proba=sorted(enumerate(proba), key=lambda tup : tup[1])
        novel=[0,0]
        n=0
        i=0
        while n < limit and 1-proba[i][1] >= self._thr:
            newX.append(list(self._Xu[proba[i][0],:]))
            newY.append(1)
            if not proba[i][0] in transferred:
                novel[1]+=1
            transferred.append(proba[i][0])
            self._Ypatients[proba[i][0]]=1
            n+=1
            i+=1
        i=len(proba)-1
        n=0
        while n < limit and proba[i][1] >= self._thr:
            newX.append(list(self._Xu[proba[i][0],:]))
            newY.append(0)
            if not proba[i][0] in transferred:
                novel[0]+=1
            transferred.append(proba[i][0])
            self._Ypatients[proba[i][0]]=0
            n+=1
            i-=1
        return np.array(newX), np.array(newY) , novel
    
    def _plot_importances(self, pdf):
        return
    def _plot_proba(self,pdf):
        return
    
    def get_patient_labelled_data(self):
        X=[]
        Y=[]
        for p in range(len(self._names_u)):
            if self._Ypatients[p]!=-1:
                Y.append(self._Ypatients[p])
                X.append(self._Xu[p,:])
        return np.array(X), np.array(Y)
     
    def get_patient_unlabelled_data(self):
        X=[]
        for p in range(len(self._names_u)):
            if self._Ypatients[p]==-1:
                X.append(self._Xu[p,:])
        return np.array(X)
           
    def _write_results(self):
        pdf=PdfPages("{}/graphs.pdf".format(self._out_dir))
        self._plot_importances(pdf)
        self._plot_proba(pdf)
        pdf.close()
        with open("{}/feature_importances.tsv".format(self._out_dir), "w") as ofs:
            ofs.write("Generation\tFeatureName\tImportance\n")
            for gen in range(len(self._feature_importances)):
                for feat in range(len(self._features)):
                    ofs.write("{}\t{}\t{:.3f}\n".format(gen, self._features[feat], self._feature_importances[gen][feat] ))
        with open("{}/patients_probabilities.tsv".format(self._out_dir), "w") as ofs:
            ofs.write("Generation\tPatientName\tProbability_{}\tProbability_{}\n".format(self._class_names[0], self._class_names[1]))
            for gen in range(len(self._patients_proba)):
                for p in range(len(self._names_u)):
                    ofs.write("{}\t{}\t{}\t{}\n".format(gen,self._names_u[p], self._patients_proba[gen][p], 1- self._patients_proba[gen][p] ))
        with open("{}/cells_probabilities.tsv".format(self._out_dir), "w") as ofs:
            ofs.write("Generation\tCellName\tProbability_{}\tProbability_{}\n".format(self._class_names[0], self._class_names[1]))
            for gen in range(len(self._cell_proba)):
                for p in range(len(self._names_l)):
                    ofs.write("{}\t{}\t{}\t{}\n".format(gen,self._names_l[p], self._cell_proba[gen][p], 1- self._cell_proba[gen][p] ))
        
