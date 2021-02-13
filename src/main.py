#!/usr/bin/env python3
import argparse
import sys
import os
import textwrap
import logging
logger = logging.getLogger("Cell2Patients")
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
     A cell-to-patient machine learning transfer approach uncovers novel basal-like breast  \  
     cancer prognostic markers amongst alternative splice variants.
    '''), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-c", "--celllines", action="store", help="Matrix of Cell lines data", required=True, type=str, dest='matrice_cellLines')
    parser.add_argument("-p", "--patients", action="store", help="Matrix of Patients data", required=True, type=str, dest='matrice_patients')
    parser.add_argument("-o", "--output", action="store", help="Output folder", default="./results/", type=str, dest='output')
    parser.add_argument("-t", "--threshold", action="store", help="Patients with a score above the threshold will be assigned to one of the two classes. 0.51 to 0.99.", required=False, default="0.6" , type=float, dest='threshold')
    parser.add_argument("-i", "--increment-rate", help="Select the number of best patients to add with a rate of increment_rate * iteration.", type=int, default="10")
    parser.add_argument("-n", "--nestimators", action="store", help="Number of trees in the random forest classifier", required=False, type=int, default="1000", dest='nestimators')
    parser.add_argument("-m", "--maxdepth", action="store", help="Random forest max depth (-1 == None )", default="-1", type=int, dest='maxdepth')
    parser.add_argument("-M", "--max-run", action="store", help="Maximum iterations for the learning transfer", default=1000, type=int, dest='max_run')
    parser.add_argument("-r", "--randomState", action="store", help="Fix RandomState.", required=False, type=int, default=2955, dest='randomState')
    parser.add_argument("--scale", action="store_true", help="Scale the features of cell lines and patients independently", default=False, required=False)
    parser.add_argument("-b", "--boruta-iterations", action="store", help="Repeat the boruta feature selection multiple times. If random state is fixed, boruta random states will be seeded with it. ", required=False, default="10", type=int, dest='boruta')
    parser.add_argument("-s", "--boruta-elite", action="store", help="Number of times a feature has to be selected to be considered elite. If between 0 and 1, considered as fraction of the total boruta iterations.", required=False, default=0.7, type=float, dest='boruta_sel')
    parser.add_argument("--debug", action="store_true", help="Debug mode", required=False)
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose", default=False, required=False)
    
    parameters = parser.parse_args()
    verbose=parameters.verbose
    if parameters.debug:
        logger.setLevel(logging.DEBUG)
        sh=logging.StreamHandler()
        sh.setFormatter(logging.Formatter('%(asctime)s -%(levelname)s - %(message)s'))
        logger.addHandler(sh)
        import sklearn
        logger.debug("Scikit-learn version: {} ".format(sklearn.__version__))
    else:
        logging.disable(logging.WARNING)
    #### setup the environment
    max_run=parameters.max_run
    pathOutput = os.path.realpath(parameters.output)
    if parameters.threshold < 0.51 or parameters.threshold > 0.99:
        raise ValueError("Threshold must be between 0.51 and 0.99")
    if not os.path.exists(pathOutput):
        os.makedirs(pathOutput)
    
    
    maxdepth = None
    if parameters.maxdepth != -1 : 
        maxdepth = parameters.maxdepth
    if (parameters.randomState == False) :
        random_state = np.random.randint(low=0, high=10000, size=1)[0]
    else : 
        random_state = parameters.randomState
    boruta_iterations= parameters.boruta
    np.random.seed(random_state)
    if parameters.boruta == 10 and parameters.randomState != False:
        boruta_random_state=[2274,931,3891,2845,6538,7524,5051,6298,877,7403]
    else:
        boruta_random_state = np.random.randint(low=0, high=10000, size=int(parameters.boruta))
    if parameters.boruta_sel > 0 and parameters.boruta_sel < 1:
        elite_sel=int(np.ceil(parameters.boruta_sel *  len(boruta_random_state)))
    else:
        elite_sel=int(parameters.boruta_sel)
    ### Logging
    mess="Output directory: {} \n".format(pathOutput)
    mess+="Assignment threshold: {}\n".format(parameters.threshold)
    mess+="Maximum number of runs: {}\n".format(max_run)
    mess+="Increment rate: {}\n".format(parameters.increment_rate)
    mess+="Number of trees: {}\n".format(parameters.nestimators)
    mess+="Max depth: {}\n".format(parameters.maxdepth)
    mess+="Scale: {}\n".format(parameters.normalize)
    mess+="Cell line matrix: {}\n".format(parameters.matrice_cellLines)
    mess+="Patients matrix: {}\n".format(parameters.matrice_patients)
    mess+="Boruta iterations: {}\n".format(boruta_iterations)
    mess+="Boruta elite selection: {}\n".format(elite_sel)
    mess+="Random state: {}\n".format(random_state)
    mess+="Boruta random states: {}\n".format(boruta_random_state)
    
    
    if verbose:
        print(mess)
    logger.debug(mess)
    
    ### Start the job
    from Cell2Patients import Cell2Patients
    from sklearn.ensemble import RandomForestClassifier
    from sklearn import clone
    from boruta import BorutaPy
    
    
    clf = RandomForestClassifier(max_depth=maxdepth,n_estimators=parameters.nestimators,min_samples_split=2, n_jobs=-1,class_weight="balanced",random_state=random_state)
    clf_t = RandomForestClassifier(max_depth=maxdepth,n_estimators=parameters.nestimators,min_samples_split = 0.1,n_jobs=-1,class_weight="balanced",random_state=random_state)
    if verbose:
        print("Reading the data...", end="", flush=True)
    c2p = Cell2Patients(clf=clf, clf_transfer= clf_t,max_run=max_run,threshold=parameters.threshold,increment_rate=parameters.increment_rate, out_dir="{}/transfer/".format(pathOutput), verbose=verbose, normalize=parameters.normalize)
    c2p.import_data(parameters.matrice_cellLines, labelled=True)
    c2p.import_data(parameters.matrice_patients, labelled=False)
    if verbose:
        print("Done.\nRunning the cell to patient transfer learning...", flush=True)
    X, Y=c2p.run()
    if verbose:
        print("Done.", flush=True)
        mess="Found {} patients in the class {} and {} in {}".format(len(Y)-sum(Y), c2p._class_names[0],sum(Y), c2p._class_names[1] )
        logger.debug(mess)
        print(mess)
        print("Running the boruta feature selection on the new data")
    if sum(Y) == 0 or sum(Y) == len(Y):
        print("Impossible to run boruta with patients from just a single class.")
        return
    supports=[]
    ranking=[]
    for it, rs in enumerate(boruta_random_state):
        if verbose:
            print("Iteration {}...".format(it), end="", flush=True)
        boruta = BorutaPy(clone(clf_t), n_estimators='auto', verbose=0, random_state=rs)
        boruta.fit(X,Y)
        if verbose:
            logger.debug("boruta iteration {} : identified {} features.".format(it, sum(boruta.support_)))
            print("identified {} features.".format(sum(boruta.support_)), flush=True)
        supports.append(boruta.support_)
        ranking.append(boruta.ranking_)
    os.makedirs("{}/boruta/".format(pathOutput),exist_ok=True )
    with open("{}/boruta/all_features_ranking.tsv".format(pathOutput), "w") as oa:
        oa.write("Feature")
        for r in range(len(ranking)):
            oa.write("\trank_it_{}".format(r))
        oa.write("\n")
        with open("{}/boruta/elite_features.tsv".format(pathOutput), "w") as ob:
            ob.write("Feature\tSuports\n")
            for fidx, feat in enumerate(c2p._features):
                supp=0
                oa.write(feat)
                for r in range(len(ranking)):
                    oa.write("\t{}".format(ranking[r][fidx]))
                    if supports[r][fidx]:
                        supp+=1
                oa.write("\n")
                if supp >= elite_sel:
                    ob.write("{}\t{}\n".format(feat, supp))
    print("Cell2Patients finished successfully. You can find your results in {}.\n".format(pathOutput))
    
    
if __name__ == "__main__":
    main();
