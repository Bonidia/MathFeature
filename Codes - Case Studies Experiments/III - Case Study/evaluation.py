import pandas as pd
import numpy as np
import os
import argparse
import sys
from catboost import *
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_predict
from imblearn.metrics import geometric_mean_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import multilabel_confusion_matrix
from sklearn.model_selection import KFold
from catboost import CatBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_validate
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import LinearSVC
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Perceptron
from sklearn.ensemble import BaggingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import confusion_matrix
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import hamming_loss
from sklearn import metrics
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.svm import SVC
from sklearn.metrics import f1_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import matthews_corrcoef
from sklearn import preprocessing 
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import make_scorer
from sklearn.ensemble import StackingClassifier
from sklearn.impute import SimpleImputer
import warnings
warnings.filterwarnings("ignore")


def header(foutput):
	file = open(foutput, 'a')
	file.write("qParameter,Classifier,ACC,std_ACC,SE,std_SE,F1,std_F1,BACC,std_BACC,kappa,std_kappa,gmean,std_gmean")
	file.write("\n")
	return
	
	
def save_measures(classifier, foutput, scores):
	file = open(foutput, 'a')
	file.write("%s,%s,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f" % (i, classifier, scores['test_ACC'].mean(), 
	+ scores['test_ACC'].std(), scores['test_recall'].mean(), scores['test_recall'].std(), 
	+ scores['test_f1'].mean(), scores['test_f1'].std(), 
	+ scores['test_ACC_B'].mean(), scores['test_ACC_B'].std(),
	+ scores['test_kappa'].mean(), scores['test_kappa'].std(),
	+ scores['test_gmean'].mean(), scores['test_gmean'].std()))
	file.write("\n")
	return


def evaluate_model_cross(classifier, model, finput):
	#####################################
	df = pd.read_csv(finput)
	# df = pd.read_csv(finput, header=None)
	X = df[df.columns[1:(len(df.columns) - 1)]]
	print(X)
	# X = df.iloc[:, 1:-1]
	y = df.iloc[:, -1]
	# y = df['label']
	print(y)
	#####################################
	pipe = Pipeline(steps=[
		('StandardScaler', StandardScaler()),
		('clf', model)])
	# scoring = {'ACC': 'accuracy', 'recall': make_scorer(matthews_corrcoef), 'f1': 'f1', 'ACC_B': 'balanced_accuracy', 'kappa': make_scorer(cohen_kappa_score), 'gmean': make_scorer(geometric_mean_score)}
	scoring = {'ACC': 'accuracy', 'recall': 'recall', 'f1': 'f1', 'ACC_B': 'balanced_accuracy', 'kappa': make_scorer(cohen_kappa_score), 'gmean': make_scorer(geometric_mean_score)}
	kfold = KFold(n_splits=10, shuffle=True, random_state=42)
	scores = cross_validate(pipe, X, y, cv=kfold, scoring=scoring)
	save_measures(classifier, foutput, scores)
	y_pred = cross_val_predict(pipe, X, y, cv=kfold)
	conf_mat = (pd.crosstab(y, y_pred, rownames=["REAL"], colnames=["PREDITO"], margins=True))
	# conf_mat = confusion_matrix(y, y_pred)
	print(conf_mat)
	# np.savetxt("scoresACC.csv", scores['test_ACC'], delimiter=",")
	return



##########################################################################
##########################################################################
if __name__ == "__main__":
	print("\n")
	print("###################################################################################")
	print("##########              Author: Robson Parmezan Bonidia                 ###########")
	print("###################################################################################")
	print("\n")
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', help='csv format file, E.g., dataset.csv')
	parser.add_argument('-o', '--output', help='CSV format file, E.g., test.csv')
	args = parser.parse_args()
	finput = str(args.input)
	foutput = str(args.output)
	estimators = [('rf', RandomForestClassifier()),
				  ('Cat', CatBoostClassifier(logging_level='Silent')),
				  ('LR', LogisticRegression()),
				  ('AB', AdaBoostClassifier()),
				  ('KNN', KNeighborsClassifier())]
	experiments = { 
		# "GaussianNB" : GaussianNB(),
		# "DecisionTree" : DecisionTreeClassifier(criterion='gini', max_depth=2, max_leaf_nodes=2, random_state=63),
		# "GradientBoosting" : GradientBoostingClassifier(n_estimators=400, learning_rate=3.0, max_depth=1, random_state=63),
		"RandomForest" : RandomForestClassifier(n_estimators=200, random_state=63),
		# "LogisticRegression" : LogisticRegression(multi_class="multinomial", solver="lbfgs", C=5),
		# "SVM" : svm.SVC(),
		# "Bagging" : BaggingClassifier(svm.SVC(kernel='linear', C=1200, gamma=0.01)),
		# Bagging" : BaggingClassifier(CatBoostClassifier(iterations=500, thread_count=-1, logging_level='Silent')),
		# "KNN" : KNeighborsClassifier(),
		# "Adaboost" : AdaBoostClassifier(),
		# "MLP" : MLPClassifier(),
		# "Catboost" : CatBoostClassifier(thread_count=2, verbose= True),
		# "Catboost" : CatBoostClassifier(iterations=1000, thread_count=-1, logging_level='Silent'),
		# "HistGradientBoosting" : HistGradientBoostingClassifier(random_state=63),
		# "Stacking" : StackingClassifier(estimators = estimators, final_estimator = svm.SVC())
	}
	header(foutput)
	for i in np.arange(6.0, 6.1, 1.0):
		for classifier, model in experiments.items():
			# evaluate_model_holdout_tuning(classifier, model, finput)
			evaluate_model_cross(classifier, model, finput)
			# evaluate_model_holdout(classifier, model, finput, finput_two)
			# evaluate_model_holdout_multi(classifier, model, finput)
##########################################################################
##########################################################################
