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
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler
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


def evaluate_model_holdout(classifier, model, finput, finput_two):
	df1 = pd.read_csv(finput, header=None)
	train_labels = df1.iloc[:, -1]
	train = df1[df1.columns[1:(len(df1.columns) - 1)]]
	print(train)
	print(train_labels)

	df2 = pd.read_csv(finput_two, header=None)
	test_labels = df2.iloc[:, -1]
	test = df2[df2.columns[1:(len(df1.columns) - 1)]]
	print(test)
	print(test_labels)

	sc = StandardScaler()
	train = sc.fit_transform(train)
	test = sc.transform(test)

	# sm = SMOTE(random_state=42)
	sm = RandomUnderSampler(random_state=42)
	train, train_labels = sm.fit_sample(train, train_labels)

	print("Amount of train: " + str(len(train)))
	print("Amount of test: " + str(len(test)))

	clf = model
	kfold = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
	scores = cross_validate(clf, train, train_labels, cv=kfold, scoring='accuracy')
	clf.fit(train, train_labels)
	preds = clf.predict(test)
	accu = accuracy_score(test_labels, preds)
	recall = recall_score(test_labels, preds)
	f1 = f1_score(test_labels, preds)
	balanced = balanced_accuracy_score(test_labels, preds)
	gmean = geometric_mean_score(test_labels, preds)
	mcc = matthews_corrcoef(test_labels, preds)
	matriz = (pd.crosstab(test_labels, preds, rownames=["REAL"], colnames=["PREDITO"], margins=True))
	print("Classificador: %s" % (classifier))
	print("Predições %s" % (preds))
	print("Train Score (kfold=10): %s" % scores['test_score'].mean())
	print("Acurácia Teste: %s" % (accu))
	print("Recall: %s" % (recall))
	print("F1: %s" % (f1))
	print("balanced: %s" % (balanced))
	print("gmean: %s" % (gmean))
	print("MCC: %s" % (mcc))
	print("%s" % (matriz))
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
	parser.add_argument('-t', '--input_t', help='csv format file, E.g., dataset.csv')
	args = parser.parse_args()
	finput = str(args.input)
	finput_two = str(args.input_t)
	estimators = [('rf', RandomForestClassifier()),
				  ('Cat', CatBoostClassifier(logging_level='Silent')),
				  ('LR', LogisticRegression()),
				  ('AB', AdaBoostClassifier()),
				  ('KNN', KNeighborsClassifier())]
	experiments = { 
		# "GaussianNB" : GaussianNB(),
		# "DecisionTree" : DecisionTreeClassifier(criterion='gini', max_depth=2, max_leaf_nodes=2, random_state=63),
		# "GradientBoosting" : GradientBoostingClassifier(n_estimators=400, learning_rate=3.0, max_depth=1, random_state=63),
		"RandomForest" : RandomForestClassifier(n_estimators=1000, random_state=63),
		# "LogisticRegression" : LogisticRegression(multi_class="multinomial", solver="lbfgs", C=5),
		# "SVM" : svm.SVC(),
		# "Bagging" : BaggingClassifier(svm.SVC(kernel='linear', C=1200, gamma=0.01)),
		# Bagging" : BaggingClassifier(CatBoostClassifier(iterations=500, thread_count=-1, logging_level='Silent')),
		# "KNN" : KNeighborsClassifier(),
		# "Adaboost" : AdaBoostClassifier(),
		# "MLP" : MLPClassifier(),
		# "Catboost" : CatBoostClassifier(thread_count=2, verbose= True),
		# "Catboost" : CatBoostClassifier(iterations=1000, thread_count=-1, logging_level='Silent', random_state=63),
		# "HistGradientBoosting" : HistGradientBoostingClassifier(random_state=63),
		# "Stacking" : StackingClassifier(estimators = estimators, final_estimator = svm.SVC())
	}
	for i in np.arange(6.0, 6.1, 1.0):
		for classifier, model in experiments.items():
			evaluate_model_holdout(classifier, model, finput, finput_two)
##########################################################################
##########################################################################
