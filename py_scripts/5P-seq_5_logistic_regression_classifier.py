import os
import sys
import re
from sklearn.model_selection import cross_validate, RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression
from collections import Counter
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import precision_score, recall_score, precision_recall_curve, auc
import matplotlib.pyplot as plt
#load data
features = 35
x_list = []
y_list = []
with open(sys.argv[1]) as infile:
    for line in infile:
        fields = line.split("\t")
        x_list.append([float(j) for j in fields[:features]])
        y_list.append(int (fields[features].rstrip()))
#count 1's and 0' for weights
counts = Counter(y_list)
weights = {0: (counts[1]/len(y_list)), 1: (counts[0]/len(y_list))}
#standardize data
scaler = MinMaxScaler ()
x_list_norm  = scaler.fit_transform(x_list)
#Logistic Regression CrossValidate
lg = LogisticRegression(class_weight=weights, random_state=1, penalty='l2', n_jobs = -1, max_iter = 1000)
n_splits , n_repeats = 5, 5
cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=1)
lg_cv = cross_validate(lg, x_list_norm, y_list, scoring= ['precision', 'recall', 'average_precision'], cv=cv, n_jobs=-1, return_estimator = True)
f_name = sys.argv[1].split(".")
files = [f for f in os.listdir('.') if re.match(r'.*p7.*'+ f_name[-3] + r'.dat.36col$', f)]
for i in range(n_splits*n_repeats):
    precisions =[]
    recalls =[]
    aucs = []
    p7_fractions = []
    p7_totals = []
    aucs_mns_baseline = []
    for f in files:
        x_list_test = []
        y_list_test = []
        with open(f) as infile:
            for line in infile:
                fields = line.split("\t")
                x_list_test.append([float(j) for j in fields[:features]])
                y_list_test.append(int (fields[features].rstrip()))
        x_list_test_norm  = scaler.fit_transform(x_list_test)
        counts = Counter(y_list_test)
        p7_1_fraction = counts[1]/len(y_list_test)
        p7_fractions.append(p7_1_fraction)
        p7_totals.append(len(y_list_test))
        precisions.append (precision_score (y_list_test,lg_cv['estimator'][i].predict(x_list_test_norm),pos_label = 1))
        recalls.append (recall_score (y_list_test,lg_cv['estimator'][i].predict(x_list_test_norm),pos_label = 1))
        p,r,t = precision_recall_curve (y_list_test,lg_cv['estimator'][i].predict_proba(x_list_test_norm)[:,1])
        aucs.append (auc(r,p))
        aucs_mns_baseline.append (aucs[-1]-p7_fractions[-1])
    p_train = precision_score (y_list,lg_cv['estimator'][i].predict(x_list_norm),pos_label = 1)
    r_train = recall_score (y_list,lg_cv['estimator'][i].predict(x_list_norm),pos_label = 1)
    p_all_train,r_all_train,t_all_train = precision_recall_curve (y_list,lg_cv['estimator'][i].predict_proba(x_list_norm)[:,1])
    auc_train = auc(r_all_train,p_all_train)
    print (*lg_cv['estimator'][i].coef_[0],
        lg_cv['estimator'][i].intercept_[0],
        len(y_list),
        weights[0],
        lg_cv['test_precision'][i],
        lg_cv['test_recall'][i],
        p_train,
        r_train,
        auc_train,
        auc_train-weights[0],
        *p7_totals,
        *p7_fractions,
        *precisions,
        *recalls,
        *aucs,
        *aucs_mns_baseline,
        sep='\t', file = sys.stdout)
    plt.figure()
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.plot([0, 1], [weights[0], weights[0]], linestyle='--', label='No Skill')
    plt.plot(r_all_train, p_all_train, marker='.', label='Logistic')
    plt.plot(r_train, p_train, marker="o", markersize=20, markeredgecolor="red", markerfacecolor="green")
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()
    plt.title("Trained and tested on " + sys.argv[1] + "." + str(i) + "\n p = " + str(round(p_train,2)) + " r = " + str(round(r_train,2)) + " auc = " + str(round(auc_train,2)) + "\nauc-baseline = " + str(round(auc_train-weights[0],2)) + " baseline = " + str(round(weights[0],2)), fontsize=7)
    plt.savefig (sys.argv[1] + "." + str(i) + ".features.pdf")
    plt.close()