#!/usr/bin/env python

import re
import sys
import numpy as np
import math
from scipy import optimize
from scipy.stats import pearsonr


#########################
# Specify path, file names and values of experimentally determined variables
#########################

# path, where the input files are stored
path = str(sys.argv[1])

# file, which contains read counts associated with different site types in various binding reactions
reads_file_name = str(sys.argv[2])

# total concentration of RNA used in the binding reaction (nM)
total_RNA_concentration = float(sys.argv[3])

RISC = float(sys.argv[4])

KD_nosite = float(sys.argv[5])

leave_one_out_do = int(sys.argv[6])

#########################
# Get read counts of each site type in binding reactions and dilution_factors (experimental data)
# Get total concentration of each site type in the library (experimental data)
#########################
# reads_ij represents the read counts associated with site type i in a sample j
# they are stored in a tab file, which contains m rows (site types) and n+2 columns
# (n binding reactions + first column describing the site types + second column counting reads in the Input)

read_count_data=open(path+reads_file_name,'r')

line = read_count_data.readline()
dilution_factors_all = []
elements = re.split(r'\t',re.split(r'\n',line)[0])

for i in range(len(elements)):
    if i != 0 and i != 1:
        dilution_factors_all.append(float(elements[i])/100)

read_counts_all = {}
input_concentrations={}
line = read_count_data.readline()

while line:
    site_type=re.split(r'\t',line)[0]
    site_read_counts=[]
    for i in range(len(dilution_factors_all)):
        site_read_counts.append(float(re.split(r'\t',re.split(r'\n',line)[0])[i+2]))
    read_counts_all.update({site_type: site_read_counts})
    input_concentrations.update({re.split(r'\t',line)[0]: float(re.split(r'\t',line)[1])})
    line=read_count_data.readline()
read_count_data.close()

total=sum(input_concentrations[key] for key in input_concentrations.keys())
for key in input_concentrations.keys():
    input_concentrations.update({key: input_concentrations[key] / float(total) * float(total_RNA_concentration)})

#########################
# Definition of xij, which is the concentration of site type i recovered in binding reaction j (mathematical model)
#########################

def xij(a_unbound,input_concentrations,KD,i,j,background):
    # a_unbound = list of a_unbound_j for all binding reactions; has the same length as dilution_factors,
    # input_concentrations = dictionary created in a previous section (experimental data),
    # KD = dictionary of KD_i values for all site types; has m elements,
    # i is the site_type and j is the binding reaction that xij is calculated for,
    # background = total concentration of nonspecifically recovered RNA, assumed to be constant across binding reactions.
    
    S = float(0)
    L = sum(input_concentrations[i] for i in input_concentrations.keys())
    for key in input_concentrations.keys():
        S += a_unbound[j] * input_concentrations[key] / (KD[key] + a_unbound[j])
    
    value = input_concentrations[i] * (a_unbound[j] / (KD[i] + a_unbound[j]) * (1 - background / (L - S)) + background / (L - S))
    return value

#########################
# Calculation of a_unbound[j] 
#########################

# a_unbound is computed at each iteration cycle of the optimization process.

def a_unbound_j(j,dilution_factors,a,input_concentrations,KD):
    
    def equation_a_unbound_j(x):
        return dilution_factors[j]*a - x - sum(x*input_concentrations[z]/(KD[z]+x) for z in KD.keys())

    res = optimize.minimize_scalar(equation_a_unbound_j, bounds=(0, a * float(dilution_factors[j])), method='bounded')
    return res.x

#########################
# Initialization
#########################

# For syntax purposes, theta cannot be a dictionary, but can be a list.
# Therefore, theta will be arranged in a certain order and corresponding indices will be stored in index_list.
index_list = input_concentrations.keys()
index_list.remove('nosite')
index_list.append('nosite')
index_list.append('RISC_stock_concentration')
index_list.append('background')

def initalization(read_counts,dilution_factors):
    # Compute initial guess of theta
    # Each theta_i is initialized as ln of the average enrichment of the site type i in each binding reaction.
    theta_0=[]
    
    for z in range(len(index_list)-3):
        site = index_list[z]
        freq_RISC_site=[]
        for j in range(len(dilution_factors)):
            # Rj is total number of reads in a binding reaction j (experimental data)
            Rj = sum(read_counts[z][j] for z in read_counts.keys())  # total number of reads in a binding reaction j (experimental data)
            freq_RISC_site.append(float(read_counts[site][j]) / float(Rj))
            # average enrichment = average frequency in binding reaction / frequency in input library
        enrichment = np.mean(freq_RISC_site) / (input_concentrations[site] / total_RNA_concentration)
        theta_0.append(np.log(1/enrichment)) # Compute KD as 1/enrichment.
    
    theta_0.append(np.log(KD_nosite))
    
    theta_0.append(np.log(RISC))
    
    # theta_m+2 corresponds to parameter on background and is initiated at 0.1 nM.
    theta_0.append(np.log(0.1))
    
    # Now theta values are partially randomized
    # by adding to each theta a value drawn from a normal distribution with mean 0 and standard deviation 0.1
    theta_0_array=np.array([])
    for i in range(len(theta_0)):
        theta_0[i] = theta_0[i] + np.random.normal(loc=0,scale=0.1,size=None)
        theta_0_array = np.append(theta_0_array, theta_0[i])
    return theta_0_array

#########################
# Defining mathematical model of concentrations recovered for all m site types in all n binding reactions 
#########################

# The output is a dictionary, which contains predicted concentrations for all m site types in all n binding reactions
# The structure of the dictionary is similar to read_counts:
# each site type is a key and the value is a list of length n of values for different binding reactions.
# The input is a list of parameters theta, which are the goal of optimization procedure:
# KD1 = exp(theta1), KD2 = exp(theta2), ..., KD_m = exp(theta_m)
# a = exp(theta_m+1), where a is stock concentration of RISC
# background = exp(theta_m+2).

def KD_from_theta(theta):
    KD = {}
    for z in range(len(theta) - 2):
        KD.update({index_list[z]: np.exp(theta[z])})
    return KD

def model(theta):
    # Step 1: compute a_unbound. 
    # This requires knowing a (stock concentration of RISC) and KDs. Define them first from theta.
    KD = KD_from_theta(theta)
    a = np.exp(theta[len(theta) - 2])
    background = np.exp(theta[len(theta) - 1])
    
    a_unbound = []
    for j in range(len(dilution_factors)):
        a_unbound.append(a_unbound_j(j,dilution_factors,a,input_concentrations,KD))
    
    # Step 2: calculate xij for all the site types in all binding reactions. This will yield 
    # a dictionary of m keys and n-long lists of values.
    model_prediction = {}
    for i in KD.keys():
        Xi=[]
        for j in range(len(dilution_factors)):
            Xi.append(xij(a_unbound,input_concentrations,KD,i,j,background))
        model_prediction.update({i: Xi})
    return model_prediction

#########################
# Definition of fcost 
#########################

def fcost(theta_array, scalar=1):
    theta=theta_array.tolist()
    model_prediction = model(theta) # compute xij
    R=[]  # Rj is total number of reads of all m sites in a binding reaction j (experimental data)
    X=[]  # Xj is total concentration of all m sites in a binding reaction j (mathematical model)
    S=[]  # Sj is sum (over all site types) of reads_ij*ln(xij)
    for j in range(len(dilution_factors)):
        R.append(sum(read_counts[key][j] for key in read_counts.keys()))
        X.append(sum(model_prediction[i][j] for i in model_prediction.keys()))
        S.append(sum(read_counts[i][j]*np.log(model_prediction[i][j]) for i in read_counts.keys()))
    return scalar * sum(R[j]*np.log(X[j])-S[j] for j in range(len(dilution_factors)))

#########################
# Definition of gradient
#########################

def gradient(theta_array):
    theta=theta_array.tolist()
    m = len(theta)-2  # get m from len(theta)
    KD = KD_from_theta(theta)
    a = math.exp(theta[len(theta) - 2])
    background = math.exp(theta[len(theta) - 1])
    a_unbound = []
    for j in range(len(dilution_factors)):
        a_unbound.append(a_unbound_j(j,dilution_factors,a,input_concentrations,KD))
    model_prediction = model(theta)  # dictionary of xij
    der=np.array([])
    for k in range(len(theta)):
        # Define delta and interval functions.
        if k < len(theta)-2: # KD
            delta_k_m1 = 0
            delta_k_m2 = 0
            I_k_m = 1
        elif k == len(theta)-2: # RISC stock
            delta_k_m1 = 1
            delta_k_m2 = 0
            I_k_m = 0
        elif k == len(theta)-1: # background
            delta_k_m1 = 0
            delta_k_m2 = 1
            I_k_m = 0
        
        # Start the calculation of d(fcost)/d(theta_k)
        der_k = 0
        for j in range(len(dilution_factors)):
            Sj = sum(a_unbound[j] * input_concentrations[z] / (KD[z] + a_unbound[j]) for z in KD.keys())
            Fj = sum(input_concentrations[z] * KD[z] / (KD[z] + a_unbound[j])**2 for z in KD.keys())
            Rj = sum(read_counts[key][j] for key in read_counts.keys())  # total number of reads in a binding reaction j (experimental data)
            Xj = sum(model_prediction[z][j] for z in model_prediction.keys())  # total concentration of all sites in a binding reaction j (mathematical model)
            
            if k < len(theta)-2:
                site_type_k = index_list[k]
                fkj = input_concentrations[site_type_k] * KD[site_type_k] / (KD[site_type_k] + a_unbound[j])**2
            else:
                fkj = 0
            
            S = 0
            for i in range(m):
                site_type = index_list[i]  # get site type from index_list
                sij = a_unbound[j] * input_concentrations[site_type] / (KD[site_type] + a_unbound[j])
                fij = input_concentrations[site_type] * KD[site_type] / (KD[site_type] + a_unbound[j])**2
                if k == i:
                    delta_k_i = 1
                else:
                    delta_k_i = 0
                
                A = background * (input_concentrations[site_type] - sij) * delta_k_m2
                B = (-a_unbound[j] * fkj * I_k_m + a * dilution_factors[j] * delta_k_m1 * Fj) / (1 + Fj)
                C = -a_unbound[j] * fij * (delta_k_i - fkj * I_k_m / (1 + Fj))
                D = fij * a * dilution_factors[j] * delta_k_m1 / (1 + Fj)
                
                S1 = A + background * (input_concentrations[site_type] - sij) / (total_RNA_concentration - Sj) * B
                S2 = (total_RNA_concentration - Sj - background) * (C + D)
                S += (Rj / Xj - read_counts[site_type][j] / model_prediction[site_type][j]) * (S1 + S2)
            
            der_k += 1 / (total_RNA_concentration - Sj)*S
        
        der = np.append(der, der_k)
        
    return der

#########################
# Prediction of motif counts given parameters values 
#########################

def predict_concentrations(theta_array,input_concentrations,dilution_factors):
    theta=theta_array.tolist()
    m = len(theta)-2  # get m from len(theta)
    KD = KD_from_theta(theta)
    a = np.exp(theta[len(theta) - 2])
    background = np.exp(theta[len(theta) - 1])
    a_unbound = []
    predicted_concentrations = {}
    for j in range(len(dilution_factors)):
        Pr = []
        Rj = sum(read_counts[z][j] for z in read_counts.keys())  # total number of reads in a binding reaction j (experimental data)
        afree = a_unbound_j(j,dilution_factors,a,input_concentrations,KD)
        for i in KD.keys():
            temp = sum(afree * input_concentrations[z] / (KD[z] + afree) for z in KD.keys())
            term1 = background / (sum(input_concentrations[z] for z in KD.keys()) - temp)
            xij = input_concentrations[i] * (afree / (KD[i] + afree) * (1-term1) + term1)
            Pr.append(float(xij))
        Prj = sum(z for z in Pr)
        temp = []
        for m in Pr:
            temp.append(m/Prj*Rj)
        predicted_concentrations.update({dilution_factors[j]: temp})
    predicted_concentrations_transposed = {}
    c = 0
    for i in KD.keys():
        temp = []
        for j in range(len(dilution_factors)):
            temp.append(predicted_concentrations[dilution_factors[j]][c])
        predicted_concentrations_transposed.update({i: temp})
        c += 1
    return predicted_concentrations_transposed


#########################
# Leave-one-out command 
#########################

def leave_one_out(read_count_dict,index_to_remove):
    ans = {}
    for site_type in read_count_dict.keys():
        value = []
        for v in range(len(read_count_dict[site_type])):
            if v != index_to_remove:
                value.append(read_count_dict[site_type][v])
        ans.update({site_type: value})
    return ans

#########################
# Optimization (repeated noptimize times) command 
#########################

noptimize = 0
list_bounds = []
for z in range(len(index_list)-3):
    list_bounds.append((np.log(0.0001),np.log(100)))
list_bounds.append((np.log(0.1),np.log(500))) # nosite
list_bounds.append((np.log(0.1),np.log(100))) # RISC
list_bounds.append((np.log(0.005),np.log(5))) # background

output = open(reads_file_name+'_RISC'+str(RISC)+'_KDnosite'+str(KD_nosite)+'_KD_all','w')
for i in range(len(index_list)):
    output.write(str(index_list[i])+'\t')
output.write('\n')

while noptimize == 0:
    print 'Starting MLE using all binding reactions...'+'\n'
    read_counts = read_counts_all.copy()
    dilution_factors = list(dilution_factors_all)
    for n in range(1):
        theta_0_array = initalization(read_counts,dilution_factors)
        res = optimize.minimize(fcost, theta_0_array, jac = gradient, bounds = list_bounds, method='L-BFGS-B', options={'disp': False})
        if res.success == True:
            prediction = predict_concentrations(res.x,input_concentrations,dilution_factors)
            flag = 1
            for j in range(len(dilution_factors)):
                predicted_counts = []
                sequenced_counts = []
                for i in prediction.keys():
                    predicted_counts.append(prediction[i][j])
                    sequenced_counts.append(read_counts[i][j])
                corr, _ = pearsonr(predicted_counts, sequenced_counts)
                if corr < 0.90:
                    flag = 0
                    print 'Pearson r between the model and the data is below 0.90.'+'\n'
            if flag == 1:            
                noptimize=1
                print 'MLE is finished.'+'\n'
for i in range(len(index_list)):
    output.write(str(np.exp(res.x[i]))+'\t')
output.write('\n')
output.close()
print 'File '+reads_file_name+'_RISC'+str(RISC)+'_KDnosite'+str(KD_nosite)+'_KD_all'+' contains MLE results.'
    
if leave_one_out_do == 1:
    for j in range(len(dilution_factors_all)):
        output = open(reads_file_name+'_RISC'+str(RISC)+'_KDnosite'+str(KD_nosite)+'_KD_wo_'+str(dilution_factors_all[j]*100),'w')
        for i in range(len(index_list)):
            output.write(str(index_list[i])+'\t')
        output.write('\n')
        read_counts = leave_one_out(read_counts_all,j)
        dilution_factors = list(dilution_factors_all)
        print str(dilution_factors[j]*100)+' removed. Starting MLE...'+'\n'
        del dilution_factors[j]
        flag = 0
        while flag == 0:
            theta_0_array = initalization(read_counts,dilution_factors)
            res = optimize.minimize(fcost, theta_0_array, jac = gradient, bounds = list_bounds, method='L-BFGS-B', options={'disp': False})
            if res.success == True:
                prediction = predict_concentrations(res.x,input_concentrations,dilution_factors)
                flag = 1
                for z in range(len(dilution_factors)):
                    predicted_counts = []
                    sequenced_counts = []
                    for i in prediction.keys():
                        predicted_counts.append(prediction[i][z])
                        sequenced_counts.append(read_counts[i][z])
                    corr, _ = pearsonr(predicted_counts, sequenced_counts)
                    if corr < 0.90:
                        flag = 0
                        print 'Pearson r between the model and the data is below 0.90.'+'\n'
        for i in range(len(index_list)):
            output.write(str(np.exp(res.x[i]))+'\t')
        output.write('\n')
        output.close()
        print 'MLE is finished.'+'\n'
        print 'File '+reads_file_name+'_RISC'+str(RISC)+'_KDnosite'+str(KD_nosite)+'_KD_wo_'+str(dilution_factors_all[j]*100)+' contains MLE results.'+'\n'

