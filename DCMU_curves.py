# -*- coding: utf-8 -*-
"""
Created on Thu May 18 08:45:55 2023

Input: CSV file with following columns:
    Time, s
    Time, ms
    raw (1 per replicate)
    offset (1 per replicate, with Fo set at the origin)

Output: integral: dictionary for each sample containing the area above the curve for all replicates
        percentages: dictionary for each sample containing the average percentage (WT set to 100) and associated error
        Raw_full.svg: Graph with full fluorescence traces from PAM data
        Raw&averages.svg: Graph with normalised and offset traces. Panel1 raw, panel2 averages+- st err
        AreaAbove.svg: Example of the method
        PieCharts.svg: Pie charts with WT=100%
    
References:
    Malkin method - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8940271/
    Error propagation - https://www.geol.lsu.edu/jlorenzo/geophysics/uncertainties/Uncertaintiespart2.html
    
@author: sacco004
"""
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import easygui as egui
import os

def integral_above_curve(x,y,lower_xbound,upper_xbound):
    '''Calculates the area above the curve. Takes x and y values of the data 
    and defined boundaries to use for the area calculation. 
    Returns area above curve.'''
    idx = np.where((np.array(x)>=lower_xbound) & (np.array(x)<=upper_xbound))[0] #define x range
    integrated_value = sp.integrate.simpson(x=np.array(x)[idx],y=np.array(y)[idx])
    return integrated_value

def error_propagation(avg_pc,avg_a,err_a,avg_b,err_b):
    '''Calculate error associated to average compared to WT. Takes the average 
    of the sample (calculated as % of WT), raw average of sample (a), 
    error of a, raw average of WT (a), error of WT.
    Returns error after propagation.'''
    err_pc = avg_pc*np.sqrt((err_a/avg_a)**2+(err_b/avg_b)**2) 
    return err_pc

#Import files
path = str(egui.fileopenbox("Select Files to Import"))
current_dir = str(os.path.dirname('/'.join(path.split('\\'))))
os.chdir(current_dir)
base_path = os.getcwd()
files = filter(os.path.isfile, os.listdir( os.curdir))
files = [ff for ff in os.listdir('.') if os.path.isfile(ff)]
data_files = [x for x in files if (x.endswith(".csv")) and ('DCMUattempt' in x)]

#time range of fast kinetics to integrate in ms
x_range=[0,400]

#Store files in a dict and add values for mean and errors
samples=[]
df ={}
for file in data_files:
    name=file.strip('.csv').split('_')[3]   #Retrieve name of the sample
    samples.append(name)
    df[name]=pd.read_csv(file)
    df[name].drop(df[name].tail(10).index,inplace = True)
    #
    for i,offset in enumerate([x for x in df[name].columns.tolist() if 'offset' in x]):
        loc = df[name].set_index('Time, ms').index.get_loc(x_range[1])
        norm = df[name][offset]/df[name][offset][loc]
        df[name][f"norm.{i}"]= norm
    
    df[name]['avg']=df[name][[x for x in df[name].columns if 'norm' in x]].mean(axis=1)
    df[name]['SD']=df[name][[x for x in df[name].columns if 'norm' in x]].std(axis=1)
    df[name]['SE']=df[name]['SD']/np.sqrt(len([x for x in df[name].columns if 'norm' in x]))
    df[name]['upp_bound']=df[name]['avg']+df[name]['SE']
    df[name]['low_bound']=df[name]['avg']-df[name]['SE']

#Define name of samples and associated colors for graphs
color_list = ['red', 'blue','black']
colors = {}
for color,sample in zip(color_list,samples):
    colors[sample]= color
    
#Plot raw complete traces
fig, ax = plt.subplots(1,1,sharey=True,figsize=(8,8))
ax.set_ylabel('Fluorescence, norm',fontsize=14)
ax.set_xlabel('Time, ms',fontsize=14)

for i,dataframe in enumerate(df):
    raw_data = [x for x in df[dataframe].columns if "raw" in x]
    for x,raw_trace in enumerate(raw_data):
        line, = ax.plot(df[dataframe]["Time, ms"],df[dataframe][raw_trace],color=colors[dataframe])
        if x==0:
            line.set_label(samples[i]) #label only one of these plots
ax.legend()
plt.savefig('Raw_full.svg')
plt.show()
    
#Plot normalised traces
fig, ax = plt.subplots(1,2,sharey=True,figsize=(9,4.5))
ax[0].set_ylabel('Fluorescence, norm',fontsize=14)
ax[0].set_xlabel('Time, ms',fontsize=14)
ax[0].set_xlim(x_range)
ax[0].set_ylim(0,1.1)
for i,dataframe in enumerate(df):
    norm_data = [x for x in df[dataframe].columns if "norm" in x]
    for x,norm_trace in enumerate(norm_data):
        line, = ax[0].plot(df[dataframe]["Time, ms"],df[dataframe][norm_trace],color=colors[dataframe])
        if x==0:
            line.set_label(samples[i]) #label only one of these plots
ax[0].legend()
    
#Plot averages and errors
ax[1].set_xlabel('Time, ms',fontsize=14)
ax[1].set_xlim(x_range)
ax[1].set_ylim(0,1.1)
for i,dataframe in enumerate(df):
    ax[1].fill_between(df[dataframe]["Time, ms"],df[dataframe]["upp_bound"],df[dataframe]["low_bound"],color=colors[dataframe],alpha=0.3)
    ax[1].plot(df[dataframe]["Time, ms"],df[dataframe]["avg"],color=colors[dataframe])
plt.tight_layout()
plt.savefig('Raw&averages.svg')
plt.show()

#Plot example area
fig, ax = plt.subplots(figsize=(8,5))
plt.xlim(x_range)
plt.ylim(0,1)
plt.ylabel('Fluorescence, norm',fontsize=14)
plt.xlabel('Time, ms',fontsize=14)
ax.plot(df['WT']["Time, ms"],df['WT']["avg"],color="black")
ax.fill_between(df['WT']["Time, ms"],df['WT']["avg"],1,color='cyan',alpha=0.3)
plt.savefig('AreaAbove.svg')
plt.show()

#integrate below each curve and store in a dictionary
integral={}
for i,dataframe in enumerate(df):
    norm_data = [x for x in df[dataframe].columns if "norm" in x]
    integral.update({dataframe:[]})
    for x,norm_trace in enumerate(norm_data):
        a = 1/(x_range[1] - integral_above_curve(df[dataframe]["Time, ms"],df[dataframe][norm_trace],*x_range))
        integral.update({dataframe:integral[dataframe]+[a]})

#Calculate averages in % and associated errors and store in a dictionary
percentages={}
err_WT = np.std(integral['WT'])/np.sqrt(len(integral['WT']))
for i,dataframe in enumerate(df):
    avg_pc = np.mean(integral[dataframe])*100/np.mean(integral['WT']) #average percentage
    err = np.std(integral[dataframe])/np.sqrt(len(integral[dataframe]))
    err_pc = error_propagation(avg_pc,np.mean(integral[dataframe]),err,np.mean(integral['WT']),err_WT)
    percentages[dataframe]=[avg_pc,err_pc]

#Plot pie charts
fig,ax = plt.subplots(1,len(df),figsize=(3.3*len(df),3)) #Define figure space

for i,dataframe in enumerate(df):
    # Extract data for the pie chart
    values = [100 - percentages[dataframe][0], percentages[dataframe][0]]
    # Plot the pie chart
    patches, texts = ax[i].pie(values, labels=[None,None], colors=['k', colors[dataframe]])
    patches[0].set_visible(False) #hide complementary percentage
    # Set title for the subplot
    ax[i].set_title(dataframe)
plt.savefig('PieCharts.svg')
plt.show()

