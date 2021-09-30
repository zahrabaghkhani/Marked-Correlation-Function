# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 14:48:29 2021

@author: zahra
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def distance_calculator(X,Y): #between a point and an array
    
    distance = ((X[:,0]-Y[0])**2+(X[:,1]-Y[1])**2+(X[:,2]-Y[2])**2)**(0.5)
    
    return distance



#This function calculated the marked correlation function for a non-periodic box
def MCF(data,marks,rbins,random_catalog):#data = XYZ position of galaxies
    
    
    r_max = np.amax(rbins)
    
    N_random =len(random_catalog)
     
    #creating a random catalogue
    
    distance = []
    distance_random = []
    distance_DR = []
    DD = []
    RR = []
    DR= []
    mark_array=[]
    mark_array_WR = []


    WW = []
    WR = []
    
    #To make the code faster, I considered a box with the size of 2*r_max around each galaxy.
    # I just computed separarions for galaxies in this box
    for i in range(len(data)):
        
        
        #conditions = True-False arrays(Nx1) showing that which galaxies have X/Y/Z component in the range
        condition_x = (data[:,0]<data[i,0]+r_max)&(data[:,0]>data[i,0]-r_max)
        condition_y= (data[:,1]<data[i,1]+r_max)&(data[:,1]>data[i,1]-r_max)
        condition_z= (data[:,2]<data[i,2]+r_max)&(data[:,2]>data[i,2]-r_max)
    
        #In order to avoid counting a specific galaxy pair twice, I began with i+1th galaxy
        distance.append(distance_calculator(data[i+1:,:][(condition_x&condition_y&condition_z)[i+1:]],data[i]))
        
        mark_array.append(marks[i]*marks[i+1:][(condition_x&condition_y&condition_z)[i+1:]]) #m_i*m_j/m^2
    print('DD and WW')
        
    # the same procedure for the random catalog    
    for i in range(len(random_catalog)):
        
        condition_x = (random_catalog[:,0]<random_catalog[i,0]+r_max)&(random_catalog[:,0]>random_catalog[i,0]-r_max)
        condition_y= (random_catalog[:,1]<random_catalog[i,1]+r_max)&(random_catalog[:,1]>random_catalog[i,1]-r_max)
        condition_z= (random_catalog[:,2]<random_catalog[i,2]+r_max)&(random_catalog[:,2]>random_catalog[i,2]-r_max)
        
        distance_random.append(distance_calculator(random_catalog[i+1:,:][(condition_x&condition_y&condition_z)[i+1:]],random_catalog[i]))
    #flattening the list ==> distance =[[1,2,3],[4,5],[6]] ==> distance = [1,2,3,4,5,6]
    print('RR')
    for i in range(len(data)): 
        
        #conditions = True-False arrays(Nx1) showing that which galaxies have X/Y/Z component in the range
        condition_x = (random_catalog[:,0]<data[i,0]+r_max)&(random_catalog[:,0]>data[i,0]-r_max)
        condition_y= (random_catalog[:,1]<data[i,1]+r_max)&(random_catalog[:,1]>data[i,1]-r_max)
        condition_z= (random_catalog[:,2]<data[i,2]+r_max)&(random_catalog[:,2]>data[i,2]-r_max)
    
        #In order to avoid counting a specific galaxy pair twice, I began with i+1th galaxy
        distance_DR.append(distance_calculator(random_catalog[(condition_x&condition_y&condition_z)],data[i]))
        marks_WR = [marks[i]]*len(distance_DR[i])
        mark_array_WR.append(marks_WR)
    print('DR and WR')
        
    distance = [i for sublist in distance for i in sublist] 
    distance_random = [i for sublist in distance_random for i in sublist]
    distance_DR = [i for sublist in distance_DR for i in sublist]
    
    distance = np.array(distance)
    distance_random = np.array(distance_random)
    distance_DR = np.array(distance_DR)
    
    mark_array = [i for sublist in mark_array for i in sublist]
    mark_array_WR = [i for sublist in mark_array_WR for i in sublist]
    
    
    mark_array = np.array(mark_array)
    mark_array_WR=np.array(mark_array_WR)
    #computing the histogram
    for j in range(len(rbins)-1):
        condition = (distance<rbins[j+1])&(distance>rbins[j])
        DD.append(len(distance[condition]))
        
        WW.append(mark_array[condition].sum())
        
        condition = (distance_random<rbins[j+1])&(distance_random>rbins[j])
        RR.append(len(distance_random[condition]))
        
        condition = (distance_DR<rbins[j+1])&(distance_DR>rbins[j])
        DR.append(len(distance_DR[condition]))
        WR.append(mark_array_WR[condition].sum())
        
        

      #xi = DD(r)/(N(N-1)/2)/RR(r)/(Nr(Nr-1)/2) -1  
    xi = ((np.array(DD)/(len(data)*(len(data)-1)/2)) - (2*np.array(DR)/(len(data)*N_random)) +(np.array(RR)/(N_random*(N_random-1)/2)))/(np.array(RR)/(N_random*(N_random-1)/2)) 
    
 

    

        
  
        
   

    
    

        
        
        
        
  

        
    W = ((np.array(WW)/(len(data)*(len(data)-1)/2)) - (2*np.array(WR)/(len(data)*N_random)) +(np.array(RR)/(N_random*(N_random-1)/2)))/(np.array(RR)/(N_random*(N_random-1)/2)) 
    
    M = (1+W)/(1+xi)
        
    
    return(M)
#This function is used to compute Jackknife errorbar        
def jackknife(data,marks,rbins,n_random,L_box,n_jackknife):
    M_subsample = [] 
    C_ij = np.zeros((len(rbins)-1,len(rbins)-1))
    Error = []
    N_random = (len(data))*n_random
    random= np.random.uniform(0,L_box,(N_random,3))
    #first we should make the subsamples
    #for each subsample, we should compute CF using the above function
    #then we should compute the covariance matrix and errorbar
    #return everything
    for i in range(n_jackknife):
        for j in range(n_jackknife):
            for k in range(n_jackknife):
                condition_x = (data[:,0]<((i+1)*L_box/n_jackknife))&(data[:,0]>(i*L_box/n_jackknife))
                condition_y = (data[:,1]<((j+1)*L_box/n_jackknife))&(data[:,1]>(j*L_box/n_jackknife))
                condition_z = (data[:,2]<((k+1)*L_box/n_jackknife))&(data[:,2]>(k*L_box/n_jackknife))
                sub_sample = data[((condition_x&condition_y&condition_z)==False)]
                mark_sample = marks[((condition_x&condition_y&condition_z)==False)]
                
                condition_x = (random[:,0]<((i+1)*L_box/n_jackknife))&(random[:,0]>(i*L_box/n_jackknife))
                condition_y = (random[:,1]<((j+1)*L_box/n_jackknife))&(random[:,1]>(j*L_box/n_jackknife))
                condition_z = (random[:,2]<((k+1)*L_box/n_jackknife))&(random[:,2]>(k*L_box/n_jackknife))
                sub_sample_random = random[((condition_x&condition_y&condition_z)==False)]
                M_subsample.append(MCF(sub_sample,mark_sample,rbins=rbins,random_catalog=sub_sample_random))
                print('n=',k)
                
    M_subsample = np.array(M_subsample)
    #Avg_CF = np.mean(CF_subsample,axis = 0)
    n_sample = n_jackknife**3  
    for i in range(len(rbins)-1):
        for j in range(len(rbins)-1):
            C_ij[i,j] = ((n_sample-1)/n_sample)*np.sum((M_subsample[:,i]-np.mean(M_subsample[:,i]))*(M_subsample[:,j]-np.mean(M_subsample[:,j])))
            if (i == j):
                Error.append(np.sqrt(C_ij[i,j]))
    Error = np.array(Error)
    M_mean = np.mean(M_subsample,axis=0)
    
    return M_mean , C_ij, Error, M_subsample    
    