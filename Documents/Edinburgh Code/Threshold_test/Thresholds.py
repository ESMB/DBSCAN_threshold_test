#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 18:30:37 2021

@author: Mathew
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from sklearn.cluster import DBSCAN


     

# Settings
image_width=512
image_height=512
Pixel_size=103
scale=8

# Thresholds to try
# This is the distance threshold (pixels)
eps=[]          


eps.append(1)


# This is the minimum number of locs per cluster

minimum_locs=[]

minimum_locs.append(50)



filename_contains="FitResults"



# Files to analyse

pathList=[]

pathList.append(r"/Users/Mathew/Documents/Current analysis/Actin LIVE-PAINT/")
        
# Various functions

#  Generate SR image (points)
def generate_SR(coords):
    SR_plot_def=np.zeros((image_width*scale,image_height*scale),dtype=float)
    j=0
    for i in coords:
        
        xcoord=coords[j,0]
        ycoord=coords[j,1]
        scale_xcoord=round(xcoord*scale)
        scale_ycoord=round(ycoord*scale)
        # if(scale_xcoord<image_width and scale_ycoord<image_height):
        SR_plot_def[scale_ycoord,scale_xcoord]+=1
        
        j+=1
    return SR_plot_def

# This is to generate a gaussian.
def gkern(l,sigx,sigy):

    ax = np.linspace(-(l - 1) / 2., (l - 1) / 2., l)
    xx, yy = np.meshgrid(ax, ax)

    kernel = np.exp(-0.5 * (np.square(xx)/np.square(sigx) + np.square(yy)/np.square(sigy)) )

    return kernel/np.sum(kernel)

# Generate image with width = precision
def generate_SR_prec(coords,precsx,precsy):
    box_size=20
    SR_prec_plot_def=np.zeros((image_width*scale,image_height*scale),dtype=float)
    dims=np.shape(SR_prec_plot_def)
    print(dims)
    j=0
    for i in coords:

      
        precisionx=precsx[j]/Pixel_size*scale
        precisiony=precsy[j]/Pixel_size*scale
        xcoord=coords[j,0]
        ycoord=coords[j,1]
        scale_xcoord=round(xcoord*scale)
        scale_ycoord=round(ycoord*scale)
        
        
        
        sigmax=precisionx
        sigmay=precisiony
        
        
        # tempgauss=SRGaussian((2*box_size,2*box_size), (sigmax,sigmay),(box_size,box_size))
        
        tempgauss=gkern(2*box_size,sigmax,sigmay)
        
        # SRGaussian((2*box_size,2*box_size), (sigmax,sigmay),(box_size,box_size))
        
        
        
        ybox_min=scale_ycoord-box_size
        ybox_max=scale_ycoord+box_size
        xbox_min=scale_xcoord-box_size
        xbox_max=scale_xcoord+box_size 
        
        
        if(np.shape(SR_prec_plot_def[ybox_min:ybox_max,xbox_min:xbox_max])==np.shape(tempgauss)):
            SR_prec_plot_def[ybox_min:ybox_max,xbox_min:xbox_max]=SR_prec_plot_def[ybox_min:ybox_max,xbox_min:xbox_max]+tempgauss
        
        
           
        j+=1
    
    return SR_prec_plot_def

# Perform DBSCAN on the coordinates. 

def cluster(coords,eps_threshold,minimum_locs_threshold):
     db = DBSCAN(eps=eps_threshold, min_samples=minimum_locs_threshold).fit(coords)
     labels = db.labels_
     n_clusters_ = len(set(labels)) - (1 if-1 in labels else 0)  # This is to calculate the number of clusters.
     print('Estimated number of clusters: %d' % n_clusters_)
     return labels

def generate_SR_cluster(coords,clusters):
    SR_plot_def=np.zeros((image_width*scale,image_height*scale),dtype=float)
    j=0
    for i in coords:
        if clusters[j]>-1:
            xcoord=coords[j,0]
            ycoord=coords[j,1]
            scale_xcoord=round(xcoord*scale)
            scale_ycoord=round(ycoord*scale)
            # if(scale_xcoord<image_width and scale_ycoord<image_height):
            SR_plot_def[scale_ycoord,scale_xcoord]+=1
        
        j+=1
    return SR_plot_def

def generate_SR_prec_cluster(coords,precsx,precsy,clusters):
    box_size=50
    SR_prec_plot_def=np.zeros((image_width*scale+100,image_height*scale+100),dtype=float)
    SR_fwhm_plot_def=np.zeros((image_width*scale+100,image_height*scale+100),dtype=float)

    j=0
    for clu in clusters:
        if clu>-1:
       
            precisionx=precsx[j]/Pixel_size*scale
            precisiony=precsy[j]/Pixel_size*scale
            xcoord=coords[j,0]
            ycoord=coords[j,1]
            scale_xcoord=round(xcoord*scale)+50
            scale_ycoord=round(ycoord*scale)+50
            
            sigmax=precisionx
            sigmay=precisiony
            
            
            # tempgauss=SRGaussian((2*box_size,2*box_size), (sigmax,sigmay),(box_size,box_size))
            tempgauss=gkern(2*box_size,sigmax,sigmay)
            ybox_min=scale_ycoord-box_size
            ybox_max=scale_ycoord+box_size
            xbox_min=scale_xcoord-box_size
            xbox_max=scale_xcoord+box_size 
        
        
            if(np.shape(SR_prec_plot_def[ybox_min:ybox_max,xbox_min:xbox_max])==np.shape(tempgauss)):
                SR_prec_plot_def[ybox_min:ybox_max,xbox_min:xbox_max]=SR_prec_plot_def[ybox_min:ybox_max,xbox_min:xbox_max]+tempgauss
                
            tempfwhm_max=tempgauss.max()
            tempfwhm=tempgauss>(0.5*tempfwhm_max)
            
            tempfwhm_num=tempfwhm*(clu+1)
           
            
            if(np.shape(SR_fwhm_plot_def[ybox_min:ybox_max,xbox_min:xbox_max])==np.shape(tempfwhm)):
               plot_temp=np.zeros((2*box_size,2*box_size),dtype=float)
               plot_add=np.zeros((2*box_size,2*box_size),dtype=float)
               plot_temp=SR_fwhm_plot_def[ybox_min:ybox_max,xbox_min:xbox_max]
               plot_add_to=plot_temp==0
               
               plot_add1=plot_temp+tempfwhm_num
               
               plot_add=plot_add1*plot_add_to
               
               SR_fwhm_plot_def[ybox_min:ybox_max,xbox_min:xbox_max]=SR_fwhm_plot_def[ybox_min:ybox_max,xbox_min:xbox_max]+plot_add
                
                
                # (SR_fwhm_plot_def[scale_ycoord-box_size:scale_ycoord+box_size,scale_xcoord-box_size:scale_xcoord+box_size]+tempfwhm_num).where(SR_fwhm_plot_def[scale_ycoord-box_size:scale_ycoord+box_size,scale_xcoord-box_size:scale_xcoord+box_size]==0)
                # SR_tot_plot_def[scale_ycoord-box_size:scale_ycoord+box_size,scale_xcoord-box_size:scale_xcoord+box_size]=SR_tot_plot_def[scale_ycoord-box_size:scale_ycoord+box_size,scale_xcoord-box_size:scale_xcoord+box_size]+tempfwhm
            
            # SR_tot_plot_def[SR_tot_plot_def==0]=1
            labelled=SR_fwhm_plot_def
            
            SR_prec_plot=SR_prec_plot_def[50:image_width*scale+50,50:image_height*scale+50]
            labelled=labelled[50:image_width*scale+50,50:image_height*scale+50]
            
            
        j+=1
    
    return SR_prec_plot,labelled,SR_fwhm_plot_def
   


for path in pathList:
    print(path)

    # Load the fits:
    for root, dirs, files in os.walk(path):
                for name in files:
                        if filename_contains in name:
                            if ".txt" in name:
                                # if "_FitResults" not in name:
                                    resultsname = name
                                    print(resultsname)
    
    fits_path=path+resultsname
    
    loc_data = pd.read_table(fits_path,sep=",")
    

    # Extract useful data:
    coords = np.array(list(zip(loc_data['X'],loc_data['Y'])))
    precsx= np.array(loc_data['Precision (nm)'])
    precsy= np.array(loc_data['Precision (nm)'])
    xcoords=np.array(loc_data['X'])
    ycoords=np.array(loc_data['Y'])
    
    
    precs_nm=precsx
    
    plt.hist(precs_nm, bins = 50,range=[0,100], rwidth=0.9,color='#ff0000')
    plt.xlabel('Precision (nm)',size=20)
    plt.ylabel('Number of Features',size=20)
    plt.title('Localisation precision',size=20)
    plt.savefig(path+"Precision.pdf")
    plt.show()
        
    # Generate points SR (ESMB method):
    SR=generate_SR(coords)
    
    imsr = Image.fromarray(SR)
    imsr.save(path+filename_contains+'SR_points_python.tif')
    
    SR_prec=generate_SR_prec(coords,precsx,precsy)
    
    imsr = Image.fromarray(SR_prec)
    imsr.save(path+'SR_width_python.tif')
    

# Cluster analysis:
    
    for e in eps:
        for locs in minimum_locs:
            clusters=cluster(coords,e,locs)
            
            path2save=path+"e_"+str(e)+"_locs_"+str(locs)+""
            # Check how many localisations per cluster
         
            cluster_list=clusters.tolist()    # Need to convert the dataframe into a list- so that we can use the count() function. 
            maximum=max(cluster_list)+1  
            
            
            cluster_contents=[]         # Make a list to store the number of clusters in
            
            for i in range(0,maximum):
                n=cluster_list.count(i)     # Count the number of times that the cluster number i is observed
               
                cluster_contents.append(n)  # Add to the list. 
            
            if len(cluster_contents)>0:
                average_locs=sum(cluster_contents)/len(cluster_contents)
         
                plt.hist(cluster_contents, bins = 20,range=[0,200], rwidth=0.9,color='#607c8e') # Plot a histogram. 
                plt.xlabel('Localisations per cluster')
                plt.ylabel('Number of clusters')
                plt.savefig(path2save+'Localisations.pdf')
                plt.show()
                
                cluster_arr=np.array(cluster_contents)
            
                median_locs=np.median(cluster_arr)
                mean_locs=cluster_arr.mean()
                std_locs=cluster_arr.std()
                
            
                # Generate the SR image.
                SR_Clu=generate_SR_cluster(coords,clusters)
                
                imsr = Image.fromarray(SR_Clu)
                imsr.save(path2save+'SR_points_python_clustered.tif')
                
                SR_clu_prec,labelled,SR_prec_plot=generate_SR_prec_cluster(coords,precsx,precsy,clusters)
                
                

            
                imsr = Image.fromarray(SR_clu_prec)
                imsr.save(path2save+'SR_width_python_clustered.tif')
                
                imsr = Image.fromarray(labelled)
                imsr.save(path2save+'SR_fwhm_python_clustered.tif')
                
                labeltot=labelled.max()+1
            
                print('Total number of clusters in labelled image: %d'%labeltot)
                
                
                
            
            
