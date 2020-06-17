import numpy as np 
import math 
import pandas as pd
from scipy import * 

#from skimage.transform import probabilistic_hough_line as ProbHough
#from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import seaborn as sns 
from pyproj import Proj

sns.set()
np.random.seed(5)

print('starting')
def plotlines(data, maskar, col, ax):
    
    #plots the line segments contained in data[maskar]
    # in the color col 
    # converts UTM to lat/long
    
    # data - data frame with columns 
    # 'Xstart', 'Xend', 'Yend', 'Ystart'
    # which are UTM coordinates 
    
    # maskar - arraylike, type logical mask of len(data)
    #       masking of data you want to plot
    
    # col   - string or RGB tuple 
    #        color you want all lines to be plotted
    
    # ax    - object you want to plot on 
    
    temp=data.loc[maskar]
    
    for i in range(0,len(temp)):
        x1=temp['Xstart'].iloc[i]
        y1=temp['Ystart'].iloc[i]
        y2=temp['Yend'].iloc[i]
        x2=temp['Xend'].iloc[i]

        myProj = Proj("+proj=utm +zone=11, +ellps=WGS84 +datum=WGS84 +units=m")
        
        lon1, lat1 = myProj(x1, y1, inverse = True)
        lon2, lat2 = myProj(x2, y2, inverse = True)
        
        LAT = [lat1, lat2]
        LONG = [lon1, lon2]
       
        ax.plot(LONG,LAT, c=col)

def findindlines(x,y, xedges, yedges, thres):
    
    # finds histogram of the data 
    # then mask it based on threshold of counts 
    # returns a matrix of number of bins above threshold x elements in x 
    # and returns the counts of each of those clusters 
    
    # x, y - arraylike 
    #       data you wish to bin 
    
    # xedges, yedges - arraylike 
    #                 edges of bins for x and y (can be uneven/different)
    
    # thres - int   
    #       threshold of counts in bins you want to return
    
    hist, _xedges, _yedges = np.histogram2d(x, y, bins=[xedges, yedges])
    
    h = hist.T  # np.histogram2d transposes x and y, therefore, transpose the resulting array
    desired = h > thres
    #plt.pcolormesh(xedges, yedges, desired, cmap='coolwarm', ec='white', lw=2)
    hclusters=np.zeros_like(h)
    mask = np.zeros_like(x, dtype=np.bool)  # start with mask all False
    
    ncluster=np.sum(np.sum(desired))
    cluster_dat=np.zeros([ncluster,2])
    ind=np.matlib.repmat(mask, ncluster, 1)
    ccount=0
    
    for i in range(len(xedges) - 1):
        for j in range(len(yedges) - 1):
            if desired[j, i]:
                mask = np.zeros_like(x, dtype=np.bool)
                # print(f'x from {xedges[i]} to {xedges[i + 1]} y from {yedges[j]} to {yedges[j + 1]}')
                mask = np.logical_or(mask, (x >= xedges[i]) & (x < xedges[i + 1]) & (y >= yedges[j]) & (y < yedges[j + 1]))
                ind[ccount]=mask
                cluster_dat[ccount][0]=ccount
                cluster_dat[ccount][1]= h[j,i]
                hclusters[j,i]=ccount
                ccount+=1
                
                
    return cluster_dat, ind, hclusters


    
def detectlines(dikeset, areamask, plength, ttol,thres, pleaseplot=True):
    
    # takes dikeset and uses the calculated 
    # accumulator array (hist2d) to connect line segments 
    # plots this and accumulator array 
    
    # dikeset - dataframe of dike data
    
    # areamask - arraylike of type logical with len(dikeset)
    #       chooses dikes to look at within data set 
    
    # pei - arraylike 
    #       edges of bins for p 
    
    # tei - arraylike 
    #       edges of bins for theta
    
    # thres - int 
    #       threshold of counts for bin to be considered a line segment
    
    # pleaseplot - logical (Default True)
    
    
    X=dikeset[areamask]
    hough=X.drop(columns=['SDMBearing', 'seg_length','host_rock', 'Id', 'SDMAzimuth', 'Xstart', 'Xend', 'Yend', 'Ystart', 'longitude'])
    hough['latitude']=pd.to_numeric(hough['latitude'], errors='coerce')
    theta=pd.to_numeric(hough['theta'], errors='coerce')
    rho=pd.to_numeric(hough['p'], errors='coerce')
    

    

    
    clusters,lineind, hclusters= findindlines(theta, rho, tei, pei, thres)
    
    
    if pleaseplot:
        
        fig, ax=plt.subplots(1,2)
        fig.set_size_inches(11,8.5)
        h,te,pe,im=ax[1].hist2d(theta,rho, bins=[tei, pei], cmap=plt.cm.BuPu)
        plotlines(dikeset, areamask,'k', ax[0])
        c=plt.colorbar(im, label='counts, threshold='+str(thres))
        for i in range(0,len(clusters)):
            colors = plt.cm.rainbow(np.linspace(0, 1, len(clusters)))
            mask=lineind[i]
            plotlines(X,mask,colors[i],ax[0])
            
        #ax[2].pcolormesh(hclusters, cmap=plt.cm.rainbow)
        plt.show()
        
        
    return clusters, lineind

#dikeset=pd.read_csv("/home/akh/myprojects/dikes/dikeset_ptheta.csv")
dikeset=pd.read_csv("PATH/TO/FILE/")
lat1=45.01
lat2=45.5

lon1=-117.5
lon2=-117.1


masklat= (dikeset['latitude'] > lat1) & (dikeset['latitude'] < lat2)
masklong=(dikeset['longitude'] > lon1) & (dikeset['longitude'] < lon2)
masklatlong= (masklat==1) & (masklong==1)

X=dikeset[masklatlong]
hough=X.drop(columns=['SDMBearing', 'seg_length','host_rock', 'Id', 'SDMAzimuth', 'Xstart', 'Xend', 'Yend', 'Ystart', 'longitude'])
hough['latitude']=pd.to_numeric(hough['latitude'], errors='coerce')
theta=pd.to_numeric(hough['theta'], errors='coerce')
rho=pd.to_numeric(hough['p'], errors='coerce')
    
pmin=min(rho)
pmax=max(rho)

plength=2
ttol=10
thres=30
pei=np.arange(pmin, pmax, plength)

tei=np.arange(min(theta), max(theta), ttol)

cluster, ind= detectlines(dikeset, masklatlong, plength, ttol,thres, pleaseplot=True)