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
        
        lon1, lat1 = myProj(x1, y1, inverse = True)
        lon2, lat2 = myProj(x2, y2, inverse = True)
        
        LAT = [lat1, lat2]
        LONG = [lon1, lon2]
       
        ax.plot(LONG,LAT, c=col)

def findindlines(hist, x,y, xedges, yedges, thres):
    
    # finds histogram of the data 
    # then mask it based on threshold of counts 
    # returns a matrix of number of bins above threshold x elements in x 
    # and returns the counts of each of those clusters 
    
    # hist - arraylike 
    #       2d histogram of data
    
    # xedges, yedges - arraylike 
    #                 edges of bins for x and y (can be uneven/different)
    
    # thres - int   
    #       threshold of counts in bins you want to return
    
    ## Returns 
    
    # cluster_dat - array like ncluster by 6 
    # 0 cluster number 1 count in that bin 
    # 2 mean theta value 3 mean rho value 
    # 4 std theta   5 std rho 
    
    # ind - matrix of len(x) by nclusters type boolean
    #       boolean matrix of all the indices in each cluster 
    
    # hclusters - matrix like h.T 
    #       matrix with cluster labels in their respective bins 
    #       use to pcolormesh to see colors for each cluster 
    #       relative to the hist2d
    
    
    h = hist.T  # np.histogram2d transposes x and y, therefore, transpose the resulting array
    desired = h > thres
    #plt.pcolormesh(xedges, yedges, desired, cmap='coolwarm', ec='white', lw=2)
    hclusters=np.zeros_like(h)
    mask = np.zeros_like(x, dtype=np.bool)  # start with mask all False
    
    ncluster=np.sum(np.sum(desired))
    cluster_dat=np.zeros([ncluster,6])
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
                cluster_dat[ccount][2]=np.average(x[mask])
                cluster_dat[ccount][3]=np.average(y[mask])
                cluster_dat[ccount][4]=np.std(x[mask])
                cluster_dat[ccount][5]=np.std(y[mask])
                
                hclusters[j,i]=ccount
                ccount+=1
                
                
    return cluster_dat, ind, hclusters

def pttoplotline(p,t, center, length):

    m=1/(np.tan( np.deg2rad(t)))
    
    # put p back into meters
    p=p*1000
    
    b=p*(1/np.sin(np.deg2rad(t)))

    # the intercept is the point 
    # (0,b)
    # in the reference frame it is 
    # (center[0], b+center[1])

    blon,blat=myProj(center[0], center[1]+b, inverse=True)
  
    xs=np.array([blon-length/2, blon+length/2])
    y=-m*(xs-blon)+blat

    return xs, y
  
    
def detectlines(dikeset, areamask, tei, pei, rthres, pleaseplot=True):
    
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
    
    # rthres - %
    #       relative threshold of counts for bin to be considered a line segment
    
    # pleaseplot - boolean (Default True)
    #            flag to plot data 
    
    ## Returns 
    
    # clusters - array like ncluster by 6 
    #       0 cluster number 1 count in that bin 
    #       2 mean theta value 3 mean rho value 
    #       4 std theta   5 std rho 
    
    # lineind - matrix of len(x) by nclusters type boolean
    #       boolean matrix of all the indices in each cluster 
    #       (See findindlines)
    
    X=dikeset[areamask]
    theta=pd.to_numeric(X['theta'], errors='coerce')
    rho=pd.to_numeric(X['p'], errors='coerce')
    

    center=[dikeset['Xstart'].values.mean(), dikeset['Ystart'].values.mean()]
    centerlatlon= myProj(center[0],center[1], inverse = True)
    
    hist, _xedges, _yedges = np.histogram2d(theta, rho, bins=[tei, pei])
    thres=rthres*np.max(hist)
    clusters,lineind, hclusters= findindlines(hist, theta, rho, tei, pei, thres)
    length=(np.abs(lon1-lon2))*0.8
    
    if pleaseplot:
        
        fig, ax=plt.subplots(1,3)
        fig.set_size_inches(11,8.5)
        h,te,pe,im=ax[2].hist2d(theta,rho, bins=[tei, pei], cmap=plt.cm.BuPu)
        #plotlines(dikeset, areamask,'k', ax[0])
        plotlines(dikeset, areamask,'k', ax[1])
        plt.colorbar(im, label='counts, threshold='+str(thres))
        for i in range(0,len(clusters)):
            colors = plt.cm.rainbow(np.linspace(0, 1, len(clusters)))
            mask=lineind[i]
            plotlines(X,mask,colors[i],ax[0])
            
            t=clusters[i,2]
            p=clusters[i,3]
            
            xs,y=pttoplotline(p, t, center, length)
            ax[1].plot(xs,y, c=colors[i])
                
            
        tol=0.05
        ax[1].set_ylim([lat1-tol,lat2+tol])
        ax[1].set_xlim([lon1-tol,lon2+tol])

        ax[0].set_ylim([lat1-tol,lat2+tol])
        ax[0].set_xlim([lon1-tol,lon2+tol])
        plt.show()
        
        
    return clusters, lineind

dikeset=pd.read_csv("/home/akh/myprojects/dikes/dikeset_ptheta.csv")

myProj = Proj("+proj=utm +zone=11, +ellps=WGS84 +datum=WGS84 +units=m")


#dikeset=pd.read_csv("PATH/TO/FILE/")
lat1=45.01
lat2=45.5

lon1=-117.5
lon2=-117.1

# center for HT
center=[dikeset['Xstart'].values.mean(), dikeset['Ystart'].values.mean()]
centerlatlon= myProj(center[0],center[1], inverse = True)
xs=np.linspace(lon1, lon2, 10)

## set mask to wallowas only

masklat= (dikeset['latitude'] > lat1) & (dikeset['latitude'] < lat2)
masklong=(dikeset['longitude'] > lon1) & (dikeset['longitude'] < lon2)
masklatlong= (masklat==1) & (masklong==1)

X=dikeset[masklatlong]
theta=pd.to_numeric(X['theta'], errors='coerce')
rho=pd.to_numeric(X['p'], errors='coerce')
    
pmin=min(rho)
pmax=max(rho)
pmean=np.average(rho)
tmean=np.average(theta)


plength=1
ttol=5
rthres=0.8

## Try the following
# 5,5,0.8
# 10,5,0.8




pei_1=np.arange(pmin, pmax, plength)
tei_1=np.arange(min(theta), max(theta), ttol)
pstd=np.std(rho)
tstd=np.std(theta)

x=np.linspace(-3,3,10)
pei_2=pmean+x*pstd
tei_2=tmean+x*tstd

cluster, ind= detectlines(dikeset, masklatlong, tei_1, pei_1, rthres)


pei_2=np.hough=X.drop(columns=['SDMBearing', 'seg_length','host_rock', 'Id', 'SDMAzimuth', 'Xstart', 'Xend', 'Yend', 'Ystart', 'longitude'])
