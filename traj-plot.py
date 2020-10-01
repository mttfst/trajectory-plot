#!/usr/bin/env python
# coding: utf-8

# # Traj Plot - Post-processing and visualization of trajectory data.

# ## Used packages

# In[1]:


from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap#, cm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import PercentFormatter
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
from matplotlib.patches import Polygon
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
import time
import datetime
import os
import sys

#import warnings
#warnings.filterwarnings("ignore")


# ## Functions
# 
# ### rot2geo 
# Transformation of rotated coordinates to geographic coordinates.
# 
# $\lambda_{geo}$... geographic longitude 
# 
# $\varphi_{geo}$... geographic latitude 
# 
# $\lambda_{rot}$... rotated longitude 
# 
# $\varphi_{rot}$... rotated latitude 
# 
# $\lambda_{rot}^{pol}$... longitude of the rotated north pole
# 
# $\varphi_{rot}^{pol}$... latitude of the rotated north pole
# 
# $\varphi_{geo} = \dfrac{180^{\circ}}{\pi} \arcsin \left[ \sin\varphi_{rot}\,\sin\varphi_{rot}^{pol} + \cos\varphi_{rot}\,\cos\lambda_{rot}\,\cos\varphi_{rot}^{pol} \right]   $
# 
# $\lambda_{geo} = \arctan \left[ \dfrac{ \cos\varphi_{rot}\,\sin\lambda_{rot} }{ \sin\varphi_{rot}^{pol}\,\cos\varphi_{rot}\,\cos\lambda_{rot}\,  - \sin\varphi_{rot}\,\cos\varphi_{rot}^{pol}}\right] + \lambda_{rot}^{pol} + 180^{\circ} $

# In[2]:


def cm2inch(value):
    return value/2.54

def rot2geo(xin,yin):


    pollon=-170.
    pollat=40.
    xgeo=[i for i in xin]
    ygeo=[i for i in yin]

    for i in range(len(xin)):
        ygeo[i]=180./np.pi * np.arcsin(np.sin(np.pi/180.*yin[i])*np.sin(np.pi/180.*pollat) + np.cos(np.pi/180.*yin[i])*np.cos(np.pi/180.*xin[i])*np.cos(np.pi/180.*pollat))

        xgeo[i]=180./np.pi * np.arctan((np.cos(np.pi/180.*yin[i])*np.sin(np.pi/180.*xin[i]))/(np.sin(np.pi/180.*pollat)*np.cos(np.pi/180.*yin[i])*np.cos(np.pi/180.*xin[i])-np.sin(np.pi/180.*yin[i])*np.cos(np.pi/180.*pollat))) + pollon + 180.

    return(xgeo,ygeo)

def meter2degree(meter):
    degree=np.zeros(len(meter))
    for i in range(len(meter)):
        degree[i] = meter[i]/6.3710088e6 * 180./np.pi
    return degree

def degree2meter(degree):
    meter=np.zeros(len(degree))
    for i in range(len(degree)):
        meter[i] = degree[i]*6.3710088e6 * np.pi/180.
    return meter


# ## Data input
# 
# Read trajectory data from NetCDF.
# Mostly the trajectory data is split into multiple files.
# You need to specify with the variables first_file and last_file which files are read.
# 
# All trajectories from all files are collected into the np.arrays: <br>
# z,zrel,rlon,rlat

# In[3]:


print()
print('Start')
print('│')
print('├──read data')

# number of first and last file
first_file = 0 
last_file = 2

# find files
stime=time.time()
files = []
for (dirpath, dirnames, filenames) in os.walk("."):
    files.extend(filenames)


files=sorted(files)
namelist=[]
# datestring=[]
for name in files:
    if name[-3:]=='.nc' and name[:4]=='part':
        namelist.append(name)
        # datestring.append(name[5:-3])

# quick look into the files to get the number 
# and the length of trajectories

ntraj = 0
maxlen=0
for fn in namelist[first_file:last_file]:
    
    f=Dataset(fn,'r')
    
    if ntraj == 0:
        times = f.variables['time'][:]
        
    nid = f.dimensions['id']
    ntraj += nid.size
    
    ntime = f.dimensions['time']
    maxlen = max(maxlen,ntime.size)
    
    #print(nid.size,ntime.size)
    indt      = f.getncattr('output_time_step_in_sec')  
    ref_year  = f.getncattr('ref_year')
    ref_month = f.getncattr('ref_month')
    ref_day   = f.getncattr('ref_day')
    ref_h     = f.getncattr('ref_hour')
    ref_min   = f.getncattr('ref_min')
    ref_s     = f.getncattr('ref_sec')
    
    
    
    f.close()


z    = np.zeros((ntraj,maxlen))
zrel = np.zeros((ntraj,maxlen))
rlat = np.zeros((ntraj,maxlen))
rlon = np.zeros((ntraj,maxlen))
size = np.zeros(ntraj)

tid_min=0
tid_max=0

for fn in namelist[first_file:last_file]:

    sys.stdout.write('\r'+'│  │ '+fn)
    sys.stdout.flush()
    
    f=Dataset(fn,'r')
    
    inz = np.asarray(f.variables['z'][:])
    inz = np.transpose(inz)
    
    inzrel = np.asarray(f.variables['zrel'][:])
    inzrel = np.transpose(inzrel)
    
    inlat  = np.asarray(f.variables['latitude'][:])
    inlat = np.transpose(inlat)
    
    inlon  = np.asarray(f.variables['longitude'][:])
    inlon = np.transpose(inlon)
    
    insize = np.asarray(f.variables['size'][:])
    
    tid_max+=inz.shape[0]
    tlen=inz.shape[1]
    
    
    z[tid_min:tid_max,-tlen:]=inz
    z[tid_min:tid_max,:-tlen]=np.asarray([[inz[i,1]for  j in range(maxlen-tlen)]for i in range(inz.shape[0])])
    z[tid_min:tid_max,0]=inz[:,0]
    z[tid_min:tid_max,-tlen]=inz[:,1]
    
    zrel[tid_min:tid_max,-tlen:]=inzrel
    zrel[tid_min:tid_max,:-tlen]=np.asarray([[inzrel[i,0]for  j in range(maxlen-tlen)]for i in range(inz.shape[0])])
    
    rlat[tid_min:tid_max,-tlen:]=inlat
    rlat[tid_min:tid_max,:-tlen]=np.asarray([[inlat[i,0]for  j in range(maxlen-tlen)]for i in range(inz.shape[0])])
    
    rlon[tid_min:tid_max,-tlen:]=inlon
    rlon[tid_min:tid_max,:-tlen]=np.asarray([[inlon[i,0]for  j in range(maxlen-tlen)]for i in range(inz.shape[0])])
    
    size[tid_min:tid_max]=insize
    
    tid_min+=inz.shape[0]
    
    
    
    
    print(' -> '+str(int(time.time()-stime))+'s')

print('│  done '+str(int(time.time()-stime))+'s')
print('│')



#x[:y.shape[0],-y.shape[1]:]=y
#x[:y.shape[0],:-len(y)]=404
#
#x[y.shape[0]:y.shape[0]+z.shape[0],-z.shape[1]:]=z
#x[y.shape[0]:y.shape[0]+z.shape[0],:-len(z)]=404


# ## Remove unnecceary data
# #### and also sort in bins
# 
# In the model output the trajectories look like this:
# 
# ```
# 0 - 0 - 0 - 0 - 0 - 0 - 0 - 1 - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 8 - 8 - 8 - 8 - 8 - 8 
# [ - before start time - | ----- active trajectory ----- | --- dead trajetoy --- ]
# ```
# 
# At the begin of the trajectory, a time span up to one hour is filled with placeholders. 
# When the model reaches the trajectory start time, the trajectory is in its active phase. When the trajectory dies because of deposition or leaving the domain it is frozen on its last position until the end of the simulation. 
# 
# For the plot, we only need the active phase. The rest is removed from the trajectory.
# 
# 
# The first start point is searched by simple iteration through the array of the vertical coordinate (z). The first index (in the array) where z changes its value is marked as the index of the start time. 
# Since usually multiple (thousands) trajectories start at the same time,  for the next trajectory, it is checked first if it starts at the same time as the trajectory before. Otherwise, a new iteration starts. 
# 
# To find the endpoint, I look at the midpoint of the trajectory. If the z value is constant at this point, then the trajectory is already dead there. Then I look at the midpoint of the first half of the trajectory and check if it is dead there as well... and so on until I find the endpoint. 
# 
# Example:
# 
# ```
# 
# 0 - 1 - 2 - 3 - 3 - 3 - 3 - 3 - 3 - 3 - 3 - 3 - 3 - 3 - 3 - 3 - 3 - 3 - 3
#                                 | midpoint dead
#                                 
# 
# 0 - 1 - 2 - 3 - 3 - 3 - 3 - 3 - 3
#                 | midpoint dead
#                 
# 0 - 1 - 2 - 3 - 3
#         | midpoint alive       
# 
# 2 - 3 - 3
#     | midpoint is end point
# 
# ```
# 
# This algorithm is really fast and can proceed hundred thousand of trajectories in just a few seconds.
# 
# 
# #### and also sort in bins
# Since I'm iterating in this cell through all trajectory I use the possibility to sort them into bins of particle size. 
# For other applications, this might be unnecessary.

# In[4]:


print('├──find start and end points')
stime=time.time()

dimz=z.shape

                
# list that holds the pure trajetory information,
# without dead spaces before start and after dead of traj
traid=[]
trastart=[]
trastop=[]

idc=-1
first=True
bintop = [0.28,0.30,0.35,0.4,0.45,0.5,0.58,0.65,0.7,0.8,1.,1.3,1.6,2.0,2.5,3.,3.5,4.,5.,6.5,7.5,8.5,10.,12.5,15.,17.5,20.,25.,30.,32.,50.]
bintopscale=[bintop[i]*1.e-6 for i in range(len(bintop))]
bins = [[] for i in range(len(bintop))]
# start with second traj
for id in range(1,dimz[0],1):
    idc+=1
    traid.append(id)
 
    # identify startpoint of traj 
    if not first:
        pvt=trastart[idc-1] + 1   # previously start time
        
    if not first and z[id,pvt] != z[id,0] and z[id,pvt] != z[id,1] and z[id,pvt-1] == z[id,1]:
        trastart.append(pvt-1)
    else: 
        for t in range(1,dimz[1]):
            #print(t,z[id,t] , z[id,0] , z[id,1])
            if z[id,t] != z[id,0] and z[id,t] != z[id,1]:
                first=False
                trastart.append(t-1)
                break
        
    # last possible case if the identification of the start position above doesn't work
    # this happen to me one time for over 1 Million different trajectories I passed thor this algorithm
    if idc == len(trastart):
        for t in range(1,dimz[1]):
            if zrel[id,t] != zrel[id,0]:
                first=False
                trastart.append(t-1)
                break
        
    # identify position when the particle is dead
    tstop=0
    btime=trastart[idc] # bottom time
    ttime=dimz[1]-1    # top time
    nloop=0
    while ttime-btime > 10: 
        nloop+=1
        jtime=int((btime+ttime)/2)
        if z[id,jtime] != z[id,jtime-1]: # -> alive
            if z[id,jtime] == z[id,jtime+1]: # -> dead
                tstop=jtime
                break
            else: # -> alive -> increase btime
                btime=jtime
        
        else: # -> dead decrease ttime
            ttime=jtime
    
    #print(nloop,btime,ttime,trastop)
    nloop=0
    if tstop==0:
        for t in range(btime,ttime):
            nloop+=1
            if z[id,t] == z[id,t+1]:
                tstop=t
                break
                
    if tstop==0 and ttime == dimz[1]-1:
        tstop=ttime
        
    if trastart[idc]==tstop:
        trastart[idc]-=1
        
    if tstop > 0:
        if tstop == dimz[1]-1:
            trastop.append(tstop)
        else:
            trastop.append(tstop+1)
        
    # sort in bins
    for b in range(len(bins)):
        top = bintop[b] * 1.e-6
        if b == 0:
            bot = 0.
        else:
            bot = bintop[b-1] * 1.e-6

        if size[id] > bot and size[id] <= top:
            bins[b].append(idc)
            break
            
    

                
#print(tra)
#print()
#print(trastart)
print('│  └─done '+str(int(time.time()-stime))+'s')
print('│')


# ## Rotate coordinates

# In[5]:


stime=time.time()

lon,lat= rot2geo (rlon,rlat)
lon=np.asarray(lon)
lat=np.asarray(lat)

print('│  └─done '+str(int(time.time()-stime))+'s')
print('│')


# ## Define the domain for the plot

# In[6]:


# boundaries of the map


# East Brandenbug
lonmin=14.16
lonmax=14.35
latmin=52.44
latmax=52.5
 
# East Germany / West Poland
# lonmin=14.
# lonmax=20.
# latmin=50.5
# latmax=53.


# # East Europe
# lonmin=12.
# lonmax=40.
# latmin=40.
# latmax=54.

# spacing of the plotted lon and lats
dfigx=lonmax - lonmin
if dfigx < 0.5:
    spcnglon=0.05
    adlvl = 0
elif dfigx < 1.0:
    spcnglon=0.2
    adlvl=1
elif dfigx < 5.:
    spcnglon=0.5
    adlvl=2
elif dfigx < 10.:
    adlvl = 2
    spcnglon=2.
else:
    spcnglon=5.
    adlvl = 3

dfigy=latmax - latmin
if dfigy < 0.5:
    spcnglat=0.02
elif dfigy < 1.0:
    spcnglat=0.2
elif dfigy < 5.:
    spcnglat=0.5
elif dfigy < 10.:
    spcnglat=2.
else:
    spcnglat=5.
    
    
# max hight of the vertical plot    
zrelmax=2500#np.max(zrel)


# ## Define the time axis
# If you need to plot a certain time range of the data it is handy to print a list of the dates to find the right indices.

# In[7]:


refdate = datetime.datetime(ref_year,ref_month,ref_day,ref_h,ref_min,ref_s)
reftime = datetime.time(ref_h,ref_min,ref_s)

dates = []
timestemps=[]
datetimes=[]
for t in times:
    time_diff = datetime.timedelta(seconds=int(t))
    date = refdate + time_diff  
    dates.append(mpl.dates.date2num(date))
    datetimes.append(date)
    timestemps.append(date.strftime("%Y%m%d%H%M%S"))

# i=0
# for d in dates:
#    print(i,mpl.dates.num2date(d),timestemps[i])
#    i+=1


# ## Plotting 
# We create a 4x1 panel plot 
# 1. a horizontal map of the trajectories  
# 2. a deposition map 
# 3. vertical dispersion of the trajectories as a time series 
# 4. a color bar

# In[8]:


print('├──create figure')
sbtime=time.time()

# create figure
fig = plt.figure(figsize=(9,12))
gs = gridspec.GridSpec(nrows=4, ncols=1, height_ratios=[1,1,1,0.1])


# create subplots
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[1,0])
ax2 = fig.add_subplot(gs[2,0])
ax3 = fig.add_subplot(gs[3,0])

# set titles
ax0.set_title('a) Horizontal Dispersion',size='x-large')
ax1.set_title('b) Deposition',size='x-large')
ax2.set_title('c) Vertical Dispersion',size='x-large')

print('│  ├─draw map')
mtime=time.time()
m = [None] * 2

# draw maps
m[0] = Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,resolution='i',projection='gall',ax=ax0)
m[1] = Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,resolution='i',projection='gall',ax=ax1)

# draw grid
for map in m:
    map.drawparallels(np.arange(0.,180.,spcnglat),labels=[1,0,0,0],fontsize='x-large')
    map.drawmeridians(np.arange(0.,90.,spcnglon),labels=[0,0,0,1],fontsize='x-large')
    # map.drawcoastlines()
    map.drawcountries()


print('│  │ └─done '+str(int(time.time()-mtime))+'s')
print('│  │ ')

print('│  ├─plot data')
ptime=time.time()


ttot=0.
ncolor=0
nt=0


cmap = LinearSegmentedColormap.from_list(name='mycmap',colors =['cyan','lightskyblue','#9ecae1','#c6dbef','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c', '#f16913','#d94801'],N=31)
ibin=0
stime=time.time()
for bin in bins[:]:
    #print(ibin)
    ibin+=1
    color = cmap(ncolor)
    ncolor+=1
    
    #print(bin)
    trajx=[]
    trajy=[]
    trajz=[]
    ztime=[]
    depox=[]
    depoy=[]
    
    if len(bin)>0:
        for idc in bin:
            # progres bar
            if nt > 0 and np.mod(nt,1000):
                ttot=time.time()-stime
                prog=100*nt/len(traid)
                # print(ttot/nt,lasttraj-nt)
                est=(ttot/nt)*(len(traid)-nt)
                # print((tnow-stime),(lasttraj-nt))
                tunit=' s  '
                if est > 120:
                    est/=60
                    tunit=' min'
                sys.stdout.write('\r'+'     └─'+'{:6.2f}'.format(prog)+'%, estimation: '+'{:6.2f}'.format(est)+tunit)
                sys.stdout.flush()
                #stime=time.time()


            # concatenate all trajectories of the bin
            xl=lon[traid[idc],trastart[idc]:trastop[idc]].tolist()
            trajx+=xl
            trajx.append(0.)
            
            yl=lat[traid[idc],trastart[idc]:trastop[idc]].tolist()
            trajy+=yl
            trajy.append(0.)
            
            if zrel[traid[idc],trastop[idc]] == 0. :
                depox.append(lon[traid[idc],trastop[idc]])
                depoy.append(lat[traid[idc],trastop[idc]])
            
                
            zl=zrel[traid[idc],trastart[idc]:trastop[idc]].tolist()
            trajz+=zl
            trajz.append(None)

            tl=dates[trastart[idc]:trastop[idc]]
            ztime+=tl
            ztime.append(None)

            nt+=1


        # project to map
        xr,yr=m[0](trajx,trajy)
        x0,y0=m[0](0.,0.)
        
        # create 'holes' between the trajs
        xr=np.where(np.asarray(xr)==x0, None, xr )
        yr=np.where(np.asarray(yr)==y0, None, yr )
        
        xd,yd=m[1](depox,depoy)
        

        # trajectories on horizontal map
        m[0].plot(xr,yr,color=color,linewidth=.25)
        
        # scatter polt for the deposition points
        m[1].scatter(xd,yd,color=color,marker='.',s=5)#linewidth=.001)
        
        # plot_date for the time series
        ax2.plot_date(ztime,trajz,color=color,fmt='-',linewidth=.25)
            

# axis label and format
ax2.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
ax2.tick_params(labelsize='x-large')
ax2.set_ylabel('Altitude [m]',size='x-large')

# time range for the time series
ax2.set_ylim(0,10)
ax2.set_xlim(dates[102],dates[252])

# ax2.set_ylim(0)
# ax2.set_xlim(dates[82],dates[1442])

# colorbar
norm = mpl.colors.BoundaryNorm(bintop, cmap.N)
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,norm=norm,orientation='horizontal')
cb1.ax.tick_params(labelsize='x-large')
cb1.set_label('Particle size [$\mu m$]',size='x-large')



# Output
fname=datetime.datetime.now().strftime("%Y%m%d%H%M")+'.png'
plt.tight_layout()

mpl.rcParams['agg.path.chunksize'] = 10000
plt.savefig(fname,dpi=500)

print()
print('│  │ └─done plotting '+str(nt)+' trajs '+str(int(time.time()-ptime))+'s                             ')
print('│  │ ')
print('│  done '+str(int(time.time()-sbtime))+'s')
print('│')


# ## Mean trajectory
# In the next few cells I define a mean trajectory.
# I need this for [another project](https://github.com/mttfst/trajectory-cross-section).
# 
# 
# 
# I only want the trajectories that travel all the way to the east. 

# In[9]:


# define the trajs that travel far to the east -> max(rlon) >= 5
trajlens=[]
long_trajs=[]
for idc in traid:
    if rlon[traid[idc-1],trastop[idc-1]] >= 5:
        trajlens.append(trastop[idc-1]-trastart[idc-1]) #traid[idc],trastart[idc]:trastop[idc]
        long_trajs.append(idc-1) #traid[idc],trastart[idc]:trastop[idc]

print(np.percentile(trajlens,99))
maxlen=int(np.percentile(trajlens,99))

rlat_long = np.zeros((len(long_trajs),maxlen))
rlon_long = np.zeros((len(long_trajs),maxlen))
z_long    = np.zeros((len(long_trajs),maxlen))
zrel_long = np.zeros((len(long_trajs),maxlen))

i=0
starttimes=[]
print(rlat.shape)
for idc in long_trajs:
    rlat_long[i,:]=rlat[traid[idc],trastart[idc]:trastart[idc]+maxlen]
    rlon_long[i,:]=rlon[traid[idc],trastart[idc]:trastart[idc]+maxlen]
    z_long[i,:]=z[traid[idc],trastart[idc]:trastart[idc]+maxlen]
    zrel_long[i,:]=zrel[traid[idc],trastart[idc]:trastart[idc]+maxlen]
    starttimes.append(trastart[idc])
    
    i+=1


mean_start_time=int(np.median(starttimes))


# Then a calculate the mean and the standard deviation of my selection.

# In[10]:


mean=np.zeros((maxlen,4))
stda=np.zeros((maxlen,4))


for t in range(maxlen):
    mx =rlon_long[:,t]
    my =rlat_long[:,t]
    mz =z_long[:,t]
    mzr=zrel_long[:,t]

    mean[t,0]=np.mean(mx)
    mean[t,1]=np.mean(my)
    mean[t,2]=np.mean(mz)
    mean[t,3]=np.mean(mzr)

    # mean[t,0]=np.median(mx)
    # mean[t,1]=np.median(my)

    stda[t,0]=np.std(mx)
    stda[t,1]=np.std(my)
    stda[t,2]=np.std(mz)
    stda[t,3]=np.std(mzr)
    


# In the end, I write everything in an ASCII file.

# In[11]:



f = open('mean_traj.txt','w')
i=0
for l in mean:
    f.write(timestemps[mean_start_time+i]+'  '+str(l[0])+'  '+str(l[1])+'   '+str(l[2])+'   '+str(l[3])+'\n')
    i+=1
f.close()
i=0
f = open('stda_traj.txt','w')
for l in stda:
    f.write(timestemps[mean_start_time+i]+'  '+str(l[0])+'  '+str(l[1])+'   '+str(l[2])+'   '+str(l[3])+'\n')
    i+=1
f.close()

