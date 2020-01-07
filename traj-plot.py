#!/usr/bin/env python
# coding: utf-8

# In[2]:


from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap#, cm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
import os
import sys
import time
from matplotlib.ticker import PercentFormatter
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
from matplotlib.patches import Polygon
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
#import warnings
#warnings.filterwarnings("ignore")

def cm2inch(value):
    return value/2.54

def rot2geo(xin,yin):


    pollon=-170.
    pollat=40.
    xgeo=[i for i in xin]
    ygeo=[i for i in yin]

    for i in range(len(xin)):
        ygeo[i]=180./math.pi * math.asin(math.sin(math.pi/180.*yin[i])*math.sin(math.pi/180.*pollat) + math.cos(math.pi/180.*yin[i])*math.cos(math.pi/180.*xin[i])*math.cos(math.pi/180.*pollat))

        xgeo[i]=180./math.pi * math.atan((math.cos(math.pi/180.*yin[i])*math.sin(math.pi/180.*xin[i]))/(math.sin(math.pi/180.*pollat)*math.cos(math.pi/180.*yin[i])*math.cos(math.pi/180.*xin[i])-math.sin(math.pi/180.*yin[i])*math.cos(math.pi/180.*pollat))) + pollon + 180

    return(xgeo,ygeo)

def meter2grad(meter):
    grad=np.zeros(len(meter))
    for i in range(len(meter)):
        grad[i] = meter[i]/6.3710088e6 * 180./np.pi
    return grad

def grad2meter(grad):
    meter=np.zeros(len(grad))
    for i in range(len(grad)):
        meter[i] = grad[i]*6.3710088e6 * np.pi/180.
    return meter


# In[187]:


totaltime=time.time()
print()
print('Start')
print('│')
print('├──read data')
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

nfile=0
timesave=0
timediff=0
for fn in namelist[:]:
    nfile+=1
    sys.stdout.write('\r'+'│  │ '+fn)
    sys.stdout.flush()
    datestring=fn[5:-3]


    # data input
    f=Dataset(fn,'r')
    intime = f.variables['time'][:]
    inz    = f.variables['z'][:]          # [m]
    inzrel = f.variables['zrel'][:]          # [m]
    inlat  = f.variables['latitude'][:]   # [rotlat]
    inlon  = f.variables['longitude'][:]  # [rotlon]
    insize = f.variables['size'][:]  # [rotlon]
    indt   = f.getncattr('output_time_step_in_sec')
    f.close()


    if nfile == 1:
        #timesave=len(intime)
        timesave=intime[2]
        z=np.asarray(inz)
        zrel=np.asarray(inzrel)
        rlat=np.asarray(inlat)
        rlon=np.asarray(inlon)
        size=np.asarray(insize)
    else:
        #timediff=timesave-len(intime)
        timediff=int((intime[2]-timesave)/indt)
        if timediff != 0:
            dim=inz.shape

            oz    = np.zeros((timediff,dim[1]))
            ozrel = np.zeros((timediff,dim[1]))
            olat  = np.zeros((timediff,dim[1]))
            olon  = np.zeros((timediff,dim[1]))
            
            oz[:,:]   = inz[1][:]
            ozrel[:,:]= inzrel[1][:]
            olat[:,:] = inlat[0][:]
            olon[:,:] = inlon[0][:]
            

            inz[0][:]=inz[1][:]
            inzrel[0][:]=inzrel[1][:]

            inz   = np.concatenate((oz  ,np.asarray(inz)),   axis=0)
            inzrel= np.concatenate((ozrel  ,np.asarray(inzrel)),   axis=0)
            inlat = np.concatenate((olat,np.asarray(inlat)), axis=0)
            inlon = np.concatenate((olon,np.asarray(inlon)), axis=0)

        
        if np.asarray(inz).shape != z.shape:
            if np.asarray(inz).shape[0] > z.shape[0]:
                znew=np.zeros((np.asarray(inz).shape[0],z.shape[1]))
                znew[:z.shape[0],:]=z
                znew[z.shape[0]:,:]=z[-1,:]

                zrelnew=np.zeros(np.asarray(inzrel).shape)
                zrelnew[:zrel.shape[0],:]=zrel
                zrelnew[zrel.shape[0]:,:]=zrel[-1,:]

                rlatnew=np.zeros(np.asarray(inlat).shape)
                rlatnew[:z.shape[0],:]=rlat
                rlatnew[z.shape[0]:,:]=rlat[-1,:]

                rlonnew=np.zeros(np.asarray(inlon).shape)
                rlonnew[:rlon.shape[0],:]=rlon
                rlonnew[rlon.shape[0]:,:]=rlon[-1,:]

                z=np.concatenate((znew,np.asarray(inz)), axis=1)
                zrel=np.concatenate((zrelnew,np.asarray(inzrel)), axis=1)
                rlat=np.concatenate((rlatnew,np.asarray(inlat)), axis=1)
                rlon=np.concatenate((rlonnew,np.asarray(inlon)), axis=1)
                
            else:
                inznew=np.zeros((z.shape[0],inz.shape[1]))
                inznew[:inz.shape[0],:]=np.asarray(inz)
                inznew[inz.shape[0]:,:]=np.asarray(inz)[-1,:]
                
                inzrelnew=np.zeros((zrel.shape[0],inzrel.shape[1]))
                inzrelnew[:inzrel.shape[0],:]=np.asarray(inzrel)
                inzrelnew[inzrel.shape[0]:,:]=np.asarray(inzrel)[-1,:]
                
                inlatnew=np.zeros((rlat.shape[0],inlat.shape[1]))
                inlatnew[:inlat.shape[0],:]=np.asarray(inlat)
                inlatnew[inlat.shape[0]:,:]=np.asarray(inlat)[-1,:]
                
                inlonnew=np.zeros((rlon.shape[0],inlon.shape[1]))
                inlonnew[:inlon.shape[0],:]=np.asarray(inlon)
                inlonnew[inlon.shape[0]:,:]=np.asarray(inlon)[-1,:]
 
                z=np.concatenate((z,inznew), axis=1)
                zrel=np.concatenate((zrel,inzrelnew), axis=1)
                rlat=np.concatenate((rlat,inlatnew), axis=1)
                rlon=np.concatenate((rlon,inlonnew), axis=1)

        else:  
            z=np.concatenate((z,np.asarray(inz)), axis=1)
            zrel=np.concatenate((zrel,np.asarray(inzrel)), axis=1)
            rlat=np.concatenate((rlat,np.asarray(inlat)), axis=1)
            rlon=np.concatenate((rlon,np.asarray(inlon)), axis=1)
            
            
        size=np.concatenate((size,np.asarray(insize)), axis=0)
    
    
    print(' -> '+str(int(time.time()-stime))+'s')


print('│  done '+str(int(time.time()-stime))+'s')
print('│')


# In[229]:


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
for id in range(1,dimz[1],10):
    idc+=1
    traid.append(id)
 
    # identify startpoint of traj
    if not first:
        pvt=trastart[idc-1]  # previously start time
        
    if not first and z[pvt][id] != z[0][id] and z[pvt][id] != z[1][id] and z[pvt-1][id] == z[1][id]:
        trastart.append(pvt)
    else: 
        for t in range(1,dimz[0]):
            if z[t][id] != z[0][id] and z[t][id] != z[1][id]:
                first=False
                trastart.append(t)
                break
            
    # identify position when the particle is dead
    tstop=0
    btime=trastart[idc] # bottom time
    ttime=dimz[0]-1    # top time
    nloop=0
    while ttime-btime > 10: 
        nloop+=1
        jtime=int((btime+ttime)/2)
        if z[jtime,id] != z[jtime-1,id]: # -> alive
            if z[jtime,id] == z[jtime+1,id]: # -> dead
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
            if z[t,id] == z[t+1,id]:
                tstop=t
                break
                
    if tstop==0 and ttime == dimz[0]-1:
        tstop=ttime
        
    if trastart[idc]==tstop:
        trastart[idc]-=1
        
    if tstop > 0:            
        trastop.append(tstop)
        
    
    
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

#print(bins)
#print(traid)


# In[230]:


zrelmax=2500#np.max(zrel)

lonmin=12.
lonmax=40.
latmin=40.
latmax=54.

dfigx=lonmax - lonmin
if dfigx < 0.5:
    spcnglon=0.02
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


# In[234]:


print('├──create figure')
sbtime=time.time()

# create figure
fig = plt.figure(figsize=(9,12))
#gs = gridspec.GridSpec(nrows=4, ncols=1, height_ratios=[1,1,1,0.1])
gs = gridspec.GridSpec(nrows=3, ncols=1, height_ratios=[1,1,0.1])

# create subplots
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[1,0])
#ax2 = fig.add_subplot(gs[2,0])
ax3 = fig.add_subplot(gs[2,0])

# set titles
ax0.set_title('Horizontal Dispersion',size='x-large')
ax1.set_title('Deposition',size='x-large')
#ax2.set_title('Vertical Dispersion',size='x-large')

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
    map.drawcoastlines()
    map.drawcountries()
#ax2.grid()
# ax2.set_xlim(lonmin,lonmax)
#ax2.set_ylim(-0.1,zrelmax)
# ax2.xaxis.set_major_formatter(FormatStrFormatter('%5.2f°E'))
#ax2.set_ylabel('height [m]',size='x-large')

de='DEU_adm_shp/'
cz='CZE_adm_shp/'
pl='POL_adm_shp/'
# add shapefiles
shapes = [de+'DEU_adm4',de+'DEU_adm3',de+'DEU_adm2',de+'DEU_adm1',de+'DEU_adm0',
          cz+'CZE_adm1',cz+'CZE_adm0',
          pl+'POL_adm1',pl+'POL_adm0']

linewidths = [.1,.2,.3,.5,1.,.5,1.,.5,1.]
colors = ['gray','gray','gray','black','red','black','red','black','red']

#for shape,line,color in zip(shapes[adlvl:],linewidths[adlvl:],colors[adlvl:]):
#    for map in m:
#        shp = map.readshapefile(shape_path+shape, 'states', drawbounds=True,linewidth=line,color=color)

print('│  │ └─done '+str(int(time.time()-mtime))+'s')
print('│  │ ')

print('│  ├─plot data')
ptime=time.time()


ttot=0.
ncolor=0
nt=0

from matplotlib.colors import LinearSegmentedColormap
# cmap = LinearSegmentedColormap.from_list(name='mycmap',colors =['cyan', 'mediumorchid'],N=31)
# cmap = LinearSegmentedColormap.from_list(name='mycmap',colors =['cyan', 'darkorange'],N=31)
# cmap = LinearSegmentedColormap.from_list(name='mycmap',colors =['#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c', '#4d004b'],N=31)
# cmap = LinearSegmentedColormap.from_list(name='mycmap',colors =['#4292c6','#6baed6','#9ecae1','#c6dbef','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c', '#f16913'],N=31)
# cmap = LinearSegmentedColormap.from_list(name='mycmap',colors =['cyan','#6baed6','#9ecae1','#c6dbef','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c', '#f16913'],N=31)
cmap = LinearSegmentedColormap.from_list(name='mycmap',colors =['cyan','lightskyblue','#9ecae1','#c6dbef','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c', '#f16913','#d94801'],N=31)
# cmap = LinearSegmentedColormap.from_list(name='mycmap',colors =['cyan','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c', '#f16913','#d94801'],N=31)

stime=time.time()
#for bx,by,bz,bzr in bins: #range(nbin):
for bin in bins:
    color = cm.cool(ncolor)
    ncolor+=8
    #color = cmap(ncolor)
    #ncolor+=1
    
    if len(bin)>0:
        for idc in bin:
            # progres bar
            if nt > 0 and np.mod(nt,100):
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


            # plot
            trajx=rlon[trastart[idc]:trastop[idc],traid[idc]]
            trajy=rlat[trastart[idc]:trastop[idc],traid[idc]]
            trajz=zrel[trastart[idc]:trastop[idc],traid[idc]]
            
            if len(trajz) == 0:
                print(bin,len(bin))
                print(idc,traid[idc],trastart[idc],trastop[idc])
                print(zrel[trastart[idc]:trastop[idc],traid[idc]])
                print(zrel[:,traid[idc]])

            lon_r, lat_r = rot2geo (trajx,trajy)
            xr,yr=m[0](lon_r,lat_r)
            m[0].plot(xr,yr,color=color,linewidth=.25)
            if trajz[-1] < 50.:
                m[1].scatter(xr[-1],yr[-1],color=color,marker='.',s=5)#linewidth=.001)
            # if len(trajx) == len(trajz):
            # ax2.plot(lon_r,trajz,color=color,linewidth=.25)
#            ax2.plot(np.arange(len(trajz)),trajz,color=color,linewidth=.25)
            nt+=1
            #stime=time.time()

#ax2.set_xticklabels(['0','09:40','10:00','10:20','10:40','11:00'],size='x-large')
#ax2.set_yticklabels([0,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0],size='x-large')
ax0.set_xticklabels([14.18,14.2,14.22,14.24,14.26,14.28,14.3],size='x-large')
ax0.set_yticklabels(ax0.get_yticklabels(),size='x-large')
cmap = mpl.cm.cool
norm = mpl.colors.BoundaryNorm(bintop, cmap.N)
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,norm=norm,orientation='horizontal')
cb1.ax.tick_params(labelsize='x-large')
cb1.set_label('Particle size [$\mu m$]',size='x-large')

import datetime

fname=datetime.datetime.now().strftime("%Y%m%d%H%M")+'.png'
plt.tight_layout()
plt.savefig(fname,dpi=500)


print('│  │ └─done '+str(int(time.time()-ptime))+'s                             ')
print('│  │ ')
print('│  done '+str(int(time.time()-sbtime))+'s')
print('│')


# In[232]:



import datetime
print(datetime.datetime.now().strftime("%Y%m%d%H%M"))


# In[ ]:





# In[ ]:




