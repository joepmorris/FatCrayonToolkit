#!/usr/local/bin/python
from __future__ import print_function
import numpy as np
import copy as copy
import math
from copy import deepcopy
import re
import os, sys; sys.path.append(os.path.dirname(__file__)); sys.path.append("./")
from FCT import SimpleGeometry as sg
from FCT import DataReader as dr
from FCT import StressCalc as sc
from FCT.Units import *
import os
import pandas as pd

#===============================================================================
# This section contains the parameters that most users will wish to modify
#===============================================================================

# The layout file that we will write to (and then read back from)
fatCrayonBoreholes='./Option_clock.txt'; vtkDir='./'

# Length of injection/production wells
lInj=265.0-15.0;# lInj=0.75*265.0#feet
# Orientation and location of the injection well
inj={'L':lInj,'dec':15,'az':50.0,'x0':np.asarray([3962.17,   -3003,   0.0])}
# The depth along the injector we intend to stimulate
stimDepth = 195#feet

# Distance between toe of injector and the toes of the production wells
toeDistance=33.0#feet
# The user provides a list of production boreholes at orientations relative to the drift direction (defined as v12 later)
# clockAz of 0 will place a production well on the side of the injection well closest to the drift
# clockAz of 180 will place a production well on the side of the injection well in direction opposite from the drift
#clockAzs=[-100,20,140]
#clockAzs=np.asarray([-120,0,120])+180
#clockAzs=[0,90,180,270]
clockAzs=[0,90,270]

# Dimensions of the excavations:
excavations=[
    # Excavation by the injector
    {'H':7.0,'L':20.0,'W':4.0,'loc':[3962.17, -3000, 0.0]}
    ]

# We define some locations to set up monitoring and we define key target parameters for each
monStations=[]
if True:
    # Monitoring collars are in the new alcove, introduce boreholes that are angled about the injector
    # like hands on a clock, similar to the arrangement of the production boreholes
    monStations.append({
        # Stipulate the distance from the injector at the toe of the injector - you will want this to be larger than toeDistance define above for producers
        'toeDistance':toeDistance,
        # Collar location
        'x0':inj['x0'],
        # Length
        #'L':285.0,
        'L':lInj+15.0,
        # Monitoring well azimuths using the same orientations as defined above for production boreholes
        'clockAzs':[180],
        #'clockAzs':[-100,20,140],
        #'clockAzs':[0,90,180,270],
        # Name for this station
        'name':'Inj'
        })
if True:
    # The south station
    # This has traditionally been used with the north station, but could also be used with the alcove station
    monStations.append({
        # Stipulate the distance of closest approach to injector by the monitoring wells
        'monDistance':35.0,#feet
        # Stipulate the distance along the injector where the monitoring well makes closest approach
        'monInjDepth':130.0+5.0,#feet
        # Collar location
        'x0':np.asarray([3965.0,   -2864.5,   0.0]),
        # Length
        #'L':185.0,
        'L':180.0,
        # Name for this station
        'name':'S'
        })
    # Excavation by the injector
    excavations.append( {'H':7.0,'L':3.0,'W':3.0,'loc':[3963.44,   -2862.5, 0.0]} )
if True:
    monStations.append({
        # Stipulate the distance of closest approach to injector by the monitoring wells
        'monDistance':40.0,#feet
        # Stipulate the distance along the injector where the monitoring well makes closest approach
        #'monInjDepth':230.0,#feet
        'monInjDepth':210.0+20.0,#feet
        # Collar location
        'x0':np.asarray([3994.92, -2834.25, 0.0]),
        # Length
        'L':195.0,
        # Name for this station
        'name':'A'
        })

# This is where we define parameters related to the excavation dimensions and the workspace requirements:
# We need 12' of workspace
# for monitoring wells:
wspc_mon={'L':12.0,'W':5.0,'H':7.0}
# AND injection and production wells:
wspc_inj_prd={'L':12.0,'W':5.0,'H':7.0}


J={
# Joint orientations
    "1":{"strike": 155.0, "dip":55.0, "r":50.0,"h":0.2}, # Tom's JS1 - mu=0.44 most prone to slip with 2020-06-20 interp of stress
    "2":{"strike": 20.0,  "dip":27.0, "r":50.0,"h":0.2}, # Tom's JS2 - mu=0.22
    "3":{"strike": 280.0, "dip":50.0, "r":50.0,"h":0.2}, # Tom's JS3 - mu=0.20
    "4":{"strike": 265.0, "dip":70.0, "r":50.0,"h":0.2}, # Tom's JS4 - mu=0.32
    "5":{"strike": 50.0,  "dip":40.0, "r":50.0,"h":0.2}, # Bill's sulfide - mu=0.35
}

# The in situ stress:
Sh=1.83e+07; SH=3.73e+07; SV=3.6e+07; ShAzimuthDeg=24.0; ShDipDeg=28.0; Pp=4.23e+06
# Hydraulic fracture orientation:
J["0"]={"strike": 114.0 , "dip":62.0, "r":50.0,"h":0.2}

# This only apply if you have the relevant point clouds available to read in:
include_scan = not True # Do we include the (time consuming) point cloud scan of the drift?
include_rhyolite = not True

#===================================================
# Everything below here is largely automatic
#===================================================
# Origin for our model:
o=np.asarray([3957,-2829,0.0])
# The vector aligned with the injection well
vInj=sg.vectFromAzDecDeg(azDeg=inj['az'], decDeg=inj['dec'])

if True:
    # Stochastic fractures:
    treatStochasticAsDeterministic=False
    # This is where we control the stochastic natural fracture statistics
    # Add in variability only here (the joint sets were already included above)
    # Hydraulic fracture orientation based on TV4100 test results, four hydraulic fractures in Amphibolite, shallow 
    J["0"]["dstrike"]=0.0; J["0"]["ddip"]=0.0
    # Tom Doe orientations from Jeff
    J["1"]["dstrike"]=20.0e0;J["1"]["ddip"]=20.0e0 # Tom's JS1
    J["2"]["dstrike"]=20.0e0;J["2"]["ddip"]=50.0e0 # Tom's JS2
    J["3"]["dstrike"]=20.0e0;J["3"]["ddip"]=50.0e0 # Tom's JS3
    J["4"]["dstrike"]=20.0e0;J["4"]["ddip"]=50.0e0 # Tom's JS4
    J["5"]["dstrike"]=10.0e0;J["5"]["ddip"]=100.0e0 # Bill's sulfide
    for key in J.keys():
        JS=J[key]
        # These will be handled by the global stochastic properties
        JS.pop("r",None)
        JS.pop("h",None)
    # Generate many fractures within a volume
    nFrac=100
    #nFrac=400
    #nFrac=800
    # Corners of box that fractures occupy - let's keep these away from the alcove where we have deterministic fractures
    # Designs now penetrate further to the east
    xbox0=o+[ 0.0,  -170.0, -150.0]
    xbox1=o+[ 250.0, 150.0,    0.0]
    #
    #xbox0=o+[-150.0,-200.0,-200.0]
    #xbox1=o+[ 150.0, -50.0,   0.0]
    ## For design A
    #xbox0=o+[ 0.0,  -150.0, -300.0]
    #xbox1=o+[ 150.0, 150.0,    0.0]
    # For design A2
    #xbox0=o+[ 0.0,   050.0, -200.0]
    #xbox1=o+[ 200.0,-250.0,    0.0]
    # The radius of the stochastic fractures
    r0_stochastic=30.0
    h0_stochastic=0.2
else:
    treatStochasticAsDeterministic=True
    # Number of fractures to include per set
    # This is normally for generating many fractures within a volume
    # In this instance we are just going to visualize one of each, colocated at the target stimulation depth
    nFrac=1
    xbox0=inj['x0']+stimDepth*vInj
    xbox1=xbox0+[1e-3, 1e-3, 1e-3]
    # We need to add in the variability that would normally be included in the J dictionary,
    # but isn't meaningful because we are only plotting one orientation
    for key in J.keys():
        JS=J[key]
        JS["dstrike"]=1e-3
        JS["ddip"]=1e-3
    r0_stochastic=30.0
    h0_stochastic=0.2

#Place caging production wells, perfect phasing, favoring two on drift side. Perhaps not perfect phasing, but rather favor protecting the drift.
#Could parameterize the caging boreholes in terms of azimuth relative to the orientation that hits the drift
#2' spacing between wells is probably sufficient

# Orientation of the drift in this area
driftOrientDeg=0.0
print('Assuming drift is oriented at azimuth of',driftOrientDeg,'degrees')
# Vector aligned with the drift:
vDrift=sg.vectFromAzDecDeg(driftOrientDeg,0.0)

# We then get an angle from the injection well
prodAngRad=2.0*math.asin(0.5*toeDistance/inj['L'])
prodAngDeg=sg.radToDeg(prodAngRad)
print('angle from injector to producers is',prodAngDeg,'degrees to get proximity of',toeDistance,'feet from injector')

# We then get the plane of the injection well and the drift - we use this as a reference orientation
# I want a vector lying in the vInj-vDrift plane, perpendicular to vInj, pointing towards drift
# Let's call that vector 12 for 12 o'clock
v12=sg.normalize(vDrift-vInj)

# List of production wells
prods=[]
nProd=0
for clockAz in clockAzs:
    prod={'L':lInj,'x0':inj['x0']}
    # We need to set the 'dec':??,'az':?? for this production well
    # We start with a production well parallel to the injection well
    # We need to rotate vProd to be prodAngDeg away from vInj, but we need to figure axis of rotation
    # The axis of rotation is the cross product of the injection borehole and clock hand
    # We need to rotate the clock hand about vInj
    vHand = sg.rotatePoints( points=v12, axis=vInj, theta=sg.degToRad(-clockAz) )
    # We now rotate the injection well away from the production well
    # The axis of rotation is the cross product of the injection borehole and clock hand
    axis=np.cross(vInj,vHand)
    vProd = sg.rotatePoints( points=vInj, axis=axis, theta=prodAngRad )
    # Now convert this direction vector into azimuth and declination
    #print vProd
    prod['az'], prod['dec'] = sg.azDecDegFromVect( vProd )
    nProd+=1
    prod['name']='Prod %d'%(nProd)
    # And append this production well to the list
    prods.append(prod)

# Generate the monitoring wells
# For each station we generate an upper and lower monitoring well
# We know the collar location and the depth of closest approach to the injector
# as well as the distance above and below we wish to pass by
mons=[]
nMon=0
for stn in monStations:
    if 'clockAzs' in stn:
        # This monitoring station is described in terms of azimuths about the injector
        # We then get an angle from the injection well - perhaps based upon distance between toes?
        monAngRad=2.0*math.asin(0.5*stn['toeDistance']/inj['L'])
        monAngDeg=sg.radToDeg(monAngRad)
        print('angle of monitoring wells at station',stn['name'],'to injection well is',monAngDeg,'degrees to get proximity of',stn['toeDistance'],'feet from injector toe')
        nMon=0
        for clockAz in stn['clockAzs']:
            mon={'L':stn['L'],'x0':stn['x0']}
            # We need to set the 'dec':??,'az':?? for this monitoring well
            # We start with a monitoring well parallel to the injection well

            # We need to rotate vMon to be monAngDeg away from vInj, but we need to figure axis of rotation
            # The axis of rotation is the cross product of the injection borehole and clock hand
            # We need to rotate the clock hand about vInj
            vHand = sg.rotatePoints( points=v12, axis=vInj, theta=sg.degToRad(-clockAz) )
            # We now rotate the monitoring well away from the injection well
            # The axis of rotation is the cross product of the injection borehole and clock hand
            #print vInj, vHand
            axis=np.cross(vInj,vHand)
            vMon = sg.rotatePoints( points=vInj, axis=axis, theta=monAngRad )
            # Now convert this direction vector into azimuth and declination
            mon['az'], mon['dec'] = sg.azDecDegFromVect( vMon )
            nMon+=1
            mon['name']='Mon %s %d'%(stn['name'],nMon)
            # And append this monitoring well to the list
            mons.append(mon)
    else:
        monUp={'L':stn['L'],'x0':stn['x0']}
        monDown={'L':stn['L'],'x0':stn['x0']}
        # Target location along the injector
        trg=inj['x0']+stn['monInjDepth']*vInj
        # Vector from collar location to the target on the injector
        v=trg-stn['x0']
        dist=math.sqrt(np.dot(v,v))
        v=v/dist
        # We need to get the angle to rotate to get the required distance above and below at the target point
        monAngRad=2.0*math.asin(0.5*stn['monDistance']/dist)
        monAngDeg=sg.radToDeg(monAngRad)
        print('angle up/down for',stn['name'],'station is',monAngDeg,'degrees to get proximity of',stn['monDistance'],'feet from injector')
        # We now rotate this up and down to get the monitoring pair
        # We need to carefully choose what axis to rotate around
        # The axis should not be the injector unless the injector is perpendicular to the monitoring well
        # The axis is perpendicular to the monitoring well and lies in the plane of the monitoring well and the injector
        # I can remove the component of the injector that is parallel to monitoring and then renormalize
        # This is a linear combination of vIng and v, so it lies in their plane:
        axis = vInj - np.dot(vInj,v)*v
        axis = sg.normalize(axis)
        # Now rotate v down and up
        vDown=sg.rotatePoints( points=v, axis=axis, theta=-monAngRad )
        monUp['az'], monUp['dec'] = sg.azDecDegFromVect( vDown );
        nMon+=1; monDown['name']='Mon %sD'%(stn['name'])
        if nMon==1:
            monDown['name']+=' (1st)'
        vUp=sg.rotatePoints( points=v, axis=axis, theta=monAngRad )
        nMon+=1; monUp['name']='Mon %sU'%(stn['name'])
        monDown['az'], monDown['dec'] = sg.azDecDegFromVect( vUp )
        mons.append(monUp)
        mons.append(monDown)
    
# We will write the layout to the layout file first, then read it back!
fd=open(fatCrayonBoreholes,'w')
fd.write(
"""BoreType,Length_ft,Declination_deg,Azimuth_deg,Easting_ft,Northing_ft,TVD_ft,comment
""")
fd.write("I,       %.1f,    %.1f,           %.1f,       %.1f,   %.1f,   %.1f, 'Inj (2nd)'\n"%
         ( inj['L'], inj['dec'], inj['az'], inj['x0'][0],inj['x0'][1],inj['x0'][2] )
        )
nMon=0
for mon in mons:
    fd.write("M,       %.1f,    %.1f,           %.1f,       %.1f,   %.1f,   %.1f, '%s'\n"%
             ( mon['L'], mon['dec'], mon['az'], mon['x0'][0],mon['x0'][1],mon['x0'][2],mon['name'] )
            )
nProd=0
for prod in prods:
    fd.write("P,       %.1f,    %.1f,           %.1f,       %.1f,   %.1f,   %.1f, '%s'\n"%
             ( prod['L'], prod['dec'], prod['az'], prod['x0'][0],prod['x0'][1],prod['x0'][2],prod['name'] )
            )
fd.close()


print('\n\nWe are using:\nSh=%.3g; SH=%.3g; SV=%.3g; ShAzimuthDeg=%.3g; ShDipDeg=%.3g; Pp=%.3g'%(Sh, SH, SV, ShAzimuthDeg, ShDipDeg, Pp));



def getFrac(x0,e0,e1):
    v0=np.asarray(e0)-np.asarray(x0)
    v1=np.asarray(e1)-np.asarray(x0)
    n=np.cross(v0,v1)
    n=n/np.linalg.norm(n)
    print('n',n)
    dipDirection,dipAngle=sg.dipDirectionAndDipAng(n)
    return x0, (dipDirection,dipAngle+90.0)

# Deterministic fractures if you have them:
deterministicFracs={
}
# The radius of the deterministic fractures
r0_deterministic=20.0*0.025
h0_deterministic=0.2

#rBore=0.5 # Easier to see
rBore=1.0 # Easier to see
#rBore=0.2 # Need to be a bit more precise

#===============================================================================
# Most users will probably not want to modify below this point
#===============================================================================

# The maximum mu observed in nature is alleged to be 0.6 (Byerlee et al, 1978)
muSlip=0.6

nf_dip_dir_deg=[]; nf_dip_deg=[]
if not True:
    # Plot more traditional stress diagrams for stableMu/slipTendancy and criticalDeltaPp
    sc.plotCritPpStableMu(Sh,SH,SV,ShAzimuthDeg,ShDipDeg,Pp,
                          muSlip,
                          nf_dip_dir_deg, nf_dip_deg,
                          critDelPpFile="Critical_Del_Pp_noFracs.png", stableMuFile='Stable_Mu_Slip_Tendancy_noFracs.png', mohrsFile="Mohrs_Circle_noFracs.png",)

# To get reproducible, pseudorandom results, we set the seed:
np.random.seed(1)

# These are various collections of 3D objects for subsequent visualization
drift_obj=None # This would be where you import a 3D model of the drift itself
injection_obj=None;
production_obj=None;
monitoring_obj=None;
hfOpt_obj=None; # Our original design hydraulic fractures
hfLim_obj=None # Larger, hopefully largest plausible hydraulic fractures
hfNotch_obj=None # The actual notch locations
labels_obj=None
natFracs_obj=None


drift_scan=None
rhyolite=None
if include_scan:
    drift_scan=sg.pointCloudToObj(drift_scan,pd.read_csv("your_point_cloud.csv"),0.01,0.2)
if include_rhyolite:
    rhyolite=sg.pointCloudToObj(rhyolite,pd.read_csv("your_point_cloud.csv"),1.0,2.0)
    
#Import proposed geometry
cwd = os.getcwd()

# May want to correct the coordinate system being used:
correction=[0.0,0.0,0.0]

print('Table of boreholes:')
bores = dr.ProLayout(fatCrayonBoreholes,name='Layout',correction=correction,format='csv',skiprows=0)
print('')
print('Summary of boreholes:')
print('Length ranges from',np.min(bores.l),'to',np.max(bores.l),'feet')
print('Total borehole length',np.sum(bores.l),'feet')
print('Declination ranges from',np.min(bores.d),'to',np.max(bores.d),'degrees')
print('')

# Orientation of the drift
driftAzDeg=0.0
alcove_e=3961.68; alcove_n=-2833.7
print('N-S Drift azimuth (deg)',driftAzDeg)
print('south-west corner of alcove   easting',alcove_e,'ft','   northing',alcove_n,'ft')

print('%15s%15s%25s%25s'%('name/comment','azimuth(deg)','angle(deg) to drift','dist frm alcove(ft)'))
for i in range(len(bores.e)):
    angToDrift=(bores.a[i]-driftAzDeg)
    if angToDrift>90.0:
        angToDrift=180.0-angToDrift
    print('%15s%12.1f   %20.1f     %20.1f'%(bores.cmt[i],bores.a[i],angToDrift,(np.sqrt( (bores.e[i]-alcove_e)**2 + (bores.n[i]-alcove_n)**2 ))))
#quit()


def genWorkspace(xo, xf, L,W,H):
    # Generate a workspace associated with the given borehole
    if True:
        # We return a workspace described with a cylinder
        #Note syntax: sg.cylObj(x0,x1,r)
        # We want a normal in direction xo->xf
        # We add a slight pertubation to avoid singularity
        ds=xf-xo+np.asarray([1e-3,0.0,0.0])
        ds/=np.linalg.norm(ds)
        #return sg.cylObj(xo-L*ds,xo,0.5*W)
        return sg.cylObj(xo-L*ds,xo,0.5)
    else:
        # We return a workspace described with a box
        # It will have
        # * length L in the lateral direction of xo->xf
        # * width W in the lateral direction perpendiculat to xo->xf
        # * height H
        # We want a normal in direction xo->xf
        # We add a slight pertubation to avoid singularity
        ds=xf-xo+np.asarray([1e-3,0.0,0.0])
        ds[2]=0.0
        sL=(ds)/np.linalg.norm(ds)
        # But lying in the plane
        #print sL
        sL[2]=0.0
        sH=np.asarray([0.0,0.0,1.0])
        sW=np.cross(sL,sH)
        sL*=L; sH*=H; sW*=W
        wp=np.asarray([xo[0],xo[1],site.zInvert-0.01]) # Place slightly through floor for visualization
        # We will have the workspace end at the collar
        #print 'sL',sL; print 'sW',sW; print 'sH',sH; quit()
        return sg.convexFromPoints(
            np.asarray([
                wp-0.5*sW, wp+0.5*sW, wp-0.5*sW-sL, wp+0.5*sW-sL,
                wp-0.5*sW+sH, wp+0.5*sW+sH, wp-0.5*sW-sL+sH, wp+0.5*sW-sL+sH,
                ])
            ) # This turns a list of points into a polyhedron
    
#Include borehole geometries
injection_obj, production_obj, monitoring_obj = bores.appendVis(injection_obj, production_obj, monitoring_obj, rBore)

# Convert from principal to global
sigG=sc.sigGfromPrincipal(Sh,SH,SV,ShAzimuthDeg,ShDipDeg)

print("Looping over fracture sets - looking for intersections with boreholes")
for key in sorted(J.keys()): # This will loop through all of the joint sets
    # For each joint set we will create a separate dictionary that we will use to store all natural fracture information
    # In retrospect, this is probably overkill! We can write the VTK files as we go, rather than just save everything
    JS=J[key]
    print("JS",key)
    JS["nf"]={}
    NF=JS["nf"] # Handy link to the dictionary for the current natural fracture set
    NF["strike_deg"]   =[]
    NF["dip_dir_deg"]  =[]
    NF["dip_deg"]      =[]
    NF["enu"]          =[]
    NF["r"]            =[]
    NF["h"]            =[]
    NF["link"]         =[] # Store our calculation of whether this fracture links injectors and producers
    # We want the natural fractures to potentially all have different colors for subsequent VTK visualization
    NF["obj"]       =[]
    NF["stableMu"]     =[]
    NF["criticalDelPp"]=[]
    for iFrac in range(nFrac):
        strike_deg=np.random.normal(loc=JS["strike"],scale=JS["dstrike"])
        dip_deg=np.random.normal(loc=JS["dip"],scale=JS["ddip"])
        x0=np.asarray([ np.random.uniform(xbox0[0],xbox1[0]),np.random.uniform(xbox0[1],xbox1[1]),np.random.uniform(xbox0[2],xbox1[2]) ])
        dip_dir_deg=sg.strikeToDipDeg(strike_deg)
        #print JS["strike"],strike_deg,dip_dir_deg
        NF["strike_deg"].append(strike_deg)
        NF["dip_dir_deg"].append(dip_dir_deg)
        NF["dip_deg"].append(dip_deg)
        NF["enu"].append(x0)
        if "r" in JS:
            r=JS["r"]
        else:
            r=r0_stochastic
        if "h" in JS:
            h=JS["h"]
        else:
            h=h0_stochastic
        NF["r"].append(r)
        NF["h"].append(h)
        # Check if this fracture intersects both an injector and a producer
        link=0.0 # 0=no flow, 0.5=intersects one injector or producer, 1.0=intersects both
        for i in bores.iInjection:
            intersects,xInt = sg.diskIntersectSegmentEndPts(x0=bores.xo[i],x1=bores.xf[i], c0=x0,strikeDeg=strike_deg,dipDeg=dip_deg,r=r)
            #print intersects,s,aa,bb,xInt
            if intersects:
                link+=0.5
                break
        for i in bores.iProduction:
            intersects,xInt = sg.diskIntersectSegmentEndPts(x0=bores.xo[i],x1=bores.xf[i], c0=x0,strikeDeg=strike_deg,dipDeg=dip_deg,r=r)
            #print intersects,s,aa,bb,xInt
            if intersects:
                link+=0.5
                break
        NF["link"].append(link)
        #for dip_dir_deg,dip_deg,enu,r,h in zip(nf_dip_dir_deg,nf_dip_deg,nf_enu,nf_r,nf_h):
        # Get normal to the fracture in global coordinates
        nrmG=sc.nrmFromDipDirDipAngDeg(dip_direction=dip_dir_deg,dip_angle=dip_deg)
        # Now we can calculate the change in pore pressure for slip and the coefficient of friction necessary
        # prevent slipping. Conventional wisdom says we are in trouble if stableMu approaches muSlip=0.6
        criticalDelPp, stableMu = sc.criticalDelPp_stableMu(nrmG, sigG, Pp, muSlip)
        NF["obj"].append(sg.HF( r=r, x0=x0, strikeRad=sg.degToRad(sg.dipToStrikeDeg(dip_dir_deg)), dipRad=sg.degToRad(dip_deg), h=h) )
        NF["stableMu"].append(stableMu)
        NF["criticalDelPp"].append(criticalDelPp)
        # Let's remember these later for plotting
        nf_dip_dir_deg.append(dip_dir_deg); nf_dip_deg.append(dip_deg)


# Now we do the same for the deterministic fractures
for key in sorted(deterministicFracs.keys()): # This will loop through all of the joint sets
    # Get the deterministic fractures associated with this key
    frcs=deterministicFracs[key]
    # We will add these deterministic sets alongside the other sets in J
    if key in J.keys():
        print("Somehow we have an deterministic jointset with the same name as one of the stochastic ones")
        quit()
    J[key]={}
    JS=J[key]
    print("JS",key)
    JS["nf"]={}
    NF=JS["nf"] # Handy link to the dictionary for the current natural fracture set
    NF["dip_dir_deg"]  =[]
    NF["dip_deg"]      =[]
    NF["enu"]          =[]
    NF["r"]            =[]
    NF["h"]            =[]
    NF["link"]         =[] # Store our calculation of whether this fracture links injectors and producers
    # We want the natural fractures to potentially all have different colors for subsequent VTK visualization
    NF["obj"]       =[]
    NF["stableMu"]     =[]
    NF["criticalDelPp"]=[]
    for iFrac in range(len(frcs.e)):
        strike_deg=frcs.strike_deg[iFrac]
        dip_deg=frcs.dip_deg[iFrac]
        x0=np.asarray([frcs.e[iFrac],frcs.n[iFrac],frcs.u[iFrac]])
        dip_dir_deg=frcs.dip_dir_deg[iFrac]
        NF["dip_dir_deg"].append(dip_dir_deg)
        NF["dip_deg"].append(dip_deg)
        NF["enu"].append(x0)
        r=r0_deterministic
        try:
            r=frcs.r[iFrac]
        except:
            pass
        h=h0_deterministic
        NF["r"].append(r)
        NF["h"].append(h)
        # Check if this fracture intersects both an injector and a producer
        link=0.0 # 0=no flow, 0.5=intersects one injector or producer, 1.0=intersects both
        for i in bores.iInjection:
            intersects,xInt = sg.diskIntersectSegmentEndPts(x0=bores.xo[i],x1=bores.xf[i], c0=x0,strikeDeg=strike_deg,dipDeg=dip_deg,r=r)
            #print intersects,s,aa,bb,xInt
            if intersects:
                link+=0.5
                break
        for i in bores.iProduction:
            intersects,xInt = sg.diskIntersectSegmentEndPts(x0=bores.xo[i],x1=bores.xf[i], c0=x0,strikeDeg=strike_deg,dipDeg=dip_deg,r=r)
            #print intersects,s,aa,bb,xInt
            if intersects:
                link+=0.5
                break
        NF["link"].append(link)
        #for dip_dir_deg,dip_deg,enu,r,h in zip(nf_dip_dir_deg,nf_dip_deg,nf_enu,nf_r,nf_h):
        # Get normal to the fracture in global coordinates
        nrmG=sc.nrmFromDipDirDipAngDeg(dip_direction=dip_dir_deg,dip_angle=dip_deg)
        # Now we can calculate the change in pore pressure for slip and the coefficient of friction necessary
        # prevent slipping. Conventional wisdom says we are in trouble if stableMu approaches muSlip=0.6
        criticalDelPp, stableMu = sc.criticalDelPp_stableMu(nrmG, sigG, Pp, muSlip)
        NF["obj"].append(sg.HF( r=r, x0=x0, strikeRad=sg.degToRad(sg.dipToStrikeDeg(dip_dir_deg)), dipRad=sg.degToRad(dip_deg), h=h) )
        NF["stableMu"].append(stableMu)
        NF["criticalDelPp"].append(criticalDelPp)
        # Let's remember these later for plotting
        nf_dip_dir_deg.append(dip_dir_deg); nf_dip_deg.append(dip_deg)

        
# Give people a sense of orientation by adding a large, north-pointing arrow
#labels_obj=sg.mergeObj(labels_obj,
#                   sg.transObj( sg.scaleObj(sg.unitArrowY,[20.0,20.0,0.1]), o+[50.0,250.0,0.0] )
#                   )

# Visualize all the workspaces
if True:
    for i in range(len(bores.xo)):
        # We don't want to visualize workspace for TV
        print("%d '%s'"%(i,bores.cmt[i]))
        if 'Mon' in bores.cmt[i]:
            L=wspc_mon['L']; W=wspc_mon['W']; H=wspc_mon['H']
        else:
            L=wspc_inj_prd['L']; W=wspc_inj_prd['W']; H=wspc_inj_prd['H']
        labels_obj=sg.mergeObj(labels_obj, genWorkspace(bores.xo[i],bores.xf[i],L,W,H) )

if True:
    for exc in excavations:
        # Add an excavation
        H=exc['H']; L=exc['L']; W=exc['W']
        # This is the direction of L down the drift: Rotate from north to requested azimuth
        az_deg=180.0 # This should be updated to be parallel to drift
        dL = sg.rotatePoints( points=[np.asarray([0.0,1.0,0.0])], axis=np.asarray([0.0,0.0,1.0]), theta=-az_deg*np.pi/180.0 )[0]
        dW = sg.rotatePoints( points=[1.0*dL], axis=np.asarray([0.0,0.0,1.0]), theta=np.pi/2.0 )[0]
        dH = sg.rotatePoints( points=[dW], axis=np.asarray([0.0,1.0,0.0]), theta=-np.pi/2.0 )[0]
        print('L',dL); print('W',dW); print('H',dH)
        x0=np.asarray(exc['loc'])-1.0*dL-3.0*dW
        obj=None
        r=0.25
        obj = sg.mergeObj( obj, sg.cylObj(x0,x0+L*dL,r) )
        obj = sg.mergeObj( obj, sg.cylObj(x0+W*dW,x0+W*dW+L*dL,r) )
        obj = sg.mergeObj( obj, sg.cylObj(x0+H*dH,x0+H*dH+L*dL,r) )
        obj = sg.mergeObj( obj, sg.cylObj(x0+H*dH+W*dW,x0+H*dH+W*dW+L*dL,r) )

        obj = sg.mergeObj( obj, sg.cylObj(x0,x0+W*dW,r) )
        obj = sg.mergeObj( obj, sg.cylObj(x0+L*dL,x0+L*dL+W*dW,r) )
        obj = sg.mergeObj( obj, sg.cylObj(x0+H*dH,x0+H*dH+W*dW,r) )
        obj = sg.mergeObj( obj, sg.cylObj(x0+L*dL+H*dH,x0+L*dL+H*dH+W*dW,r) )

        obj = sg.mergeObj( obj, sg.cylObj(x0,x0+H*dH,r) )
        obj = sg.mergeObj( obj, sg.cylObj(x0+L*dL,x0+L*dL+H*dH,r) )
        obj = sg.mergeObj( obj, sg.cylObj(x0+W*dW,x0+W*dW+H*dH,r) )
        obj = sg.mergeObj( obj, sg.cylObj(x0+L*dL+W*dW,x0+L*dL+W*dW+H*dH,r) )
        labels_obj = sg.mergeObj( labels_obj, obj )

# List the closest approaches among boreholes
if True:
    for i in range(len(bores.xo)):
        for j in range(len(bores.xo)):
            dist=sg.shortestDistanceBetweenLineSegments(bores.xo[i],bores.xf[i], bores.xo[j],bores.xf[j])
            print("%10s %10s %.1f"%(bores.cmt[i],bores.cmt[j],dist))


if True:
    # Note that I put the drift in twice as a lazy way to pin the colormap if we don't have other objects
    objectList=[ labels_obj,  drift_obj,  drift_scan, rhyolite, drift_obj, monitoring_obj, injection_obj, production_obj]
    colorList= [ 0.75,        0.0,        0.0,        0.85,     0.2,       0.75,           0.5,           1.0]
    # Add in the fractures
    nJoints=len(J.keys())
    fct=1.0/(nJoints-1.0)
    iCnt=0
    if treatStochasticAsDeterministic:
        for key in sorted(J.keys()): # This will loop through all of the joint sets
        #if False:
            JS=J[key]
            NF=JS["nf"] # Handy link to the dictionary for the current natural fracture set
            #print 'len(NF["obj"])',len(NF["obj"])
            for i in range(len(NF["obj"])):
                objectList.append(NF["obj"][i])
                colorList.append(NF["stableMu"][i])
                #colorList.append(fct*iCnt) # color by joint set
            iCnt+=1
        
    #print colorList
    sg.writeVtk(objectList, [colorList],["Color"], vtkFile=vtkDir+"/Boreholes.vtk")

if True:
    print("Looping over fracture sets for visualization")
    allNFs={"obj":[],"stableMu":[],"criticalDelPp":[],"link":[]}
    for key in sorted(J.keys()): # This will loop through all of the joint sets
        # For each joint set we will create a separate dictionary that we will use to store all natural fracture information
        JS=J[key]
        print("JS",key)
        NF=JS["nf"] # Handy link to the dictionary for the current natural fracture set
        # Write the natural fractures with their fields
        sg.writeVtk(NF["obj"], [NF["stableMu"],NF["criticalDelPp"],NF["link"]], ["stableMu","criticalDeltaPp","link"], vtkFile=vtkDir+"/NatFracs_%s.vtk"%(key))
        # Append to the visualization for all fractures
        allNFs["obj"]+=NF["obj"]
        allNFs["stableMu"]+=NF["stableMu"]
        allNFs["criticalDelPp"]+=NF["criticalDelPp"]
        allNFs["link"]+=NF["link"]
    sg.writeVtk(allNFs["obj"], [allNFs["stableMu"],allNFs["criticalDelPp"],allNFs["link"]], ["stableMu","criticalDeltaPp","link"], vtkFile=vtkDir+"/NatFracs_All.vtk")

if True:
    # Plot more traditional stress diagrams for stableMu/slipTendancy and criticalDeltaPp
    sc.plotCritPpStableMu(Sh,SH,SV,ShAzimuthDeg,ShDipDeg,Pp,
                          muSlip,
                          nf_dip_dir_deg, nf_dip_deg,
                          critDelPpFile="Critical_Del_Pp_withFracs.png", stableMuFile='Stable_Mu_Slip_Tendancy_withFracs.png', mohrsFile="Mohrs_Circle_withFracs.png",)

print("Finished")
