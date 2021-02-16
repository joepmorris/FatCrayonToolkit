# Copyright 2018-2021 Lawrence Livermore National Security, LLC and other
# Fat Crayon Toolkit Project Developers. See the top-level COPYRIGHT file for details.
import pandas as pd
import Units as unit
import numpy as np
from scipy import interpolate
import math as math
import SimpleGeometry as sg
import StressCalc as sc
import random



def read_laser_pts(pts_file, keepFraction=0.01):
    np.random.seed(1)
    # Reads a laser scan PTS file and returns an array of points in a form that looks like a Panda dataframe
    # dt['x'][], dt['x'][], dt['x'][]
    dt = {'x':[], 'y':[], 'z':[]}
    fd = open(pts_file,'r')
    # Skip the first line, it is only the number of points:
    line = fd.readline()
    for line in fd:
        col = line.split()
        if random.random()<keepFraction:
            dt['x'].append(float(col[0]))
            dt['y'].append(float(col[1]))
            dt['z'].append(float(col[2]))
            #print dt; quit()
    return dt

# Read in a fracture pick file
class FracPick:
    def __init__(self,file,survey=None):
        dt=pd.read_excel(file, skiprows=1,
                         names=["depth_m","dip_dir_deg","dip_deg","ap_mm","type"],
                         usecols=range(5),
                         index_col=None,
                         na_values=['NA'])
        #print dt["depth_m"]
        # If requested we now need to read in the corresponding borehole to map from depth to coordinates
        if survey!=None:
            self.e=[]; self.n=[]; self.u=[]
        self.depth_ft=[]; self.dip_dir_deg=[]; self.dip_deg=[]; self.ap_mm=[]; self.type=[]
        for i in range(len(dt["depth_m"])):
            self.depth_ft.append( float(dt["depth_m"][i])/unit.ft )
            self.dip_dir_deg.append( float(dt["dip_dir_deg"][i]) )
            self.dip_deg.append( float(dt["dip_deg"][i]) )
            self.ap_mm.append( float(dt["ap_mm"][i]) )
            self.type.append( int(dt["type"][i]) )
            if survey!=None:
                #print float(dt["depth_m"][i])/unit.ft
                enu=survey.enu_at_md(float(dt["depth_m"][i])/unit.ft)
                self.e.append(enu[0])
                self.n.append(enu[1])
                self.u.append(enu[2])
        if survey!=None:
            self.e=np.asarray(self.e)
            self.n=np.asarray(self.n)
            self.u=np.asarray(self.u)
        self.depth_ft=np.asarray(self.depth_ft)
        self.dip_dir_deg=np.asarray(self.dip_dir_deg)
        self.dip_deg=np.asarray(self.dip_deg)
        self.ap_mm=np.asarray(self.ap_mm)
        self.type=np.asarray(self.type)
        #print self.depth_ft
        # Azimuth increases in clockwise direction
        # Strike has dip direction on right
        # So a dip direction at 90 degrees (east) would have a strike direction of 0 (north)
        self.strike_deg=(self.dip_dir_deg-90)
        self.strikeRad=(self.dip_dir_deg-90)*(math.pi/180.0)
        self.dipRad=self.dip_deg*(math.pi/180.0)

    def write(self):
        print "%10s\t%10s\t%10s\t%10s\t%10s"%("i","depth_ft","dip_dir_deg","strike_deg","dip_deg")
        for i in range(len(self.depth_ft)):
            print "%10d\t%10g\t%10g\t%10g\t%10g"%(i,self.depth_ft[i],self.dip_dir_deg[i],self.strike_deg[i],self.dip_deg[i])
        
    # This does not belong here. We should only be providing helpful functions for reading data
    def appendVis(self, obj, r, h, sigG, Pp):
        # Append visualization of the natural fractures to the given object
        for i in range(len(self.e)):
            #obj=sg.mergeObj( obj, sg.HF( r=r, x0=[self.e[i],self.n[i],self.u[i]], strikeRad=self.strikeRad[i], dipRad=self.dipRad[i], h=h) )
            # Size reflects aperture
            #obj=sg.mergeObj( obj, sg.HF( r=0.1+0.1*self.ap_mm[i], x0=[self.e[i],self.n[i],self.u[i]], strikeRad=self.strikeRad[i], dipRad=self.dipRad[i], h=h) )
            # Size reflects slip tendency
            # Get normal to the fracture
            nrmG=sc.nrmFromDipDirDipAngDeg(dip_direction=self.dip_dir_deg[i],dip_angle=self.dip_deg[i])
            criticalDelPp, stableMu = sc.criticalDelPp_stableMu(nrmG, sigG, Pp)
            obj=sg.mergeObj( obj, sg.HF( r=0.1+0.1*self.ap_mm[i], x0=[self.e[i],self.n[i],self.u[i]], strikeRad=self.strikeRad[i], dipRad=self.dipRad[i], h=h) )
            #obj=sg.mergeObj( obj, sg.HF( r=r, x0=[self.e[i],self.n[i],self.u[i]], strikeRad=self.strikeRad[i], dipRad=self.dipRad[i], h=0.1+0.1*self.ap_mm[i] ) )
            #obj=sg.mergeObj( obj, sg.HF( r=r, x0=[self.e[i],self.n[i],self.u[i]], strikeRad=self.strikeRad[i], dipRad=self.dipRad[i], h=h*min(0.5,0*self.ap_mm[i])) )
        #print obj
        return obj
        
class Survey:
    def __init__(self,surveyFile, name="gyro_filename", correction=[0.0,0.0,0.0], format=0, match="", skiprows=0):
        print "Reading",surveyFile
        # Sometimes the survey is systematically off
        self.correction=np.array(correction)
        self.name=name
        if (format==0):
            dt=pd.read_excel(surveyFile, skiprows=skiprows,
                         #          A       B       C       D     E     F      G         H         I      J          K
                         names=["station","Dip","Azimuth","GHS","GTF","DLS","Easting","Northing","TVD","UpDown","LeftRight",
                                #    L           M            N              O                     P
                                "Shortfall","GravField","Drift check","Accumulated roll","Accumulated roll (rev)",
                                #    Q              R                        S                   T           U
                                "Roll to station","Roll to station (rev)","Time to station","Temperature","Battery",
                                #    V              W               X               Y            Z
                                "Station acc.","Station rate","Station Quality","Max. rate","Motion Quality"],
                         index_col=None,
                         na_values=['NA'])
        elif (format==1):
            dt=pd.read_excel(surveyFile, skiprows=skiprows,
                         #          A       B       C       D          E         F
                         names=["station","Dip","Azimuth","Easting","Northing","TVD"],
                         index_col=None,
                         na_values=['NA'])
        elif (format==2):
            dt=pd.read_excel(surveyFile, skiprows=skiprows,
                         #          A       B       C       D       E         F        G
                         names=["station","Dip","Azimuth","DLS","Easting","Northing","TVD"],
                         index_col=None,
                         na_values=['NA'])
        elif (format=="match_dxyz"):
            # We will scan through the file and look for lines that match the string we were given
            # For this format, we assume
            #  0       1         2     3 4 5
            #HoleID,Depth(ft),Depth(m),x,y,z
            x=[]; y=[]; z=[]; md=[]
            with open(surveyFile) as fd:
                for line in fd:
                    col=line.split(',')
                    if col[0]==match:
                        x.append(float(col[3]))
                        y.append(float(col[4]))
                        z.append(float(col[5]))
                        md.append(float(col[1]))
            # Now we pack the data we collected into a format that matches what we've gotten from other formats
            dt={"Easting":x, "Northing":y, "TVD":z, "station":md}
            #print dt; quit()
        else:
            print "failed"
            quit()
        #print df.keys()
        #print dt["Easting"]
        self.e=dt["Easting"]
        self.n=dt["Northing"]
        self.u=dt["TVD"]
        self.md=dt["station"]
        #print df[""]
        # The list of slices corresponding to multiple runs (typically going in and out)
        self.i0=[]; self.i1=[]
        i=0
        while(True):
            # Search for a number
            while(True):
                try:
                    val=float(self.e[i])
                except:
                    i+=1
                    continue
                self.i0.append(i)
                break
            # Search for the end of numbers
            while(True):
                try:
                    val=float(self.e[i])
                except:
                    break
                if math.isnan(val):
                    #print i,val,"is Nan?"
                    break
                i+=1
                if (i>=len(self.e)):
                    break
            self.i1.append(i)
            i+=1
            if (i>=len(self.e)):
                break
        # We now have a sequence of slices into the data we were given
        #print "i0",self.i0
        #print "i1",self.i1
        #for i in range(len(self.i0)):
        #    print self.i0[i],self.i1[i]
        #    print self.e[self.i0[i]:self.i1[i]]

        xx=[]; yy=[]
        for ii in range(self.i0[0],self.i1[0]):
            xx.append(float(self.md[ii]))
            yy.append(float(ii))
        xx=np.asarray(xx); yy=np.asarray(yy)
        nL=len(xx)
        #for i in range(nL):
        #    print 'xx,yy',xx[i],yy[i]
        #print "about to interp",len(xx),len(yy)
        self.interpMdToI=interpolate.interp1d(xx,yy)
        self.interpMdToI_min=xx[0]
        self.interpMdToI_max=xx[-1]
        #for md in [45.0,46.0,47.0,48.0,49.0,50.0]:
        #    print md,self.interpMdToI(md),self.enu_at_md(md)
        #quit()

    def enu_at_md(self,md):
        if md>self.interpMdToI_max:
            print md,">",self.interpMdToI_max,"in",self.name,"setting to end of survey"
            md=0.999*self.interpMdToI_max
        ii=self.interpMdToI(md)
        i=int(ii)
        fct=ii-i
        return (1.0-fct)*self.x(i)+fct*self.x(i+1)
        
    def x(self, j):
        return np.asarray( [float(self.e[j]),float(self.n[j]),float(self.u[j])] )+self.correction

    def collar(self):
        return self.x( self.i0[0] )

    def compare(self, x0,v,survey_index=0):
        # Compare the prescribed line with our own
        fd=open(self.name+"_compare.txt",'w')
        print "md_feet dist_feet"
        for j in range(self.i0[survey_index],self.i1[survey_index]):
            x=self.x(j)
            dv=x-x0 # Vector to point x from the prescribed collar (x0)
            d=np.sqrt(np.dot(dv,dv)) # Distance from collar
            s=np.dot(dv,v) # Projected distance along the design borehole
            m=np.sqrt( d*d - s*s ) # Perpendicular distance back to the design borehole
            # Note: s is not measured depth. Measured depth will be farther than s
            fd.write("%g %g\n"%(self.md[j],m))
        fd.close()
    
    def appendVis(self, obj, r, survey_indices=[]):
        # Append visualization of the survey to the given object
        if survey_indices==[]:
            survey_indices=range(len(self.i0))
        for i in survey_indices:
            for j in range(self.i0[i],self.i1[i]-1):
                #x0=np.asarray( [float(self.e[j]),float(self.n[j]),float(self.u[j])] )
                x0=self.x(j)
                #x1=np.asarray( [float(self.e[j+1]),float(self.n[j+1]),float(self.u[j+1])] )
                x1=self.x(j+1)
                obj=sg.mergeObj( obj, sg.cylObj( x0=x0, x1=x1, r=r ) )
        return obj

class NatFracsSurvey:
    def __init__(self,surveyFile):
        dt=pd.read_excel(surveyFile, skiprows=0,
                         names=["ID","Easting","Northing","Drifts","Distance","DipDir","Inclination","Comment","Strike"],
                         usecols=range(9),
                         index_col=None,
                         na_values=['NA'])
        self.e=[]; self.n=[]; self.strikeRad=[]; self.dipRad=[]
        for i in range(len(dt["Distance"])):
            #print dt["Distance"][i]
            self.e.append( float(dt["Easting"][i]) )
            self.n.append( float(dt["Northing"][i]) )
            self.strikeRad.append( float(dt["Strike"][i])*math.pi/180.0 )
            self.dipRad.append( float(dt["Inclination"][i])*math.pi/180.0 )
        
    def appendVis(self, obj, r, u, h):
        # Append visualization of the natural fractures to the given object
        for i in range(len(self.e)):
            obj=sg.mergeObj( obj, sg.HF( r=r, x0=[self.e[i],self.n[i],u], strikeRad=self.strikeRad[i], dipRad=self.dipRad[i], h=h) )
        #print obj
        return obj
        


# Yet another format of fracture pick file
# This one has the format:
#   0       1       2         3            4             5             6      7        8       9             10   11   12
#wellID,Depth[ft],Pole Trend,Pole  Plunge,Dip Vector Tr,Dip Vector Pl,Strike,DipAngle,FracSet,EquivRadius[m],x_ft,y_ft,z_ft
# I believe that Dip Vector Tr,Dip Vector Pl are the ones we are after
class CDFN_Summary:
    def __init__(self,file,match):
        self.e=[]; self.n=[]; self.u=[]
        self.dip_dir_deg=[]; self.dip_deg=[]
        with open(file) as fd:
            for line in fd:
                col=line.split(',')
                #print col; quit()
                if col[0]==match:
                    self.e.append(float(col[10]))
                    self.n.append(float(col[11]))
                    self.u.append(float(col[12]))
                    self.dip_dir_deg.append(float(col[4]))
                    self.dip_deg.append(float(col[5]))
        self.e=np.asarray(self.e)
        self.n=np.asarray(self.n)
        self.u=np.asarray(self.u)
        self.dip_dir_deg=np.asarray(self.dip_dir_deg)
        self.dip_deg=np.asarray(self.dip_deg)
        # Azimuth increases in clockwise direction
        # Strike has dip direction on right
        # So a dip direction at 90 degrees (east) would have a strike direction of 0 (north)
        self.strike_deg=(self.dip_dir_deg-90)
        self.strikeRad=(self.dip_dir_deg-90)*(math.pi/180.0)
        self.dipRad=self.dip_deg*(math.pi/180.0)


# This will attempt to be smart and read column headings to get what we want
class DFN_CSV:
    def __init__(self,surveyFile,match):
        dt=pd.read_csv(surveyFile)
        mask = np.asarray([dt["Well_ID"]==match])[0]
        self.e=np.asarray(dt["x"])[mask]
        self.n=np.asarray(dt["y"])[mask]
        self.u=np.asarray(dt["z"])[mask]
        # Azimuth increases in clockwise direction
        # Strike has dip direction on right
        # So a dip direction at 90 degrees (east) would have a strike direction of 0 (north)
        if True:
            # According to Hari: Azimuth in these files is the dip direction. 
            print('As per Hari, assuming azimuth is dip direction')
            self.dip_dir_deg=np.asarray(dt["Azimuth_deg"])[mask]
            self.strike_deg=self.dip_dir_deg-90
        else:
            print('Assuming azimuth is dip direction')
            self.strike_deg=np.asarray(dt["Azimuth_deg"])[mask]
            self.dip_dir_deg=self.strike_deg+90

        self.dip_deg=np.asarray(dt["Dip_deg"][mask])
        self.strikeRad=(self.dip_dir_deg-90)*(math.pi/180.0)
        self.dipRad=self.dip_deg*(math.pi/180.0)
        
class DFN_Arrays:
    def __init__(self,x,y,z,dip_dir_deg,dip_deg,r=None):
        self.e=np.asarray(x)
        self.n=np.asarray(y)
        self.u=np.asarray(z)
        if r is not None:
            self.r=np.asarray(r)
        # Azimuth increases in clockwise direction
        # Strike has dip direction on right
        # So a dip direction at 90 degrees (east) would have a strike direction of 0 (north)
        self.dip_dir_deg=np.asarray(dip_dir_deg)
        self.strike_deg=self.dip_dir_deg-90
        self.dip_deg=np.asarray(dip_deg)
        self.strikeRad=(self.dip_dir_deg-90)*(math.pi/180.0)
        self.dipRad=self.dip_deg*(math.pi/180.0)
    
class ProLayout:
    def __init__(self,surveyFile, name="gyro_filename", correction=[0.0,0.0,0.0], format=0, skiprows=0):
        print "Reading",surveyFile
        # Sometimes the survey is systematically off
        self.correction=np.array(correction)
        self.name=name
        if (format==0):
            dt=pd.read_excel(surveyFile, skiprows=skiprows,
                         #          A       B       C       D     E     F      G         H         I
                         names=["BoreType","Length","Declination","Azimuth","Easting","Northing","TVD","Length_m","ID"],
                         index_col=None,
                         na_values=['NA'])
            self.e=np.asarray(dt["Easting"],dtype=float)
            self.n=np.asarray(dt["Northing"],dtype=float)
            self.u=np.asarray(dt["TVD"],dtype=float)
            self.y=dt["BoreType"]
            self.l=np.asarray(dt["Length"],dtype=float)
            self.d=np.asarray(dt["Declination"],dtype=float)
        elif (format=='csv'):
            dt=pd.read_csv(surveyFile, skiprows=skiprows)
            #print dt.head()
            print dt
            #print 'read in',dt["comment"]
            self.e=np.asarray(dt["Easting_ft"],dtype=float)
            self.n=np.asarray(dt["Northing_ft"],dtype=float)
            self.u=np.asarray(dt["TVD_ft"],dtype=float)
            self.y=dt["BoreType"]
            self.l=np.asarray(dt["Length_ft"],dtype=float)
            self.d=np.asarray(dt["Declination_deg"],dtype=float)
            try:
                self.d_err=np.asarray(dt["Declination_Error"],dtype=float)
            except:
                self.d_err=0.0*self.d
            self.a=np.asarray(dt["Azimuth_deg"],dtype=float)
            try:
                self.a_err=np.asarray(dt["Azimuth_Error"],dtype=float)
            except:
                self.a_err=0.0*self.d
            self.cmt=dt["comment"]
        else:
            print "failed"
            quit()

        # If the Declination or Azimuth errors are non-zero, we will introduce additional error boreholes
        nBores = len(self.e)
        for i in range(nBores):
            if self.d_err[i]>0.0 or self.a_err[i]>0.0:
                # We have been asked to add in some uncertain orientations
                for dd,aa in [ [-1.0,0.0], [1.0,0.0], [0.0,-1.0], [0.0,1.0] ]:
                    self.e=np.append(self.e,self.e[i])
                    self.n=np.append(self.n,[self.n[i]])
                    self.u=np.append(self.u,[self.u[i]])
                    self.y=np.append(self.y,[self.y[i]])
                    self.l=np.append(self.l,[self.l[i]])
                    self.d=np.append(self.d,[self.d[i]+dd*self.d_err[i]])
                    self.d_err=np.append(self.d_err,[self.d_err[i]])
                    self.a=np.append(self.a,[self.a[i]+aa*self.a_err[i]])
                    self.a_err=np.append(self.a_err,[self.a_err[i]])
                    self.cmt=np.append(self.cmt,[self.cmt[i]])
        #Collar
        self.xo = np.transpose(np.asarray([self.e,self.n,self.u]))+self.correction
        
        #Toe
        ff = []
        self.iInjection=[]
        self.iProduction=[]
        for j in range(0,len(self.e)):
            azn = float(self.a[j])
            dec = float(self.d[j])
            vN = np.asarray([math.sin(azn*math.pi/180.0)*math.cos(-dec*math.pi/180.0),
                         math.cos(azn*math.pi/180.0)*math.cos(-dec*math.pi/180.0),math.sin(-dec*math.pi/180.0)])
            #ff += [self.xo[j] + float(self.l[j])*vN + self.correction]
            ff += [self.xo[j] + float(self.l[j])*vN] # Correction is already in xo
            if self.y[j] == 'I':
                self.iInjection.append(j)
            elif self.y[j] == 'P':
                self.iProduction.append(j)
        self.xf = np.asarray(ff)
    
    def appendVis(self, injection_obj, production_obj, monitoring_obj, r): #Add layout to geometry
        #print '/n***** running appendVis'
        for j in range(0,len(self.e)):
            if self.y[j] == 'I':
                injection_obj = sg.mergeObj(injection_obj,sg.cylObj(x0=self.xo[j],x1=self.xf[j],r=r))
            elif self.y[j] == 'P':
                production_obj = sg.mergeObj(production_obj,sg.cylObj(x0=self.xo[j],x1=self.xf[j],r=r))
            elif self.y[j] == 'M':
                monitoring_obj = sg.mergeObj(monitoring_obj,sg.cylObj(x0=self.xo[j],x1=self.xf[j],r=r))
            else:
                print 'error: appendVis geometry type unknown'
            #print 'j = %i' %(j)
        return injection_obj, production_obj, monitoring_obj
