# Copyright 2018-2021 Lawrence Livermore National Security, LLC and other
# Fat Crayon Toolkit Project Developers. See the top-level COPYRIGHT file for details.
""" Classes and routines for generating 3D objects
"""
import math
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from copy import deepcopy
import sys
import random
import re

# If you have placed the other modules in another directory, this can be useful to add that location to the path
#import os, sys; sys.path.append(os.path.dirname(__file__)); import Units as unit

def dipDirectionAndDipAng(tangent):
    # Given the tangent to a curve, convert into dip direction and dip
    e=tangent[0]
    n=tangent[1]
    up=tangent[2]
    # If we rotate compass to align with math coords:
    #        W
    #        |
    #    S-------N
    #        |
    #        E
    x=n
    y=-e
    thetaMath=math.atan2(y,x)
    thetaCompass=-thetaMath*180.0/math.pi
    dipDirection=thetaCompass
    # Dip angle is the amount we are dipping from horizontal
    # We chose orientation such that up is -ve
    dipAngle=math.atan2( -up, math.sqrt( e*e + n*n ) )*180.0/math.pi
    return dipDirection,dipAngle
    
def dipToStrikeDeg(dip_dir_deg):
    # Definitions published by Wikipedia: https://en.wikipedia.org/wiki/Strike_and_dip
    # One technique is to always take the strike so the dip is 90 deg to the right of the strike, in which case the redundant letter following the dip angle is omitted (right hand rule, or RHR).
    #strike_rad=dip_dir_radians-0.5*np.pi
    strike_deg=dip_dir_deg-90.0
    return strike_deg

def strikeToDipDeg(strike_deg):
    # Definitions published by Wikipedia: https://en.wikipedia.org/wiki/Strike_and_dip
    # One technique is to always take the strike so the dip is 90 deg to the right of the strike, in which case the redundant letter following the dip angle is omitted (right hand rule, or RHR).
    #strike_rad=dip_dir_radians-0.5*np.pi
    #strike_deg=dip_dir_deg-90.0
    dip_dir_deg=strike_deg+90.0
    return dip_dir_deg

def degToRad(deg):
    return deg*np.pi/180.0

def radToDeg(rad):
    return rad*180.0/np.pi

# Return a tangent normal given azimuth and declination in degrees
def vectFromAzDecDeg(azDeg, decDeg):
    v=np.asarray([0.0,1.0,0.0]) # Due North
    v=rotatePoints( [v], np.asarray([1,0,0]), -degToRad(decDeg) )[0] # Rotate decDeg degrees sub-horizontal
    v=rotatePoints( [v], np.asarray([0,0,1]), -degToRad(azDeg) )[0] # Rotate azDeg degrees clockwise of North
    return v

# Return azimuth and dec from a tangent normal
def azDecDegFromVect(v):
    # math.atan2(y, x): Return atan(y / x), in radians. The result is between -pi and pi
    return (
        (90.0-radToDeg(math.atan2(v[1],v[0])))%360.0, # Azimuth has different sign convention than math convention
        -radToDeg(math.atan2(v[2],math.sqrt(v[0]*v[0]+v[1]*v[1])))
        )

# Return a normalized vector
def normalize(v):
    return v/math.sqrt(np.dot(v, v))

def writeStlObject(points,simplices,fd):
    for simplex in simplices:
        # I'm not sure we need to calculate a normal...
        fd.write("facet normal 0.0 0.0 0.0\n")
        fd.write("outer loop\n")
        for iPt in simplex:
            #print iPt,simplex
            fd.write("vertex %g %g %g\n"%(points[iPt][0],points[iPt][1],points[iPt][2]))
        fd.write("endloop\n")
        fd.write("endfacet\n")
    
def writeStlFile(points,simplices,stlFile,name="stlObject"):
    fd=open(stlFile,'w')
    fd.write("solid %s\n"%(name))
    writeStlObject(points,simplices,fd)
    fd.write("endsolid\n");

def writeObjectsStlFile(objects,stlFile,name="stlObject"):
    fd=open(stlFile,'w')
    fd.write("solid %s\n"%(name))
    #for object in objects:
    (points,simplices)=objects
    writeStlObject(points,simplices,fd)
    fd.write("endsolid\n");

def writeVtk(objectListIn,scalarsIn,scalarNames,vtkFile,name="vtkObjects"):
    
    # Remove empty object lists

    #print 'len(objectListIn)',len(objectListIn),'len(scalarsIn)',len(scalarsIn),'len(scalarsIn[0])',len(scalarsIn[0])
    if False:
        # I don't get it, but this seems to be misbehaving now:
        # In retrospect, it isn't clear it ever worked for the case where we have more than one scalar!
        objectList=[objectListIn[i] for i in range(len(objectListIn)) if objectListIn[i] is not None]
        scalars=[[scalarsIn[0][i] for i in range(len(objectListIn)) if objectListIn[i] is not None]]
    if False:
        # This is for a more general case with multiple scalars
        # But it doesn't seem to work
        objectList=[objectListIn[i] for i in range(len(objectListIn)) if objectListIn[i] is not None]
        #scalars=[[scalarsIn[:][i]] for i in range(len(objectListIn)) if objectListIn[i] is not None]
        scalars=[scalarsIn[:][i] for i in range(len(objectListIn)) if objectListIn[i] is not None]
    if True:
        # This works, but is not Pythonic
        objectList=[]
        scalars=[]
        for i in range(len(scalarsIn)):
            scalars.append([])
        for i in range(len(objectListIn)):
            if objectListIn[i] is not None:
                objectList.append(objectListIn[i])
                for iS in range(len(scalarsIn)):
                    scalars[iS].append(scalarsIn[iS][i])
        
    #print objectList
    #print scalars
    #print 'len(objectList)',len(objectList),'len(scalars)',len(scalars)
    
    fd=open(vtkFile,'w')

    nPtsObj=[]
    nPts=0
    nTri=0
    nObj=len(objectList)
    for pts,simps in (objectList):
        nPtsObj.append(len(pts))
        nPts+=len(pts)
        nTri+=len(simps)
    nShift=[0]*nObj
    for iShift in range(nObj-1):
        nShift[iShift+1]=nShift[iShift]+nPtsObj[iShift]
    
    fd.write("# vtk DataFile Version 2.0\n")
    fd.write("%s\n"%(name))
    fd.write("ASCII\n")
    fd.write("DATASET UNSTRUCTURED_GRID\n")
    fd.write("POINTS %d float\n"%(nPts))
    for pts,simps in (objectList):
        for pt in (pts):
            fd.write("%g %g %g\n"%(pt[0],pt[1],pt[2]))
    
    fd.write("CELLS %d %d\n"%(nTri,(1+3)*nTri))
    iObj=0
    #col=[]
    for pts,simps in (objectList):
        for tri in (simps):
            fd.write("3 %d %d %d\n"%(tri[0]+nShift[iObj],tri[1]+nShift[iObj],tri[2]+nShift[iObj]))
            #col.append(colorList[iObj])
        iObj+=1
            
    fd.write("CELL_TYPES %d\n"%(nTri))
    for i in range(nTri):
        fd.write("5 ") # http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf (see Fig. 2)
        if (i%10==9):
            fd.write("\n")
    fd.write("\n")

    fd.write("CELL_DATA %d\n"%(nTri))


    # Repeat as many of these as you want to define data on the tris
    #for colorList, scalarName in scalars,scalarNames:
    #print "check",len(scalars),scalarNames
    for iCol in range(len(scalars)):
        colorList=scalars[iCol]; scalarName=scalarNames[iCol]
        fd.write("SCALARS "+scalarName+" float 1\n")
        fd.write("LOOKUP_TABLE default\n")
        iObj=0
        i=0
        for pts,simps in (objectList):
            for tri in (simps):
                fd.write("%g "%(colorList[iObj])); i+=1
                if (i%10==9):
                    fd.write("\n")
            iObj+=1
    fd.write("\n")
    fd.close()
    
    
def simplicesFromPoints(points):
    hull=ConvexHull(points)
    return hull.simplices

def convexFromPoints(points):
    return ( points, simplicesFromPoints(points) )

def delaunayFromPoints(points):
    return ( points, Delaunay(points).simplices )

# A non-object
emptyObject=None

# Merging two objects requires a shift in the indices
def mergeObj(obj1, obj2):
    if (obj1 is None):
        return obj2
    if (obj2 is None):
        return obj1
    if (obj1==emptyObject):
        return obj2
    if (obj2==emptyObject):
        return obj1
    return (
        np.vstack( (obj1[0],obj2[0]) ),
        np.vstack( (obj1[1],obj2[1]+len(obj1[0])) )
    )

def mergeObjects(objects):
    nObj=len(objects)
    merged=np.asarray(deepcopy(objects[0]))
    nShift=0
    for i in range(nObj-1):
        #print i
        nShift+=len(objects[i][0])
        merged[0]=np.vstack( (merged[0],objects[i+1][0]) )
        merged[1]=np.vstack( (merged[1],objects[i+1][1]+nShift) )
    return merged
    
# Some useful objects
unitCubePts=np.asarray([
    [-0.5,-0.5,-0.5],
    [ 0.5,-0.5,-0.5],
    [-0.5, 0.5,-0.5],
    [ 0.5, 0.5,-0.5],
    [-0.5,-0.5, 0.5],
    [ 0.5,-0.5, 0.5],
    [-0.5, 0.5, 0.5],
    [ 0.5, 0.5, 0.5]
    ])
Cube=convexFromPoints(unitCubePts)
unitWedgePts=np.asarray([
    [-0.5,-0.5,-0.5],
    [ 0.5,-0.5,-0.5],
    [ 0.0, 0.5,-0.5],
    [-0.5,-0.5, 0.5],
    [ 0.5,-0.5, 0.5],
    [ 0.0, 0.5, 0.5]
    ])
unitWedge=convexFromPoints(unitWedgePts)

def diskObj(r, h, n=50):
    dTh=2*math.pi/n
    pts=[]
    for i in range(n):
        x=r*math.cos(i*dTh); y=r*math.sin(i*dTh)
        pts.append( [x,y,-0.5*h] )
        pts.append( [x,y, 0.5*h] )
    pts=np.asarray(pts)
    return convexFromPoints(pts)

# This was published on: https://en.wikipedia.org/wiki/Regular_dodecahedron
# Golden ratio
gr=(1.0+math.sqrt(5.0))/2.0
radiusOneSpherePts=np.asarray([
    [-1,-1,-1],[ 1,-1,-1], [-1, 1,-1],[ 1, 1,-1], [-1,-1, 1],[ 1,-1, 1], [-1, 1, 1],[ 1, 1, 1],
    [0,-1/gr,-gr],[0, 1/gr,-gr],[0,-1/gr, gr],[0, 1/gr, gr],
    [-1/gr,-gr,0],[ 1/gr,-gr,0],[-1/gr, gr,0],[ 1/gr, gr,0],
    [-gr,0,-1/gr],[-gr,0, 1/gr],[ gr,0,-1/gr],[ gr,0, 1/gr]
    ])
radiusOneSphereObj=convexFromPoints(radiusOneSpherePts)

def randSpherePtsFromGaussians(n,rad):
    np.random.seed(1)
    pts=[]
    for i in range(n):
        u = np.random.normal(0,1)
        v = np.random.normal(0,1)
        w = np.random.normal(0,1)
        d= math.sqrt(u*u+v*v+w*w)
        pts.append( rad*np.asarray([u,v,w])/d )
    return pts


def cylObj(x0, x1, r, n=10, lengthSum=None):
    sphere0=(r*radiusOneSpherePts)
    sphere0[:,0]+=x0[0]; sphere0[:,1]+=x0[1]; sphere0[:,2]+=x0[2];
    sphere1=(r*radiusOneSpherePts)
    sphere1[:,0]+=x1[0]; sphere1[:,1]+=x1[1]; sphere1[:,2]+=x1[2];
    pts=np.vstack( (sphere0, sphere1) )
    #print lengthSum
    #if (lengthSum != None):
    try:
        lengthSum[0]+=np.sqrt( np.dot((x1-x0),(x1-x0)) )
    except:
        pass
    #print lengthSum
    return convexFromPoints(pts)

# Set up a unit arrow pointing in y-direction
pts1=deepcopy(unitCubePts); pts1[:,1]-=0.5
pts2=deepcopy(unitWedgePts); pts2[:,0]*=2.0; pts2[:,1]+=0.5
unitArrow1=convexFromPoints(pts1)
unitArrow2=convexFromPoints(pts2)
unitArrowY=mergeObj(unitArrow1,unitArrow2)

def extrudePoints(points, disp):
    """
    Return a list of points including the initial points and extruded end
    """
    farEnd=deepcopy(points)
    farEnd[:,0]+=disp[0]
    farEnd[:,1]+=disp[1]
    farEnd[:,2]+=disp[2]
    return np.vstack( (points,farEnd) )

def transObj(object, disp):
    """
    Translate an object
    """
    return (object[0]+disp,object[1])

def scaleObj(object, scale):
    """
    Scale an object
    """
    return (object[0]*scale,object[1])


def getPoints(object):
    """
    Return the list of points/vertices
    """
    return object[0]

def getNPoly(object):
    """
    Return the number polygons in the object
    """
    return len(object[1])

def getPolyPoints(object, i):
    """
    Return the list of points for polygon i
    """
    return object[0][object[1][i]]

def getPolyNormal(object, x0, i):
    """
    Return the normal to polygon i that points away from x0
    """
    pts=getPolyPoints(object,i)
    # Note that this normal could be badly behaved if aVec and bVec are close to parallel
    aVec=pts[2]-pts[0]
    bVec=pts[1]-pts[0]
    nVec=np.cross(aVec,bVec)
    nVec=nVec/np.linalg.norm(nVec)
    # Check if our normal is pointing away from x0
    #print 'nVec',nVec
    #print 'pts[0]',pts[0]
    #print 'x0',x0
    if np.dot( nVec, pts[0]-x0 ) > 0.0:
        return nVec
    else:
        return -1.0*nVec

def getPolyArea(object, i):
    """
    Return the area of the polygon i
    """
    pts = getPolyPoints(object, i)
    area=0.0
    for j in range(1,len(pts)-1):
        #print 'j',j
        vtmp = np.cross(pts[j]-pts[0],pts[j+1]-pts[0])
        area += 0.5*np.sqrt(np.dot(vtmp,vtmp))
    return area


# http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
def rotationMatrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def rotatePoints(points, axis, theta):
    rot=rotationMatrix(axis,theta)
    #return rot*points
    return np.transpose( np.dot(rot, np.transpose( points ) ) )

def rotateTensor(tensor, axis, theta):
    #http://www.continuummechanics.org/stressxforms.html
    rot=rotationMatrix(axis,theta)
    return np.dot(rot,np.dot(tensor,np.transpose(rot)))
    
def rotateObj(object, axis, theta):
    rot=rotationMatrix(axis,theta)
    return ( np.transpose( np.dot(rot, np.transpose(object[0])) ) , object[1])

# This was published by http://geomalgorithms.com/a05-_intersect-1.html
def intersectionOfLineAndPlane(lineX,lineS,planeX,planeN):
    V0=np.asarray(planeX)
    n=np.asarray(planeN)
    P0=np.asarray(lineX)
    u=np.asarray(lineS)
    sI=( np.dot( n, (V0-P0) ) )/( np.dot( n,u ) )
    return P0+sI*u
def distOfIntersectionOfLineAndPlane(lineX,lineS,planeX,planeN):
    V0=np.asarray(planeX)
    n=np.asarray(planeN)
    P0=np.asarray(lineX)
    u=np.asarray(lineS)
    sI=( np.dot( n, (V0-P0) ) )/( np.dot( n,u ) )
    return sI,P0+sI*u

def shortestDistanceBetweenLineSegments( xio,xif, xjo,xjf ):
    # Calculate tangents to the line segments
    p1=xio; p2=xif
    p3=xjo; p4=xjf
    # The Python code in this function is based upon C++ code developed by Dan Sunday.
    # The original C++ code had the following request:
    # // Copyright 2001 softSurfer, 2012 Dan Sunday
    # // This code may be freely used and modified for any purpose
    # // providing that this copyright notice is included with it.
    # // SoftSurfer makes no warranty for this code, and cannot be held
    # // liable for any real or imagined damage resulting from its use.
    # // Users of this code must verify correctness for their application.
    
    u = p1 - p2;
    v = p3 - p4;
    w = p2 - p4;

    a = np.dot(u,u);
    b = np.dot(u,v);
    c = np.dot(v,v);
    d = np.dot(u,w);
    e = np.dot(v,w);
    D = a*c - b*b;
    sD = D;
    tD = D;

    SMALL_NUM = 0.00000001;
    
    # compute the line parameters of the two closest points
    if (D < SMALL_NUM):  # the lines are almost parallel
        sN = 0.0;       # force using point P0 on segment S1
        sD = 1.0;       # to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    else:                # get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0):   # sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        elif (sN > sD):# sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
    
    if (tN < 0.0):            # tc < 0 => the t=0 edge is visible
        tN = 0.0;
        # recompute sc for this edge
        if (-d < 0.0):
            sN = 0.0;
        elif (-d > a):
            sN = sD;
        else:
            sN = -d;
            sD = a;
    elif (tN > tD):       # tc > 1 => the t=1 edge is visible
        tN = tD;
        # recompute sc for this edge
        if ((-d + b) < 0.0):
            sN = 0;
        elif ((-d + b) > a):
            sN = sD;
        else:
            sN = (-d + b);
            sD = a;
    
    # finally do the division to get sc and tc
    if(abs(sN) < SMALL_NUM):
        sc = 0.0;
    else:
        sc = sN / sD;
    
    if(abs(tN) < SMALL_NUM):
        tc = 0.0;
    else:
        tc = tN / tD;
    
    # get the difference of the two closest points
    dP = w + (sc * u) - (tc * v);
    distance = np.linalg.norm(dP);

    return distance


# Generate a convex hull from points - Not good in general because drifts are not always convex
def pointCloudToConvexPolyhedron(drift_scan, dt, keepFraction=0.01):
    np.random.seed(1)
    nPoints=len(dt['x'])
    pts = []
    for i in range(nPoints):
        #print "%.1f"%((100.0*i)/nPoints)
        if (i%100==0):
            sys.stdout.write("Scan progress %d%%   \r" % ((100.0*i)/nPoints) )
            sys.stdout.flush()
        if random.random()<keepFraction:
            pts.append( [dt['x'][i],dt['y'][i],dt['z'][i]] )
    drift_scan=mergeObj( drift_scan, convexFromPoints( pts ) )
    return drift_scan

# Generate a Delaunay triangular mesh from points - Not good on its own because it will not handle concavity
# However, we can prune the larger triangles to recover a resonable representation of a complex tunnel
def pointCloudToDelaunay(drift_scan, dt, keepFraction=0.01):
    np.random.seed(1)
    nPoints=len(dt['x'])
    pts = []
    for i in range(nPoints):
        #print "%.1f"%((100.0*i)/nPoints)
        if (i%100==0):
            sys.stdout.write("Scan progress %d%%   \r" % ((100.0*i)/nPoints) )
            sys.stdout.flush()
        if random.random()<keepFraction:
            pts.append( [dt['x'][i],dt['y'][i],dt['z'][i]] )
    drift_scan=mergeObj( drift_scan, delaunayFromPoints( pts ) )
    return drift_scan

# Remove large triangles from an object - we assume we are given a bunch of triangles
# We drop any triangle with at least one side larger than L
def pruneLargeTriangles(obj, L):
    pts, simps_in = obj
    simps_out = []
    nTri = len(simps_in)
    for i in range(nTri):
        if (i%100==0):
            sys.stdout.write("Pruning progress %d%%   \r" % ((100.0*i)/nTri) )
            sys.stdout.flush()
        tri = simps_in[i]
        if np.linalg.norm( np.asarray(pts[tri[1]])-np.asarray(pts[tri[0]]) ) > L:
            continue
        if np.linalg.norm( np.asarray(pts[tri[2]])-np.asarray(pts[tri[1]]) ) > L:
            continue
        if np.linalg.norm( np.asarray(pts[tri[0]])-np.asarray(pts[tri[2]]) ) > L:
            continue
        # If we made it this far, the triangle is small enough to keep
        simps_out.append( tri )
    return (pts, simps_out)

def pointCloudToCubes(drift_scan, dt, keepFraction=0.01, size=0.2):
    np.random.seed(1)
    #print 'reading from',scanFile
    #dt=pd.read_csv(scanFile,usecols=[0,1,2],names=['x','y','z'])
    #dt=pd.read_csv(scanFile)
    #print dt['x'][:10]
    #pts=[]
    nPoints=len(dt['x'])
    #pntObj=sg.scaleObj(sg.radiusOneSphereObj,[0.1,0.1,0.1])
    pntObj=scaleObj(Cube,[size,size,size])
    for i in range(nPoints):
        #print "%.1f"%((100.0*i)/nPoints)
        if (i%100==0):
            sys.stdout.write("Scan progress %d%%   \r" % ((100.0*i)/nPoints) )
            sys.stdout.flush()
        if random.random()<keepFraction:
            #pts.append( [dt['x'][i],dt['y'][i],dt['z'][i]] )
            # Note that this repeated merging is not at all efficient, but it's acceptable performance for the number of points we want to keep
            drift_scan=mergeObj(drift_scan,
                                   transObj( pntObj, [dt['x'][i],dt['y'][i],dt['z'][i]] )
                                   )

    print 'Scan complete'
    #pts=np.asarray(pts)
    #drift_scan=sg.convexFromPoints(pts) # This turns a list of points into a polyhedron
    return drift_scan
    
def grepPointCloudToCubes(drift_scan, filename, grepString, keepFraction=0.01, size=0.2):
    np.random.seed(1)
    dt={'x':[],'y':[],'z':[]}
    fd=open(filename,'r')
    for line in fd:
        m = re.match(grepString, line)
        if m:
            x = float(m.group(1))
            y = float(m.group(2))
            z = float(m.group(3))
            dt['x'].append(x)
            dt['y'].append(y)
            dt['z'].append(z)
    fd.close()
    nPoints=len(dt['x'])
    #pntObj=sg.scaleObj(sg.radiusOneSphereObj,[0.1,0.1,0.1])
    pntObj=scaleObj(Cube,[size,size,size])
    for i in range(nPoints):
        #print "%.1f"%((100.0*i)/nPoints)
        if (i%100==0):
            sys.stdout.write("Scan progress %d%%   \r" % ((100.0*i)/nPoints) )
            sys.stdout.flush()
        if random.random()<keepFraction:
            #pts.append( [dt['x'][i],dt['y'][i],dt['z'][i]] )
            # Note that this repeated merging is not at all efficient, but it's acceptable performance for the number of points we want to keep
            drift_scan=mergeObj(drift_scan,
                                   transObj( pntObj, [dt['x'][i],dt['y'][i],dt['z'][i]] )
                                   )
    print 'Scan complete'
    #pts=np.asarray(pts)
    #drift_scan=sg.convexFromPoints(pts) # This turns a list of points into a polyhedron
    return drift_scan
    
# From math.stackexchange.com find-shortest-distance-between-lines-in-3d
def shortestDistanceBetweenLines(a,b, c,d):
    # a=origin of first line
    # b=tangent to first line
    # c=origin of second line
    # d=tangent to second line
    print "a",a
    print "b",b
    print "c",c
    print "d",d

    # t=path length along first line
    # s=path length along second line

    e=a-c

    A = -np.dot(b,b)*np.dot(d,d) + np.dot(b,d)*np.dot(b,d)

    # A=0 if the lines are parallel
    
    s = ( -np.dot(b,b)*np.dot(d,e) + np.dot(b,e)*np.dot(d,b) )/A
    t = (  np.dot(d,d)*np.dot(b,e) - np.dot(b,e)*np.dot(d,b) )/A

    dvect=e+b*t-d*s

    dist=np.sqrt( np.dot( dvect, dvect ) )

    return dist

# Place a radial hydraulic fracture of radius r at x0
def HF(r,x0, strikeRad, dipRad, h=0.5):
    # start with a disk
    disk=diskObj(r,h)
    disk=rotateObj(disk,[0.0,1.0,0.0],dipRad)
    disk=rotateObj(disk,[0.0,0.0,1.0],-strikeRad)
    disk=transObj(disk,x0)
    return disk

# Look for intersection of ellipses with a line segment
# A line segment is described by:
#   - Origin l0
#   - Tangent tVec
#   - Length l
# An ellipse is described by:
#   - Center e0
#   - Major axis aVec
#   - Minor axis bVec
# Note: These choices of defininition make the calculations more efficient
# First we locate the intersection of the line with the plane of the ellipse
# Then we check if the point of intersection is BOTH
#   - Inside the ellipse
#   - Inside the line segment
def ellipseIntersectSegmentVect(l0, tVec, l, e0, aVec, bVec):
    # Obtain normal to the ellipse (note that we must non-dimensionalize)
    aMag=np.linalg.norm(aVec)
    aVecN=aVec/aMag
    bMag=np.linalg.norm(bVec)
    bVecN=bVec/bMag
    nVec=np.cross(aVecN,bVecN)
    # Location of intersection along the length of the line segment:
    s,xInt=distOfIntersectionOfLineAndPlane(lineX=l0,lineS=tVec,planeX=e0,planeN=nVec)
    # These will later be set to the distances along the a- and b-axes
    aa=np.Inf; bb=np.Inf
    if s<0.0 or s>l:
        # The point of intersection lies outside the line segment
        return False,s,aa,bb,xInt
    # The intersection is within the line segment
    # Is the intersection inside the ellipse?
    # Need to obtain the projection of the point of intersection onto the aVec and bVec directions
    aa=np.dot((xInt-e0),aVecN)
    bb=np.dot((xInt-e0),bVecN)
    if aa*aa/(aMag*aMag)+bb*bb/(bMag*bMag) > 1.0:
        # We are outside the ellipse
        return False,s,aa,bb,xInt
    # If we made it this far we are inside the ellipse and the line segment
    return True,s,aa,bb,xInt

def ellipseIntersectSegmentEndPts(x0, x1, e0, aVec, bVec):
    l0=x0
    tVec=x1-x0
    l=np.linalg.norm(tVec)
    tVec=tVec/l
    return ellipseIntersectSegmentVect(l0, tVec, l, e0, aVec, bVec)

def diskIntersectSegmentEndPts(x0, x1, c0, strikeDeg, dipDeg, r):
    # Convert our disk into an equivalent ellipse
    aVec, bVec=strikeDipToAxes(strikeDeg,dipDeg,r)
    intersects,s,aa,bb,xInt=ellipseIntersectSegmentEndPts(x0, x1, c0, aVec, bVec)
    return intersects,xInt

# Convert strike and dip into major and minor axes of an ellipse
def strikeDipToAxes(strikeDeg,dipDeg,radius):
    strikeRad=degToRad(strikeDeg)
    dipRad=degToRad(dipDeg)
    aVec=rotatePoints( [0.0,radius,0.0], [0.0,0.0,1.0], -strikeRad )
    bVec=rotatePoints( [0.0,0.0,radius], aVec, dipRad+0.5*np.pi )
    return aVec, bVec
