# Copyright 2018-2021 Lawrence Livermore National Security, LLC and other
# Fat Crayon Toolkit Project Developers. See the top-level COPYRIGHT file for details.
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
import SimpleGeometry as sg

# Calulate the normal and shear stresses on a plane
def Sn_tau(
        nrmP, # Normal
        sigP # Stress tensor
        ):
    # Let's now get the traction on the fault segment
    # The traction is the force exerted on a plane:
    # t_i = sig_ji n_j
    t=np.zeros([3])
    for i in range(3):
        for j in range(3):
            t[i] += sigP[j][i]*nrmP[j]

    # We can now decompose the traction into the normal and shear components

    # The normal component of the traction is then:
    Sn=0.0;
    for i in range(3):
        Sn += t[i]*nrmP[i]

    # The shear component is the traction less the normal
    tauV=np.zeros([3])
    for i in range(3):
        tauV[i]=t[i]-Sn*nrmP[i]
    tau=np.sqrt(tauV[0]*tauV[0]+tauV[1]*tauV[1]+tauV[2]*tauV[2])
    return Sn,tau

# Calculate the mu you would need for this fracture orientation to be stable
# Also calculate the maximum change in pore pressure before we slip
# Note that this was translated from an old Perl script:
# geoscience_style_stereogram_given_stresses
# Ideally we would make the syntax more Pythonic
def criticalDelPp_stableMu(
        nrmP, # Normal
        sigP, # Stress tensor
        Pp, # Pore pressure
        muSlip # Value of friction coefficient where slip would typically occur (often taken to be 0.6)
        ):

    # Let's now get the traction on the fault segment
    # The traction is the force exerted on a plane:
    # t_i = sig_ji n_j
    t=np.zeros([3])
    for i in range(3):
        for j in range(3):
            t[i] += sigP[j][i]*nrmP[j]

    # We can now decompose the traction into the normal and shear components

    # The normal component of the traction is then:
    Sn=0.0;
    for i in range(3):
        Sn += t[i]*nrmP[i]

    # The shear component is the traction less the normal
    tauV=np.zeros([3])
    for i in range(3):
        tauV[i]=t[i]-Sn*nrmP[i]
    tau=np.sqrt(tauV[0]*tauV[0]+tauV[1]*tauV[1]+tauV[2]*tauV[2])

    # Critical pressure to reactivate fault slip
    Pc=Sn-tau/muSlip
    # Change in pore pressure required for activation
    Pcp=Pc-Pp
    #if (Pcp<0.0):
    #    # We are inherently unstable
    #    Pcp=float("nan")

    # Coefficient of friction required for stability
    # Note: This is also known as slip tendency
    # Either way, compare with a critical slip value (such as 0.6) to decide how close to slip we are
    #if (Pp>Sn):
    #    # We have failed in tension!
    #    mu=float("nan")
    #    #mu=100.0
    #else:
    mu=tau/(Sn-Pp)
    #if (mu>muSlip):
    #    mu=float("nan")

    return Pcp, mu

def nrmFromDipDirDipAngDeg(dip_direction,dip_angle):
        # The projection of the normal to the plane is determined by the dip_direction
        # dip_direction=0 is north -> nrmxH=0, nrmyH=1
        # dip_direction=90 is east -> nrmxH=1, nrmyH=0
        nrmxH=np.sin(np.pi*dip_direction/180.0)
        nrmyH=np.cos(np.pi*dip_direction/180.0)
        # The lateral components of the normals are corrected for the dip angle
        nrmx=nrmxH*np.sin(np.pi*dip_angle/180.0)
        nrmy=nrmyH*np.sin(np.pi*dip_angle/180.0)
        # The vertical
        nrmz=np.cos(np.pi*dip_angle/180.0)
        return np.asarray([nrmx,nrmy,nrmz])
    

def sigGfromPrincipal(Sh,SH,SV,ShAzimuthDeg,ShDipDeg):
    # Stresses in the principal stress directions
    # We have been given the azimuth of the minimum stress, ShAzimuthDeg, to compare with 90 (x-dir)
    # That is Sh==Sxx, SH=Syy, SV=Szz in the principal stress coord system
    sigP=np.asarray([
        [Sh,  0.0, 0.0],
        [0.0, SH,  0.0],
        [0.0, 0.0, SV ]
        ])
    # Want to test rotating stress into global coords
    if True:
        # We have been given the azimuth of the minimum stress, ShAzimuthDeg, to compare with 90 (x-dir)
        # and we rotate about z-axis
        deltaShAz=(ShAzimuthDeg-90.0)*np.pi/180.0
        sigG=deepcopy(sigP)
        # We have been given the dipping angle of the minimum stress, ShDipDeg, to compare with 0 horizontal
        # So we rotate about y-axis
        ShDip=-ShDipDeg*np.pi/180.0
        sigG=sg.rotateTensor(sigG,[0.0,1.0,0.0],-ShDip)
        #print;print "sigP",sigP; print "sigG",sigG;
        sigG=sg.rotateTensor(sigG,[0.0,0.0,1.0],-deltaShAz)
        #print;print "sigP",sigP; print "sigG",sigG; quit()
    return sigG


def plotCritPpStableMu(
        Sh,SH,SV,ShAzimuthDeg,ShDipDeg,Pp, # Stress state
        muSlip, # Rock property (typically 0.6)
        nf_dip_dir_deg_in=[], nf_dip_deg_in=[], # Natrual fracture orientations
        critDelPpFile="Critical_Del_Pp.png", stableMuFile='Stable_Mu_Slip_Tendency.png', mohrsFile="Mohrs_Circle.png",
        plotTitle=None,
        stressMin=0.0, stressMax=15.0e6):
    nf_dip_dir_deg=[]; nf_dip_deg=[]
    for i in range(len(nf_dip_deg_in)):
        #print 'zzz',nf_dip_deg_in[i]
        # Check if the fracture dip angles go past 90 degrees
        if nf_dip_deg_in[i] > 180.0:
            nf_dip_deg.append( nf_dip_deg_in[i]-180.0 )
            nf_dip_dir_deg.append(nf_dip_dir_deg_in[i])
        elif nf_dip_deg_in[i] > 90.0:
            nf_dip_deg.append( 90.0-(nf_dip_deg_in[i]-90.0) )
            nf_dip_dir_deg.append(nf_dip_dir_deg_in[i]+180.0)
        elif nf_dip_deg_in[i] < -90.0:
            nf_dip_deg.append( nf_dip_deg_in[i]+180.0 )
            nf_dip_dir_deg.append(nf_dip_dir_deg_in[i])
        elif nf_dip_deg_in[i] < 0.0:
            nf_dip_deg.append( -nf_dip_deg_in[i] )
            nf_dip_dir_deg.append(nf_dip_dir_deg_in[i]+180.0)
        else:
            nf_dip_deg.append(nf_dip_deg_in[i])
            nf_dip_dir_deg.append(nf_dip_dir_deg_in[i])
    # Get stress state in x,y,z coords from principal stresses
    sigG=sigGfromPrincipal(Sh,SH,SV,ShAzimuthDeg,ShDipDeg)
    # Set up plot arrays
    nRad=100
    dip_angle_deg = np.asarray(range(nRad+1))*(90.0/nRad)
    nTheta=200
    dip_dir_radians = np.asarray(range(nTheta+1))*(2.0*np.pi/nTheta)
    png_dpi=600
    # Published by Wikipedia: https://en.wikipedia.org/wiki/Strike_and_dip
    # One technique is to always take the strike so the dip is 90 deg to the right of the strike, in which case the redundant letter following the dip angle is omitted (right hand rule, or RHR).
    strike_radians=dip_dir_radians-0.5*np.pi
    criticalDelPpG=np.zeros([nRad,nTheta])
    stableMuG=np.zeros([nRad,nTheta])
    for i in range(nRad):
        for j in range(nTheta):
            # Need to convert mathematical theta and r to:
            # * Usable dip and azimuth
            # * Usable normals
            # I've selected plotting options such that theta and azimuth align
            # dip direction is 0 in y-dir
            # dip angle is the radius
            dip_angle=dip_angle_deg[i]
            dip_direction=180.0*dip_dir_radians[j]/np.pi
            # Convert the x and y into a normal to the fracture using dip angle and dip direction
            nrmG=nrmFromDipDirDipAngDeg(dip_direction,dip_angle)
            nrmx,nrmy,nrmz=nrmG[0],nrmG[1],nrmG[2]
            criticalDelPpG[i,j], stableMuG[i,j] = criticalDelPp_stableMu(np.asarray([nrmx,nrmy,nrmz]), sigG, Pp, muSlip)
            #print dip_angle,dip_direction,stableMuG[i,j]

    ax = plt.subplot(111, projection='polar')
    # These will make the angle theta in the plot the same as azimuth
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    #print dip_angle_deg
    plt.pcolormesh(dip_dir_radians,dip_angle_deg,np.ma.masked_where(np.isnan(stableMuG),stableMuG), vmin=0.0, vmax=muSlip,cmap='jet')
    # Plot the natural fractures (if any) we have been given
    plt.plot(np.asarray(nf_dip_dir_deg)*np.pi/180.0,np.asarray(nf_dip_deg[:]),'k*')
    ax.grid(True)
    plt.colorbar()
    if plotTitle==None:
        ax.set_title("Stable Mu (aka Slip tendency) (dip direction, dip angle)", va='bottom')
    else:
        ax.set_title(plotTitle+" (dip direction, dip angle)", va='bottom')
    #plt.show()
    plt.savefig(stableMuFile, format='png',dpi=png_dpi)
    plt.close()

    ax = plt.subplot(111, projection='polar')
    # These will make the angle theta in the plot the same as azimuth
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    # For the lower hemisphere plot we rotate the plot 180 degrees = pi radians
    plt.pcolormesh(dip_dir_radians+np.pi,dip_angle_deg,np.ma.masked_where(np.isnan(stableMuG),stableMuG), vmin=0.0, vmax=muSlip,cmap='jet')
    # Plot the natural fractures (if any) we have been given
    # For the lower hemisphere plot we rotate the plot 180 degrees = pi radians
    plt.plot((np.asarray(nf_dip_dir_deg)+180)*np.pi/180.0,np.asarray(nf_dip_deg[:]),'k*')
    ax.grid(True)
    plt.colorbar()
    if plotTitle==None:
        ax.set_title("Slip tendency - lower hemisphere pole plot", va='bottom')
    else:
        ax.set_title(plotTitle+" (dip direction, dip angle)", va='bottom')
    #plt.show()
    plt.savefig(stableMuFile[:-4]+"_lower_hemisphere_pole_plot.png", format='png',dpi=png_dpi)
    plt.close()


    
    ax = plt.subplot(111, projection='polar')
    # These will make the angle theta in the plot the same as azimuth
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    plt.pcolormesh(dip_dir_radians,dip_angle_deg,np.ma.masked_where(np.isnan(criticalDelPpG),criticalDelPpG), vmin=stressMin, vmax=stressMax,cmap='jet')
    # Plot the natural fractures (if any) we have been given
    plt.plot(np.asarray(nf_dip_dir_deg)*np.pi/180.0,np.asarray(nf_dip_deg[:]),'k*')
    ax.grid(True)
    plt.colorbar()
    if plotTitle==None:
        ax.set_title("Critical Pressure Change (dip direction, dip angle)", va='bottom')
    else:
        ax.set_title(plotTitle+" (dip direction, dip angle)", va='bottom')
    #plt.show()
    plt.savefig(critDelPpFile, format='png',dpi=png_dpi)
    plt.close()

    # And a Mohr's circle would be handy
    # https://en.wikipedia.org/wiki/Mohr%27s_circle#/media/File:Mohr_Circle.svg
    ax = plt.subplot(111)
    # Note that we need to sort the stresses
    # There is no guarantee that SV>SH, for example
    s=[Sh,SH,SV]
    #print 's',s
    s.sort()
    s1=s[2]; s2=s[1]; s3=s[0]
    ax.set_autoscale_on(False)
    ax.set_xlim([0.0,s1])
    ax.set_xlabel('Normal Stress')
    ax.set_ylim([0.0,muSlip*s1])
    ax.set_ylabel('Shear Stress')
    ax.set_aspect(1.0)
    #print 's',s
    ax.add_artist( plt.Circle( ( 0.5*(s1+s3), 0.0), 0.5*(s1-s3), edgecolor='k', facecolor='#ffaaaa' ) )
    ax.add_artist( plt.Circle( ( 0.5*(s2+s3), 0.0), 0.5*(s2-s3), edgecolor='k', facecolor='#aaffaa' ) )
    ax.add_artist( plt.Circle( ( 0.5*(s1+s2), 0.0), 0.5*(s1-s2), edgecolor='k', facecolor='#aaaaff' ) )
    ax.plot([Pp,Pp+s1],[0.0,muSlip*s1])
#    for bore in fracs.keys():
#        print bore
#        frcs=fracs[bore]
#        plt.plot(frcs.dip_dir_deg[:]*np.pi/180.0,frcs.dip_deg[:],'k*')
    # Plot the natural fractures (if any) we have been given
    Sn=np.zeros([len(nf_dip_dir_deg)])
    tau=np.zeros([len(nf_dip_dir_deg)])
    i=0
    for dip_dir_deg,dip_deg in zip( nf_dip_dir_deg, nf_dip_deg ):
        # Get normal to the fracture in global coordinates
        nrmG=nrmFromDipDirDipAngDeg(dip_dir_deg,dip_deg)
        Sn[i],tau[i]=Sn_tau(nrmG, sigG)
        i+=1
    plt.plot(Sn,tau,'k*')
    if plotTitle==None:
        ax.set_title("Mohr's Circle", va='bottom')
    else:
        ax.set_title(plotTitle, va='bottom')
    #plt.show()
    plt.savefig(mohrsFile, format='png',dpi=png_dpi)
    plt.close()
    
