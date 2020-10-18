# -*- coding: utf-8 -*-
# a basic script to read in the dicom file from the Lynx and analyze it- loops through all files that end in dcm in the folder
# in general, all varialble names start with f- file, m-matrix, i-integer, s-string
# sara st. james 
############################
# TO DO list
## have a place to enter the offset X-Y as reported by Lynx.  ## GET THE OFFSETS SORTED OUT :/ IT STINKS RIGHT NOW :/
## DOUBLE CHECK THE X'S AND THE Y'S .. MAKE SURE IT ALL MAKES SENSE, DUDE. 
import pydicom
import scipy.ndimage
from  scipy.optimize import curve_fit
import os
import tkinter, tkinter.filedialog
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import csv
from matplotlib.patches import Ellipse, Rectangle
vIdealLocations=[60,220,380,540]
# define some simple functions that will be used later!
def getKeyX(item):
    return item[0]   
def getKeyY(item):
    return item[1]
def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
# select the files.. all in a single directory 
root = tkinter.Tk()
print('Please select a directory with the Lynx .dcm files.')
root.withdraw()    
dirname = tkinter.filedialog.askdirectory(initialdir="./",title='Please select a directory')
if len(dirname ) > 0:
    print("You chose %s" % dirname) 
# convert dirname to string from unicode
os.chdir(dirname)
files = [f for f in os.listdir('.') if f.endswith('.dcm')] # finds all the dicom files in the folder that you choose.
filesOffset = [f for f in os.listdir('.') if f.endswith('.csv')]
### GET THE OFFSETS  ##########################################################
mOffset=[]
with open(filesOffset[0]) as csvfile:
    readOffset=csv.reader(csvfile, delimiter=',')
    fOffsetcsv=list(readOffset)
    lOffsetX=fOffsetcsv[1]
    lOffsetY=fOffsetcsv[0] 
npOffsetX=np.array(lOffsetX)
npOffsetY=np.array(lOffsetY)

iOffsetX=npOffsetX[1].astype(np.float)
iOffsetY=-1*npOffsetY[1].astype(np.float)

print("iOffsetX is %f" % iOffsetX)
print("iOffsetY is %f" % -iOffsetY)

#### now loop through the dicom files from lynx and analyze each one ########
#### all of the results will be reported in Lynx.csv
mDataXYSxSy=[]
mDataXY=[]
iCounter=0 
iNumFiles=len(files)
while iCounter < iNumFiles:
    sFilename=files[iCounter]
    fLynx=pydicom.read_file(sFilename)
    #set up the output image
    fig=plt.figure(figsize=(10,10))
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax4 = plt.subplot(224)
    iSize = 600 #width and height of image in pixels
    ax3 = plt.subplot(223)
    iThreshold=200 # the threshold for finding the points
    iBuffer=25
    mPix=fLynx.pixel_array
    # start by thresholding the image
    mThreshold=np.copy(mPix) # copy the image from Lynx
    mThreshold[mThreshold<iThreshold]=0  # set everything below threshold to zero
    labeled_image, number_of_objects = scipy.ndimage.label(mThreshold)
    peak_slices = scipy.ndimage.find_objects(labeled_image)
    mMaximum=[]
    mAmplitude=[]
    mSigma=[]
    mCentroidsPixels=[]
    mCentroidTest=[]
    mCentroidSigma=[]
    mCentroids=[]
    mSpots=[]
    mFWHMx=np.zeros((2,16), dtype=float)
    mFWHMy=np.zeros((2,16), dtype=float)
    mSigmax=np.zeros((2,16), dtype=float)
    mSigmay=np.zeros((2,16), dtype=float)
    iIndex=0
    for peak_slice in peak_slices:
        dy,dx  = peak_slice
        x,y = dx.start, dy.start
        #cx,cy = centroid(mThreshold[peak_slice])
        cx,cy=scipy.ndimage.measurements.center_of_mass(mThreshold[peak_slice])
        cx1=round(cx+x, 0)
        cy1=round(cy+y, 0)  
        cxx,cyy=scipy.ndimage.measurements.maximum_position(mThreshold[peak_slice])
        cx2=round(cxx+x, 0)
        cy2=round(cyy+y, 0)
        vYn=np.linspace(cy1-iBuffer,cy1+iBuffer,51)
        vXn=np.linspace(cx1-iBuffer,cx1+iBuffer,51)
        vProfilex=mPix[cx1-iBuffer:cx1+iBuffer+1,cy1]
        vProfiley=mPix[cx1,cy1-iBuffer:cy1+iBuffer+1]
        cy_mm=(cy+y)*0.5
        cx_mm=(cx+x)*0.5
        mCentroids.append((cx_mm,cy_mm))
        mCentroidsPixels.append((cx+x, cy+y))
        # Executing curve_fit on noisy data
        vOptx,vCovx= curve_fit(gauss,vXn,vProfilex,[300,cx1,2])
        vOpty,vCovy= curve_fit(gauss,vYn,vProfiley,[300,cy1,2]) 
        mMaximum.append((cx2,cy2))
        mAmplitude.append((vOptx[0],vOpty[0]))
        mSigma.append((vOptx[2],vOpty[2]))
        mSigmax[1, iIndex]=vOptx[2]*0.5
        mSigmay[1,iIndex]=vOpty[2]*0.5    
        mSigmax[0,iIndex]=iIndex+1
        mSigmay[0,iIndex]=iIndex+1
     
        mSpots.append((cy+y,cx+x,vOptx[2]*0.5,vOpty[2]*0.5)) 
        mCentroidSigma.append((cy+y,cx+x,vOptx[2],vOpty[2]))
        mCentroidsPixels.append((cy+y, cx+x))
        iIndex=iIndex+1        
        
    ###########################################################################
    # now sort the spots so that the order is consistent
    mSpotsSortedX=sorted(mSpots,key=getKeyX)
    # now sort according to cy
    mSpotsSortedY0=np.array(sorted(mSpotsSortedX[0:4],key=getKeyY))
    mSpotsSortedY1=np.array(sorted(mSpotsSortedX[4:8],key=getKeyY))
    mSpotsSortedY2=np.array(sorted(mSpotsSortedX[8:12],key=getKeyY))
    mSpotsSortedY3=np.array(sorted(mSpotsSortedX[12:16],key=getKeyY))
    mSpotsSorted=np.zeros([16,6])
    # determine the deviations from the expected cx  ##########################
    mSpotsSorted[0:4,0:4]=mSpotsSortedY0[:,:]
    mSpotsSorted[4:8,0:4]=mSpotsSortedY1[:,:]
    mSpotsSorted[8:12,0:4]=mSpotsSortedY2[:,:]
    mSpotsSorted[12:16,0:4]=mSpotsSortedY3[:,:]
    ###########################################################################
    mSpotsSorted[0:4,4]=mSpotsSortedY0[:,0]-59
    mSpotsSorted[4:8,4]=mSpotsSortedY1[:,0]-219
    mSpotsSorted[8:12,4]=mSpotsSortedY2[:,0]-379
    mSpotsSorted[12:16,4]=mSpotsSortedY3[:,0]-539
    # determine the deviations from the expected cy ###########################
    mSpotsSorted[0,5]=mSpotsSortedY0[0,1]-59
    mSpotsSorted[1,5]=mSpotsSortedY0[1,1]-219
    mSpotsSorted[2,5]=mSpotsSortedY0[2,1]-379
    mSpotsSorted[3,5]=mSpotsSortedY0[3,1]-539
    mSpotsSorted[4,5]=mSpotsSortedY1[0,1]-59
    mSpotsSorted[5,5]=mSpotsSortedY1[1,1]-219
    mSpotsSorted[6,5]=mSpotsSortedY1[2,1]-379
    mSpotsSorted[7,5]=mSpotsSortedY1[3,1]-539
    mSpotsSorted[8,5]=mSpotsSortedY2[0,1]-59
    mSpotsSorted[9,5]=mSpotsSortedY2[1,1]-219
    mSpotsSorted[10,5]=mSpotsSortedY2[2,1]-379
    mSpotsSorted[11,5]=mSpotsSortedY2[3,1]-539
    mSpotsSorted[12,5]=mSpotsSortedY3[0,1]-59
    mSpotsSorted[13,5]=mSpotsSortedY3[1,1]-219
    mSpotsSorted[14,5]=mSpotsSortedY3[2,1]-379
    mSpotsSorted[15,5]=mSpotsSortedY3[3,1]-539
    mSpotsSorted[:,5]=mSpotsSorted[:,5]*0.5 # convert to mm
    mSpotsSorted[:,4]=mSpotsSorted[:,4]*0.5  # convert to mm
    # define offset as the average shift
    mSpotsSorted[:,0]=mSpotsSorted[:,0] +iOffsetX
    mSpotsSorted[:,1]=mSpotsSorted[:,1] +iOffsetY
    mSpotsSorted[:,4]=mSpotsSorted[:,4] +iOffsetX
    mSpotsSorted[:,5]=mSpotsSorted[:,5] +iOffsetY
    # calc the uniformity of the spots & mean sigma  ##########################
    mSpotsUniformity=np.zeros([16,6])
    mSpotsUniformity=mSpotsSorted[:,2]/mSpotsSorted[:,3]
    iSpotsUniformityMean=np.mean(mSpotsUniformity)
    iSpotsMeanSigmaX=np.mean(mSpotsSorted[:,2])
    iSpotsMeanSigmaY=np.mean(mSpotsSorted[:,3])    
    iSpotsMeanDiffX=np.mean(mSpotsSorted[:,4])
    iSpotsMeanDiffY=np.mean(mSpotsSorted[:,5])
    # print "max diff x = %f" % np.max(mSpotsSorted[:,4])
    # print "min diff x = %f" % np.min(mSpotsSorted[:,4])
    vSpotsDiffXAbs=np.abs(mSpotsSorted[:,4])
    vSpotsDiffYAbs=np.abs(mSpotsSorted[:,5])
    iSpotsMaxDiffX=np.amax(vSpotsDiffXAbs)
    iSpotsMaxDiffY=np.amax(vSpotsDiffYAbs)
    #Now make the plots:
    for ax in (ax1,ax2,ax3,ax4): ax.clear()
    ax1.set_title('Original image')
    cax1=ax1.imshow(mPix,origin='lower')
    plt.colorbar(cax1,ax=ax1)
    iSpotCounter=0 
    for a in vIdealLocations:
        for b in vIdealLocations: 
            ax2.plot(b- iOffsetY,a- iOffsetX,'wo')
            ax2.annotate(str(iSpotCounter), xy=(a+5,b+5), fontsize=10, color='White')
            iSpotCounter=iSpotCounter+1
    for x,y in mCentroidsPixels:
        ax2.plot(y,x,'b.')
        mDataXY.append((x,y))    
    for x,y,sx,sy in mCentroidSigma:
        mDataXYSxSy.append((x,y, sx,sy))
    ax2.set_title('Thresholded image')
    cax2=ax2.imshow(mThreshold,origin='lower', cmap=cm.jet)
    plt.colorbar(cax2,ax=ax2)
    ax2.set_xlim(0,600)
    ax2.set_ylim(0,600)   
    ax3.set_title('Centroid Location Differences. Red= X . Blue= Y')
    ax3.plot(mSpotsSorted[:,4],'rx')
    ax3.plot(mSpotsSorted[:,5],'bx')
    ax3.set_xlabel('Spot Number')
    ax3.set_ylabel('Delta (mm)')
    ax3.axis([-0.5,17,-3,3])
    ax4.set_title('Sigma (mm)')
    ax4.plot(mSpotsSorted[:,2],'r.')
    ax4.plot(mSpotsSorted[:,3],'bx')
    ax4.set_xlabel('Spot Number')
    ax4.set_ylabel('Sigma (mm)')
    ax4.legend( ["Sigma x", "Sigma y"], loc= 'best')
    ax4.axis([0,17, 2, 8])
    ###########################################################################
    sSigmaX=mSpotsSorted[:,2].tolist()
    sSigmaY=mSpotsSorted[:,3].tolist()
    sCentroidX=mSpotsSorted[:,0].tolist()
    sCentroidY=mSpotsSorted[:,1].tolist()
    sCentroidDiffX=mSpotsSorted[:,4].tolist()
    sCentroidDiffY=mSpotsSorted[:,5].tolist()
    #################  prepare the output  ####################################
    iIndexResults=0;
    lTodaysResultsO=[]
    lTodaysResultsO.append(str(files[iCounter]))
    lTodaysResultsO.append('MEAN_SIGMA_X')
    lTodaysResultsO.append(iSpotsMeanSigmaX)
    lTodaysResultsO.append('MEAN_SIGMA_Y')
    lTodaysResultsO.append(iSpotsMeanSigmaY)
    lTodaysResultsO.append('SIGMA_UNIFORMITY')
    lTodaysResultsO.append(iSpotsUniformityMean)
    lTodaysResultsO.append('MEAN_DIFF_X') 
    lTodaysResultsO.append(iSpotsMeanDiffX)
    lTodaysResultsO.append('MAX_DIFF_X')
    lTodaysResultsO.append(iSpotsMaxDiffX)
    lTodaysResultsO.append('MEAN_DIFF_Y') 
    lTodaysResultsO.append(iSpotsMeanDiffY)
    lTodaysResultsO.append('MAX_DIFF_Y')
    lTodaysResultsO.append(iSpotsMaxDiffY)   
    lTodaysResultsO.append('SIGMA_X')
    while iIndexResults<= 15:
        lTodaysResultsO.append(sSigmaX[iIndexResults])
        iIndexResults=iIndexResults+1
    iIndexResults=0
    lTodaysResultsO.append('SIGMA_Y')
    while iIndexResults<= 15:
        lTodaysResultsO.append(sSigmaY[iIndexResults])
        iIndexResults=iIndexResults+1
    iIndexResults=0
    lTodaysResultsO.append('CENTROID_X')
    while iIndexResults<= 15:
        lTodaysResultsO.append(sCentroidX[iIndexResults])
        iIndexResults=iIndexResults+1
    iIndexResults=0
    lTodaysResultsO.append('CENTROID_Y')
    while iIndexResults<= 15:
        lTodaysResultsO.append(sCentroidY[iIndexResults])
        iIndexResults=iIndexResults+1
    iIndexResults=0
    lTodaysResultsO.append('DIFF_X')
    while iIndexResults<= 15:
        lTodaysResultsO.append(sCentroidDiffX[iIndexResults])
        iIndexResults=iIndexResults+1
    iIndexResults=0
    lTodaysResultsO.append('DIFF_Y')
    while iIndexResults<= 15:
        lTodaysResultsO.append(sCentroidDiffY[iIndexResults])
        iIndexResults=iIndexResults+1    
    with open('Lynx.csv', 'ab') as f:                                    
        writer = csv.writer(f, delimiter=',')                                                       
        writer.writerows([lTodaysResultsO])
    plt.savefig(str(files[iCounter])+'.pdf') # this saves a pdf
    iCounter=iCounter+1
print ('all is well in QA land...')

# create the plots of sigma and spot location for all layers analyzed
plt.figure(figsize=(20,20)) 
  
for x,y, sx, sy in mDataXYSxSy:
    currentAxis = plt.gca()
    plt.gca().set_aspect('equal', adjustable='box')
    uniformityFC='aliceblue'
    uniformityEC='blue'
    alphaFC=0.1
    if sx/sy> 1.1:
        uniformityFC='red'
        alphaFC=0.3
        uniformityEC='red'
    if sx/sy< 0.9:
        uniformityFC='red'
        uniformityEC='red'
        alphaFC=0.3
    currentAxis.add_patch(Ellipse((y, x), 2.355*sy, 2.355*sx, angle=0 , facecolor= uniformityFC , edgecolor= uniformityEC, alpha=alphaFC))   
plt.title('Measured Spot Size (FWHM) = Ellipses.   Red Ellipse = Uniformity > 10 %.')
for a in vIdealLocations:
        for b in vIdealLocations:   
            plt.plot(a-iOffsetY,b-iOffsetX,'k.')
plt.xlabel('Pixels (1 pixel = 0.5 mm)')
plt.ylabel('Pixels (1 pixel = 0.5 mm)')
plt.savefig('AllSigmas.pdf')       

plt.figure(figsize=(20,20))
plt.gca().set_aspect('equal', adjustable='box')
for a in vIdealLocations:
        for b in vIdealLocations: 
            currentAxis = plt.gca()
            currentAxis.add_patch(Rectangle((a-4.25-iOffsetY, b-4.25-iOffsetX), 8.5, 8.5, facecolor="white", edgecolor='red'))   
            
for x,y in mDataXY:
    plt.plot(y,x,'b.', markersize=1)
plt.title('Centroid locations. Red Box = 4 mm Box around ideal location')
plt.xlabel('Pixels (1 pixel = 0.5 mm)')
plt.ylabel('Pixels (1 pixel = 0.5 mm)')
plt.savefig('AllCentroids.pdf')






