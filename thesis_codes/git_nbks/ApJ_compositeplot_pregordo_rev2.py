import pyfits
import numpy as n
#import matplotlib as plt
import matplotlib.pyplot as py
import sys
import os
import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d
import re

'''
July 17th: edited for apj letter, both on same plot
Haven't started edits yet, wait to hear comments on plots to include
Jan 24th, 2016: To reduce white space I make plots in three's, one for each region. I think I will
have to manually enter region titles to insert into bodies of plots N, Sf, S.
May have to make a new version of the program with three subplots, once decide if showing three
Figure out first how to customize and put names inside plots
Dec 7th, 2015:
Shortened and compressed spectra to remove Hbeta and OII for aesthetic purposes, edited in paper
Version with smoothed spectra for aesthics for paper. I think I'll add ivar too, or normalize dynamic range of plots
around median value.
Try scipy.ndimage.filters.gaussian_filter1d() or  scipy.ndimage.filters.gaussian_filter(lambda, sigma) where sigma is in px
px scale is 0.32 A/px restframe, so 0.96Ang smoothed by 3px, I think
Also need to move region cid text to annotate in bottom corner of plot, as Brian's suggested


Aug 28 2015: Edited plot aesthetics with spacing and title, haven't run yet

Editing so number of cats in each composite goes in the titles
'''
#run in folder with .fits ()
#then the targetdir will work
#PLOT .fits file from co-add program coadd_xray_spectra(cat, outspec, normalize=1) using pyfits with key spectral lines
#APJ PAPER VERSION. changes: set max BC plots wrt BC max 3.5, RS norm wrt 1.8, cut off at 5300, title, label height, labels only at top

outfolder='/Users/alison/DLS/deimosspectracatalog/spectrastack/briansprog/preElgordo/Pre_elGordo_paper/cats/cats_groupby/'

#targetdir='/Users/alison/DLS/deimosspectracatalog/Apj_MB/spectranotes/EWcats'
#so don't have to read in from each folder separately

#PROBLEM: this is run in .fits directory, not dir with cats and lens
#could just run this once while plotting from inside folder
#targetdir='/Users/alison/DLS/deimosspectracatalog/ApJ_MB/cats/comp_cats'
#    fitsfile0 ='/Users/alison/DLS/deimosspectracatalog/ApJ_MB/coadds/fits/'+cat0

#targetdir='.'

#folder='1Mpc'
#folder='p5Mpc_shifted'
#folder='Rvir_radialg'
#folder='Rvir_radialg2'
#folder='Rvir_radialg3'
folder='Rvir_radialg4'

#For cut without assigned, reduced range removes this
#labels=['merger back','infall B','core B','merger front','infall A','core A'] #conpare 1+13 to see how diff they are
#bins3=['1','2','3','4','5','6'] #This is to assign new bin numbers for coadd program if don't want to use string names


#labels=['merger back','infall B','core B','merger front','infall A','core A','unassigned'] #conpare 1+13 to see how diff they are
#bins3=['1','2','3','4','5','6','7'] #This is to assign new bin numbers for coadd program if don't want to use string names

labels=['Merger Back','Merger Front'] #conpare 1+13 to see how diff they are
bins3=['1','4'] #This is to assign new bin numbers for coadd program if don't want to use string names
cats0=['coadd_cat_gqz_IDL_1.fits','coadd_cat_gqz_IDL_4.fits'] #override to just plot merger back and front

dic=dict(zip(bins3,labels))

#labels=['Merger Back','Merger Front']
#nums=['1','4']

#cats0=filelist[1:] #These will be fits files. I think 1 starts on right bc axis flipped and monotonic increase
#nums=pd.Series(cats0).apply(lambda x: x.split('_')[4].split('.')[0]).values #extract STRING array of numbers for labeling of bin...in plot of composites by bin

'''
{'1': 'merger back',
 '2': 'infall B',
 '3': 'core B',
 '4': 'merger front',
 '5': 'infall A',
 '6': 'core A'}
'''

targetdir='/Users/alison/DLS/deimosspectracatalog/spectrastack/briansprog/preElgordo/Pre_elGordo_paper/cats/cats_groupby/'+folder
'''
filelist = os.listdir(targetdir+'/coadd_fits/')
cats0=filelist[1:] #These will be fits files. I think 1 starts on right bc axis flipped and monotonic increase
nums=pd.Series(cats0).apply(lambda x: x.split('_')[4].split('.')[0]).values #extract STRING array of numbers for labeling of bin...in plot of composites by bin
'''
#nums correspond to dic[num] to access corresponding string of merger name
#Wastch out because num is zero indexed but string of first entry is '1', do don't loop
#for num in nums: print dic[num],num

#df_list = [pd.read_table(targetdir+'/'+filename, sep=r"\s*") for filename in filelist[1:]]
#open and print 'incat=[', then loop, then ']', then 'outcat=[' loop ']'
'''
cats0_out=[]
for i in enumerate(filelist):
#    print re.sub('[\.\s]','_',re.sub('[(\],tab]', '',filelist[i[0]] ) )+ str(i[0]+1)+'.tab,'
#    file_new=re.sub('[\.\s]','_',re.sub('[(\],tab]', '',filelist[i[0]] ) )+ str(i[0]+1)+'.tab'
    #print re.sub('[\.\s]','_',re.sub('[(\],]', '',filelist[i[0]] ) )+ str(i[0]+1)+'.tab,'
    print 'This is for making IDL script template to paste for incat .tab, outcat .fits, then maybe in/outplots for composite plot'
    #Could also just use this to rename files
    file_new=re.sub('[\.\s]','_',re.sub('[(\],]', '',filelist[i[0]] ) )+ str(i[0]+1)+'.tab'
    cats0_out.append(file_new)
print 'IDL incat paste template .tab', cats0
print 'IDL outcat paste template .fits, will be read into composite plot and used for output savefile png', cats0_out

'''
#maybe zip with enumerate for plotting label
#maybe use the pdf pages?

#tuple of title and length (if add prefix of population include this)
#catlens=[( len(pd.read_table(filename, sep=r"\s*")), filename)
#        for filename in cat0 ] #works
#note these also include the red population files
#use this for title and folder ext
#these are same for fits and tabs
#if iterate can run this from outside
#BACK THIS UP BEFORE START EDITING filelist[1:-1]
#pop=['ALL_color', 'ALL_region', 'BlueBright', 'BlueBrightMassive', 'BlueMassive']

#One more nest, get inside folder then read ?
#Match each in list of filenames? For now just print out lengths of each and name

#pops=['All_color', 'BlueBright', 'BlueBrightMassive', 'BlueMassive']
#titles=['North','South','Superfield','NorthEast']
#cats0=['Comp_NorthRed.fits','Comp_SouthRed.fits','Comp_RemainderRed.fits', 'Comp_NorthEastRed.fits']
#cats1=['Comp_NorthBlue.fits','Comp_SouthBlue.fits','Comp_RemainderBlue.fits','Comp_NorthEastBlue.fits']
#cat0='coadd_cat_gqz_IDL(96.49,96.546].fits'

#cat0='coadd_cat_gqz_IDL(96.546,96.603].fits'
#cat0='coadd_cat_gqz_IDL(96.603,96.659].fits'
#cat0='coadd_cat_gqz_IDL(96.659,96.716].fits'
#cat0='coadd_cat_gqz_IDL(96.716,96.772].fits'
#cat0='coadd_cat_gqz_IDL(96.772,96.829].fits'

#temp add NorthGreen.fits
#titles=['North','South','Remainder','NorthEast','North Green']
#cats0=['Comp_NorthRed.fits','Comp_SouthRed.fits','Comp_RemainderRed.fits', 'Comp_NorthEastRed.fits','Comp_NorthRed.fits']
#cats1=['Comp_NorthBlue.fits','Comp_SouthBlue.fits','Comp_RemainderBlue.fits','Comp_NorthEastBlue.fits','Comp_NorthGreen.fits']

#USE THIS FOR FULL REGION BOTH COLORS
#May have to make a new version of the program with three subplots, once decide if showing three
#pops=['ALL_region']
#titles=['North vs South','NorthEast v Superfield','North vs Superfield','South vs Superfield']
#cats0=['Comp_North.fits','Comp_NorthEast.fits','Comp_North.fits','Comp_South.fits']
#cats1=['Comp_South.fits','Comp_Remainder.fits','Comp_Remainder.fits','Comp_Remainder.fits']

#titles=['North vs South','NorthEast v Remainder','North vs Remainder','South vs Remainder','North vs North Green']
#cats0=['Comp_North.fits','Comp_NorthEast.fits','Comp_North.fits','Comp_South.fits','Comp_North.fits']
#cats1=['Comp_South.fits','Comp_Remainder.fits','Comp_Remainder.fits','Comp_Remainder.fits','Comp_Northg.fits']

print cats0
#for bina,binb in zip(dic['merger back'],dic['merger front']):
#for num in nums: print dic[num],num
for cat0,cat1,num,num1 in zip(['coadd_cat_gqz_IDL_1.fits'],['coadd_cat_gqz_IDL_4.fits'],'1','4'):
    hdulist = pyfits.open(targetdir+'/coadd_fits/'+cat0)
    tbdata = hdulist[1].data #data portion of fits table
    lamb0 = tbdata.field('lambda')[0]
    spec0 = tbdata.field('spec')[0]
    ivar0 = tbdata.field('ivar')[0] #Could do plus or minus shaded ivar, fill between
    spec0 = gaussian_filter1d(spec0,7) #Smoothing for aesthetics, I think sig is Ang so 0.32A/px,3px 0.96A, or 0.21A/px in galaxy restframe
    totn0 = tbdata.field('totn')[0] #Make a subplot. For some reason the max for 14 spec is 13? Check...
    totn0 = gaussian_filter1d(totn0,7)
    hdulist.close()

    hdulist1 = pyfits.open(targetdir+'/coadd_fits/'+cat1)
    tbdata1 = hdulist1[1].data #data portion of fits table
    lamb1 = tbdata1.field('lambda')[0]
    spec1 = tbdata1.field('spec')[0]
    ivar1 = tbdata1.field('ivar')[0] #Could do plus or minus shaded ivar, fill between
    spec1 = gaussian_filter1d(spec1,7) #Smoothing for aesthetics, I think sig is Ang so 0.32A/px,3px 0.96A, or 0.21A/px in galaxy restframe
    totn1 = tbdata1.field('totn')[0] #Make a subplot. For some reason the max for 14 spec is 13? Check...
    totn1 = gaussian_filter1d(totn1,7)
    hdulist.close()
    
    #fitsfile1 ='/Users/alison/DLS/deimosspectracatalog/ApJ_MB/coadds/fits/'+cat1
    #fitsfile1 ='/Users/alison/DLS/deimosspectracatalog/ApJ_MB/coadds/'+cat1
    '''
    fitsfile1 ='/Users/alison/DLS/deimosspectracatalog/ApJ_MB/coadds/fits/'+pop+'/'+cat1
    #HAVE TO STRIP .fits for .tab
    #cat1=n.str.rstrip(str(cat1), '.fits') + '.tab'
    cat1=cat1.replace('.fits','.tab')
    len1 ='/Users/alison/DLS/deimosspectracatalog/ApJ_MB/cats/comp_cats/'+pop+'/'+cat1
    catlen1=len(pd.read_table(len1, sep=r"\s*"))
    #print title,cat0,catlen0,cat1,catlen1
    numlist.append(catlen0)
    numlist.append(catlen1)
    catlist.append(cat0)
    catlist.append(cat1)
    
    hdulist = pyfits.open(fitsfile1)
    tbdata = hdulist[1].data #data portion of fits table
    lamb1 = tbdata.field('lambda')[0]
    spec1 = tbdata.field('spec')[0]
    ivar1 = tbdata.field('ivar')[0]
    spec1 = gaussian_filter1d(spec1,3) #Smoothing for aesthetics, I think sig is Ang so 0.32A/px,3px 0.96A, or 0.21A/px in galaxy restframe
    #totn1 = tbdata.field('totn')[0]
    hdulist.close()
    '''
    
    fig=py.figure(2)
    ax=fig.add_subplot(111)
    #py.subplot(111)

    #py.title(title+': '+pt+' population'+'\n'+title1, fontsize='12')
    py.ylabel(r'$f_{\lambda}$ (Arbitrary Units)',fontsize=22)
    ax.set_xlabel(r'$\lambda\ (\AA)$',fontsize=22)
    setmax=4400
    setmin=3700

    #the problem is that the RS spectra cut off super high so make this min the same min for x
    #or maybe I should keep this and it doesn't matter where it cuts off for now
    py.xlim(setmin, setmax)
    #ax.set_yticks([0.5,1.0,1.5,2.0,2.5,3.0])
    py.tick_params(axis='y',labelleft='off',left='off',right='off') 

    #Merger Back, Merger Front
    py.plot(lamb0[spec0>0], (spec0[spec0>0]/n.median(spec0[spec0>0]))-0.5, 'k')    
    py.plot(lamb1[spec1>0], (spec1[spec1>0]/n.median(spec1[spec1>0]))+1.0, 'k')    
    ymin,ymax=0.0,4.0
    #py.ylim()=ymin,ymax
    xmin,xmax=3700,4400
    ax.set_ylim((0.0,4.0))
    ax.set_xlim((3700,4400))
    print py.ylim()
    print py.xlim()
    wave=[3727.61, 3869, 3933.667, 3968.472, 4101.74, 4305, 4340.47]
    label=['[OII] [NeIII]', ' ','  ','CaII H&K','H$\delta$       G-Band',' ','H$\gamma$']    
    height=[0.2,0.2,0.2,0.2,0.2,0.2,0.2]
    height=(n.asarray(height)-0.01).tolist()
    for i in n.arange(len(wave)):
        py.vlines(wave[i],
                ymin,
                ymax,
                linestyle='dashed', linewidth=0.6, color='red' )
        py.text(wave[i]+10,
                ymax-height[i],        
                #            ymax-height[i]+0.1,
                label[i], color='k', fontsize='18')
    py.text(setmax-310, ymin+1.3, 'Merger Back', color='k', fontsize='19')    
    py.text(setmax-310, ymin+3.05, 'Merger Front', color='k', fontsize='19')    

    py.text(setmax-620, ymin+1.36, r'<M$_{U}$> -20.80', color='k', fontsize='15')    
    py.text(setmax-620, ymin+3.08, r'<M$_{U}$> -20.49', color='k', fontsize='15')    

#    py.text(setmax-620, ymin+1.15, r'<SFR> 0.95$\pm$0.24', color='k', fontsize='15')    
#    py.text(setmax-620, ymin+2.89, r'<SFR> 0.36$\pm$0.09', color='k', fontsize='15')    

    py.text(setmax-620, ymin+2.89, r'<SFR> 1.13$\pm$0.28', color='k', fontsize='15')    
    py.text(setmax-620, ymin+1.15, r'<SFR> 2.98$\pm$0.75', color='k', fontsize='15')    

#    py.text(setmax-620, ymin+1.15, r'<SFR> 1.06$\pm$0.12', color='k', fontsize='15')    
#    py.text(setmax-620, ymin+2.89, r'<SFR> 0.42$\pm$0.22', color='k', fontsize='15')    

    py.text(setmax-620, ymin+0.94, r'23 Galaxies', color='k', fontsize='15')    
    py.text(setmax-620, ymin+2.70, r'17 Galaxies', color='k', fontsize='15')    

    ymin,ymax=0.0,4.0
    xmin,xmax=3700,4400
    ax.set_ylim((0.0,4.0))
    ax.set_xlim((3700,4400))
    py.xticks([3800,3900,4000,4100,4200,4300], fontsize = 18)
    print py.ylim()
    print py.xlim()
    py.gcf().subplots_adjust(left=0.08)
    py.gcf().subplots_adjust(bottom=0.13)
    #py.annotate(title, xy=(2, 1), xytext=(3, 1.5), arrowprops=dict(facecolor='black', shrink=0.05), )
    cat0=cat0.replace('.fits','')
    py.savefig(targetdir+'/plots/Dual_MFMB_v2'+cat0)
    #py.close()
    #print 'pop list = ',numlist
    #print 'cat list = ',catlist
    
    

    #hdulist.close()

