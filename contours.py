'''
These are a few functions that may be helpful to an observational cosmologist, but there are many applications. My goal was to take irregular 2D contours and determine which objects from a data set were inside. I found the "point in polygon" method to work well, on even the most irregularly shaped or sparsely defined contours. The functions are specific to my catalog formats, but you can easily edit them, or I will generalize them once I have more time (as a grad student extra time is unlikely). The first four functions work together in sequence. The last two are more primitive and work for rectangular regions. The great part is that the resulting DataFrame allows for easy selection by cross-section (for ex. df.xs(['North'],level=['region_c']) ) for plotting/analysis. Just make the resulting 'region_c' a MultiIndex using df.set_index(drop=False, append=True, inplace=True). I have a ton of other functions that work together in this way and will publish them eventually in an ipython notebook. They are all ideal for the dirty work necessary to prepare data for scientific analysis.

The functions are not as elegant as they could be because I had to quickly rewrite them after I upgraded matplotlib and pandas, which of course broke everything (as the horror story of upgrading usually goes). I am posting them in hopes of helping someone who may have had a similar issue when they upgraded, to save her/him lots of hair-pulling.

Enjoy!
Alison Mansheim
UC Davis Physics (Cosmology) PhD canditate
asmansheim at ucdavis dot edu
'''

import numpy as np
import matplotlib.pyplot as py
import pandas as pd

def make_cont(file1, file2=None):
    '''
    Takes file(s) with arbitrary chunks of coordinates, makes DataFrame grouped by assigned contour number
    Use as input for cont_plot() and ppoly()
    grouped=make_cont(file1,file2=file2)
    Inputs:
    file1 format is ds9 .con: alpha delta coordintes without a header
    file2 is optional second contour file to append onto end
    file2 currenty CSV (easy from DS9 region file, change to .con delete everything but middle comma)
    Outputs:
    Option to pickle
    Outputs grouped contours (pandas grouby object) ready for plotting and/or point in poly functions below
    '''
    import pickle

    conts=pd.read_table(file1, sep='\s',names=['alpha','delta'])
    cutnan=np.where(conts.alpha.isnull() )[0]
    nconts=len(np.split(conts.index, conts.alpha.loc[cutnan].index))
    dfs=[conts.loc[np.split(conts.index, conts.alpha.loc[cutnan].index)[i]] for i in range(nconts)] 
    cont=[c.dropna() for c in dfs if not isinstance(c, np.ndarray)]
    cont.sort(key=lambda c: len(c))
    cont=pd.concat(cont, names=['contour'],keys=np.arange(len(cont)) )
    n=cont.index.max()[0]+1 #need this if 'contour' starts at 1, not 0
    cont=cont.reset_index().drop(['level_1'],axis=1).set_index('contour',drop=False)
    if (file2!=None):
        contR=pd.read_csv(file2, names=['alpha','delta'])       
        #contR=pd.read_csv(file2, names=['alpha','delta'],header=False) #use if header    
        #contR=pd.read_table(file2, names=['alpha','delta']) #use if not csv
        A=np.repeat(n,len(contR)) #adds nth contour
        contR.index=A
        contR.index.names=['contour']
        contR=contR.reset_index().set_index('contour',drop=False)
        cont=pd.merge(cont,contR,how='outer',on=['alpha','delta','contour'],right_index=True,left_index=True)
    grouped=cont.groupby(level=0)
    print 'pickling: cont.pickle'
    filename = 'cont.pickle'
    usepath = '/User/'
    f = open(userpath+filename,'wb')
    pickle.dump(cont,f)
    f.close()
    return grouped

def ppoly(group,X,zgroup=None):
    '''
    REPLACED with Path, contains_points after nxutils depreciated points_inside_poly() in MPL upgrade
    Takes groupby object of contours returned by make_cont()
    Determines contained objects from input cat X for each contour
    Option to make plot of each contour with contained objects
    Read output region_select() function with user defined, desired contours
    Isodensity_contours=grouped.apply(ppoly,X)
    Input:
    X=data[['objid','alpha','delta']]
    X=X.reset_index().set_index('objid',drop=False) #make ID multiindex
    Output:
    DataFrame of contours and df.N with number of objs in each
    '''
    from matplotlib.path import Path
    path = Path(group[['alpha','delta']])
    mask = path.contains_points(X[['alpha','delta']])
    Y=X[mask]
    Y['N']=len(X[mask])
    Y['contour']=group.name
    print 'contour N=',len(group),', #',group.name,'contains:',len(X[mask]),'gals out of',len(X)
    if zgroup!=None:
        fig = py.figure()
        py.plot(X[mask].alpha,X[mask].delta,'ro',ms=2,alpha=0.3)
        py.plot(group.alpha,group.delta,'b-',label=str(group.name))
        py.gca().invert_xaxis()
        py.xlabel('RA')
        py.ylabel('DEC')
        py.tick_params(axis='both', which='major', labelsize=10)
        py.ticklabel_format(useOffset=False)
        title='Contour '+str(group.name)+' contains: '+str(len(X[mask]))+' galaxies'
        py.title(title)
        py.show()
        filename = 'contour_'+str(group.name)+zgroup
        py.savefig(filename)
    return Y

def cont_plot(grouped,annotate=True):
    '''
    Plots groupby objects of contours returned by make_cont()
    cont_plot(grouped,annotate=False)
    cont_plot(grouped)
    Option to annotate numbers contours for further grouping by region in reg_contour()
    '''
    for a,b in grouped:
        py.plot(b.alpha,b.delta,'k-',label=str(a))
        if (annotate==True):
            py.annotate(str(a), xy=(np.array(b.alpha.values)[0], np.array(b.delta.values)[0]), size=20)
    py.gca().invert_xaxis()
    py.xlabel('RA')
    py.ylabel('Dec')
    py.tick_params(axis='both', which='major', labelsize=10)
    py.ticklabel_format(useOffset=False)

    
def reg_contour(data_rs,Iso_weight):
    '''
    Takes DataFrame of contours with contained objects from ppoly()
    cat_all=reg_contour(cat_all,Iso_weight_all)
    Must specify which contours to select in function below
    Output:
    DataFrame with column for region assignment based on contours.
    Below are example contour selections and labels.
    Select from DataFrame easily by making region Multiindex and taking cross-section with .xs :)
    '''
    R=[49]    
    N1c=[48]
    S1c=[47]
    NE=[45,43]
    data_rs['region_c']='Remainder' #initialize with biggest contour if multiple overlaps
    list_all=[R]
    list_reg=['Remainder']
    for reg,name in zip(list_all,list_reg):
        list_dfs=[Iso_weight.iloc[Iso_weight.index.get_level_values('contour')==a] for a in reg]
        reg_df=pd.concat(list_dfs)
        data_rs['region_c'][data_rs.objid.isin(reg_df.objid) ]=name
    list_all=[N1c, S1c, NE]
    list_reg=['North', 'South', 'NorthEast']
    for reg,name in zip(list_all,list_reg):
        list_dfs=[Iso_weight.iloc[Iso_weight.index.get_level_values('contour')==a] for a in reg]
        reg_df=pd.concat(list_dfs)
        data_rs['region_c'][data_rs.objid.isin(reg_df.objid) ]=name
    return data_rs

def func_reg(data_rs, params):
    '''
    This is an older function that does assignments for rectangular regions, unrotated
    Input:
    DataFrame, params file with min and max region bounds, read in below within region_mask()
    Output:
    DataFrame with region assignments
    '''
    reg = strcid( params.cid )
    alpha_min=pd.Series(np.repeat(np.array([params.alpha_min]),len(data_rs)))
    alpha_max=pd.Series(np.repeat(np.array([params.alpha_max]),len(data_rs)))
    delta_min=pd.Series(np.repeat(np.array([params.delta_min]),len(data_rs)))
    delta_max=pd.Series(np.repeat(np.array([params.delta_max]),len(data_rs)))
    condition = ( ( (data_rs.alpha >= alpha_min)&(data_rs.alpha <= alpha_max) ) &  ( (data_rs.delta >= delta_min)&(data_rs.delta <= delta_max) ) )
    data_rs.loc[condition,'region'] = reg
    return data_rs

def region_mask(data_rs):
    file='/Users/params.tab'
    params=pd.read_table(file, skiprows=1)
    data_rs['region']=pd.Series(np.repeat(np.array(['Remainder']),len(data_rs)))
    data_rs = func_reg(data_rs, params[params['N']==0])
    print data_rs.groupby(['region']).size()
    return data_rs
