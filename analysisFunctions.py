import numpy as np
import os
import pickle
import pandas as pd
import scipy.interpolate as interpolate
from scipy.stats import bootstrap
import random
import galpy.util.bovy_coords as bc
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from matplotlib.colors import *
from matplotlib.ticker import *
import matplotlib

from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
import astropy.units as un


from gaia_dr3_photometric_uncertainties import EDR3_Photometric_Uncertainties


u = EDR3_Photometric_Uncertainties.Edr3LogMagUncertainty('./gaia_dr3_photometric_uncertainties/LogErrVsMagSpline.csv')

gCalSamples = pd.read_csv('./calSamples/gCalSamples.csv', header=None)
gCalDistSamples = pd.read_csv('./calSamples/gDistCalSamples.csv', header=None)

bCalSamples = pd.read_csv('./calSamples/bCalSamples.csv', header=None)
bCalDistSamples = pd.read_csv('./calSamples/bDistCalSamples.csv', header=None)

rCalSamples = pd.read_csv('./calSamples/rCalSamples.csv', header=None)
rCalDistSamples = pd.read_csv('./calSamples/rDistCalSamples.csv', header=None)

def query(file, rewrite=False):
    return gaiaQuery(file, rewrite)

def calculate(group, rv=-1):
    calcVarG(group)
    calcVarBP(group)
    calcVarRP(group)
    
    return calcVar90(group, rv)

def age(group, distance=True, band='overall', rv=-1):
    varGage, varGageErr, varBPage, varBPageErr, varRPage, varRPageErr = var90Age(group, distance=distance, rv=rv)
    bestAge, bestAgeErr = combineAge(varGage, varGageErr, varBPage, varBPageErr, varRPage, varRPageErr)
    
    if band == 'G': return varGage, varGageErr
        
    if band == 'RP': return varRPage, varRPageErr
        
    if band == 'BP': return varBPage, varBPageErr
    
    if band == 'overall': return bestAge, bestAgeErr
    
    return  bestAge, bestAgeErr, varGage, varGageErr, varBPage, varBPageErr, varRPage, varRPageErr


def analyze(group, band='G'):
    plotXYZ(group, band)
    plotHist(group, band)
    
    return
    

############### HELPER FUNCTIONS BELOW ####################

def filters(group):
    
    # snr cuts 
    bpsnrcut = group['bpsnr'] > 20
    rpsnrcut = group['rpsnr'] > 20
    gsnrcut = group['gsnr'] > 30
    

    # plx cut
    plxcut = group['plx'] / group['plx_err'] > 20
    
    
    # white dwarf cut
    Mg = group['gmags'] - 5*(np.log10(group['dist']) - 1)
    wd = Mg < 10
    wd2 = group['bp-rp'] > 1
    wdcut = wd + wd2 # signals a not WD
    
    # color cut
    colorcut = group['bp-rp'] < 2.5
    
    
    return bpsnrcut * rpsnrcut * gsnrcut  * plxcut  * wdcut * colorcut



def gaiaQuery(file, rewrite=False):
    filename = os.path.split(file)[-1].split(".")[0]
    if os.path.exists(f"./{filename}_gaiaResults.pkl") and not rewrite:
        print("Gaia query results already exist for this file. If you would like to override the previous query, rerun with parameter 'rewrite=True'")
        return pd.read_pickle(f"./{filename}_gaiaResults.pkl")
    
    try:
        FFinfo = pd.read_csv(file)[['Gaia DR3', 'Vr(pred)','Vr(obs)']]
         
    except:
        try: 
            FFinfo = pd.read_csv(file)[['Gaia DR3']]
        except:
            print(f"{file} not found or is in the wrong format")
            return None
    
    bp_rp = []

    gmags = []
    gfluxs = []
    gfluxerrs = []
    gnums = []


    bmags = []
    bfluxs = []
    bfluxerrs = []
    bnums = []


    rmags = []
    rfluxs = []
    rfluxerrs = []
    rnums = []


    ras = []
    decs = []
    pmras = []
    pmdecs = []
    ids = []

    dist = []
    bpsnr = []
    rpsnr = []
    gsnr = []
    ruwe = []
    plx = []
    plxerr = []
    
    rvObs = []
    rvPred = []


    for star in FFinfo.values:
        try:
            r = Gaia.cone_search(f'{star[0]}', radius = 5*un.arcsec, table_name = "gaiadr3.gaia_source")
            r = r.get_results()
            
            # grabs target closest to the search location
            i = np.where(r['dist'] == min(r['dist']))[0][0]
            
            if int(r[i]['phot_bp_n_obs']) != 0 and int(r[i]['phot_rp_n_obs']) != 0:
            
                # put everything in the tables
                bp_rp.append(float(r[i]['bp_rp']))


                gmags.append(float(r[i]['phot_g_mean_mag']))        
                gfluxs.append(float(r[i]['phot_g_mean_flux']))
                gfluxerrs.append(float(r[i]['phot_g_mean_flux_error']))
                gnums.append(int(r[i]['phot_g_n_obs']))


                bmags.append(float(r[i]['phot_bp_mean_mag']))
                bfluxs.append(float(r[i]['phot_bp_mean_flux']))
                bfluxerrs.append(float(r[i]['phot_bp_mean_flux_error']))
                bnums.append(int(r[i]['phot_bp_n_obs']))

                rmags.append(float(r[i]['phot_rp_mean_mag']))
                rfluxs.append(float(r[i]['phot_rp_mean_flux']))
                rfluxerrs.append(float(r[i]['phot_rp_mean_flux_error']))
                rnums.append(int(r[i]['phot_rp_n_obs']))


                ras.append(float(r[i]['ra']))
                decs.append(float(r[i]['dec']))
                pmras.append(float(r[i]['pmra']))
                pmdecs.append(float(r[i]['pmdec']))


                dist.append(1000/float(r[i]['parallax']))
                bpsnr.append(float(r[i]['phot_bp_mean_flux_over_error']))
                rpsnr.append(float(r[i]['phot_rp_mean_flux_over_error']))
                gsnr.append(float(r[i]['phot_g_mean_flux_over_error']))
                ruwe.append(float(r[i]['ruwe']))
                plx.append(float(r[i]['parallax']))
                plxerr.append(float(r[i]['parallax_error']))

                ids.append(f'Gaia DR3 {r[i]["source_id"]}')
                
                try:
                    rvObs.append(star[2])
                    rvPred.append(star[1])
                except:
                    rvObs.append(0)
                    rvPred.append(0)
                    

        except:
            continue
            
            
    headers = ['Gaia DR3', 'ra', 'dec', 'dist', 'plx', 'plx_err', 'pmra', 'pmdec', 'Vr(obs)', 'Vr(pred)', 'ruwe',  'gmags', 'gfluxs', 'gfluxerrs', 'gnums', 'bpmags', 'bpfluxs', 'bpfluxerrs', 'bpnums', 'rpmags', 'rpfluxs', 'rpfluxerrs', 'rpnums', 'bp-rp', 'bpsnr', 'rpsnr', 'gsnr']
    data =    [ids       ,  ras,  decs,   dist,   plx,    plxerr,  pmras,  pmdecs,   rvObs, rvPred, ruwe,    gmags,   gfluxs,   gfluxerrs,   gnums,    bmags,    bfluxs,    bfluxerrs,    bnums,    rmags,    rfluxs,    rfluxerrs,    rnums,   bp_rp,   bpsnr,   rpsnr,   gsnr]

    tempDict = {}
    for i in range(len(headers)):
        tempDict[headers[i]] = data[i]


    df = pd.DataFrame(tempDict)
    
    pickle.dump(df, open(f'./{filename}_gaiaResults.pkl', 'wb'))
    
    print(f"Gaia query results saved to './{filename}_gaiaResults.pkl'")
    
    return df



def calcVarG(group):
    
    nobs = np.arange(np.min(group['gnums'])-1,np.max(group['gnums'])+1,1)
    gn = u.estimate('g',nobs=nobs)
    mag = gn['mag_g']
    
    sig_est = np.zeros(np.size(group['gmags']))


    for i in np.arange(0,np.size(group['gmags'])):
        est = gn[f'logU_{group["gnums"][i]:d}']
        f = interpolate.interp1d(mag,est,fill_value="extrapolate")
        sig_est[i] = (f(group['gmags'][i]))
    
    try:
        group.insert(len(group.columns), "log10(sigma_g_nobs)", sig_est)
    except:
        for v in range(len(group["log10(sigma_g_nobs)"])):
            group.__getitem__("log10(sigma_g_nobs)").__setitem__(v, sig_est[v])
    
    varindx = np.log10(((2.5/np.log(10))*group['gfluxerrs'] / group['gfluxs'])) - (group['log10(sigma_g_nobs)'])
    
    try:
        group.insert(len(group.columns), "varG", varindx)
    except:
        for v in range(len(group["varG"])):
            group.__getitem__("varG").__setitem__(v, varindx[v])
    
    return



def calcVarRP(group):
    
    nobs = np.arange(np.min(group['rpnums'])-1,np.max(group['rpnums'])+1,1)
    gn = u.estimate('rp',nobs=nobs)
    mag = gn['mag_rp']
    
    sig_est = np.zeros(np.size(group['rpmags']))


    for i in np.arange(0,np.size(group['rpmags'])):
        est = gn[f'logU_{group["rpnums"][i]:d}']
        f = interpolate.interp1d(mag,est,fill_value="extrapolate")
        sig_est[i] = (f(group['rpmags'][i]))
    
    try: 
        group.insert(len(group.columns), "log10(sigma_rp_nobs)", sig_est)
    except:
        for v in range(len(group["log10(sigma_rp_nobs)"])):
            group.__getitem__("log10(sigma_rp_nobs)").__setitem__(v, sig_est[v])
    
    varindx = np.log10(((2.5/np.log(10))*group['rpfluxerrs'] / group['rpfluxs'])) - (group['log10(sigma_rp_nobs)'])
    
    try:
        group.insert(len(group.columns), "varRP", varindx)
    except:
        for v in range(len(group["varRP"])):
            group.__getitem__("varRP").__setitem__(v, varindx[v])
    
    return


def calcVarBP(group):
    
    nobs = np.arange(np.min(group['bpnums'])-1,np.max(group['bpnums'])+1,1)
    gn = u.estimate('bp',nobs=nobs)
    mag = gn['mag_bp']
    
    sig_est = np.zeros(np.size(group['bpmags']))


    for i in np.arange(0,np.size(group['bpmags'])):
        est = gn[f'logU_{group["bpnums"][i]:d}']
        f = interpolate.interp1d(mag,est,fill_value="extrapolate")
        sig_est[i] = (f(group['bpmags'][i]))
    
    try:
        group.insert(len(group.columns), "log10(sigma_bp_nobs)", sig_est)
    except:
        for v in range(len(group["log10(sigma_bp_nobs)"])):
            group.__getitem__("log10(sigma_bp_nobs)").__setitem__(v, sig_est[v])
    
    varindx = np.log10(((2.5/np.log(10))*group['bpfluxerrs'] / group['bpfluxs'])) - (group['log10(sigma_bp_nobs)'])
    
    try:
        group.insert(len(group.columns), "varBP", varindx)
    except:
        for v in range(len(group["varBP"])):
            group.__getitem__("varBP").__setitem__(v, varindx[v])
    
    return



def ninety(arr,axis=None):
    return np.nanpercentile(arr,90,axis=axis)


def calcVar90(group, rv = -1):
    
    if rv == False or None:
        rv = -1
    
    if rv == True:
        rv = 5
    
    if rv != -1:
        try:
            rvCut = abs(group['Vr(obs)'] - group['Vr(pred)'] ) < rv
        except:
            rvCut = abs(group['dist']) >= 0 
            print('Error when trying to impliment RV cut. Reverting to no cut.')
    else:
        rvCut = abs(group['dist']) >= 0 # workaround for no rv cuts
        
        
    Gperc90 = np.nanpercentile(group[f"varG"][filters(group) * rvCut ], 90)
    Gper90Err = (bootstrap((group[f"varG"][filters(group) * rvCut ],),ninety)).standard_error

    Bperc90 = np.nanpercentile(group[f"varBP"][filters(group)* rvCut ], 90)
    Bper90Err = (bootstrap((group[f"varBP"][filters(group)* rvCut ],),ninety)).standard_error

    Rperc90 = np.nanpercentile(group[f"varRP"][filters(group)* rvCut ], 90)
    Rper90Err = (bootstrap((group[f"varRP"][filters(group)* rvCut ],),ninety)).standard_error
    
    return Gperc90, Gper90Err, Bperc90, Bper90Err, Rperc90, Rper90Err



def var90Age(group, distance=True, rv=-1):
    
    iterations = 10000
    distVal = 0
    
    Gperc90, Gper90Err, Bperc90, Bper90Err, Rperc90, Rper90Err = calcVar90(group, rv=rv)
    
    if type(distance) != bool:
        if distance == None:
            distance = True
        elif type(distance) == float or type(distance) == int:
            distVal = distance
            distance = True
        
    
    
    if not distance:
        
        # G-band
        trials = []
    
        for i in range(iterations):

            randNum = random.randint(0,len(gCalSamples)-1)

            mSample, bSample, fSample = gCalSamples.iloc[randNum]
            perSample = np.random.normal(Gperc90, Gper90Err)

            logageSample = perSample * mSample + bSample

            trials.append(logageSample)

        trials = np.random.normal(trials,np.random.normal(0.178,0.024,np.size(trials)))


        medA = np.median(trials)
        stdA = np.std(trials)
        medAges = [10**medA, np.asarray([(10**medA - 10**(medA-stdA)), (10**(medA+stdA) - 10**medA)])]

        a = 10**np.percentile(trials, [16, 50, 84])
        b = np.diff(a)
        c = 10**np.percentile(trials, 50)
        
        varGage, varGageErr = c, b
        
        
        
        # Bp-band
        trials = []

        for i in range(iterations):
            randNum = random.randint(0,len(bCalSamples)-1)

            mSample, bSample, fSample = bCalSamples.iloc[randNum]
            perSample = np.random.normal(Bperc90, Bper90Err)

            logageSample = perSample * mSample + bSample

            trials.append(logageSample)

        trials = np.random.normal(trials,np.random.normal(0.177,0.025,np.size(trials)))


        medA = np.median(trials)
        stdA = np.std(trials)
        medAges = [10**medA, np.asarray([(10**medA - 10**(medA-stdA)), (10**(medA+stdA) - 10**medA)])]

        a = 10**np.percentile(trials, [16, 50, 84])
        b = np.diff(a)
        c = 10**np.percentile(trials, 50)

        varBPage, varBPageErr = c, b
        
        
        
        # RP-band
        trials = []

        for i in range(iterations):
            randNum = random.randint(0,len(rCalSamples)-1)

            mSample, bSample, fSample = rCalSamples.iloc[randNum]
            perSample = np.random.normal(Rperc90, Rper90Err)

            logageSample = perSample * mSample + bSample

            trials.append(logageSample)

        trials = np.random.normal(trials,np.random.normal(0.141,0.025,np.size(trials)))

        
        medA = np.median(trials)
        stdA = np.std(trials)
        medAges = [10**medA, np.asarray([(10**medA - 10**(medA-stdA)), (10**(medA+stdA) - 10**medA)])]

        a = 10**np.percentile(trials, [16, 50, 84])
        b = np.diff(a)
        c = 10**np.percentile(trials, 50)

        varRPage, varRPageErr = c, b

    
    
    
    
    if distance:
        if distVal == 0:
            distVal = np.nanmedian(group['dist'])
        
        
        
        # G-Band
        trials = []

        for i in range(iterations):
            randNum = random.randint(0,len(gCalDistSamples)-1)

            mSample, bSample, fSample, dSample = gCalDistSamples.iloc[randNum]
            perSample = np.random.normal(Gperc90, Gper90Err)

            logageSample = perSample * mSample + dSample * distVal + bSample

            trials.append(logageSample)

        trials = np.random.normal(trials,abs(np.random.normal(0.155,0.025,np.size(trials))))

        medA = np.median(trials)
        stdA = np.std(trials)
        medAges = [10**medA, np.asarray([(10**medA - 10**(medA-stdA)), (10**(medA+stdA) - 10**medA)])]

        a = 10**np.percentile(trials, [16, 50, 84])
        b = np.diff(a)
        c = 10**np.percentile(trials, 50)

        varGage, varGageErr = c, b
        
        
        # BP-band
        trials = []

        for i in range(iterations):
            randNum = random.randint(0,len(bCalDistSamples)-1)

            mSample, bSample, fSample, dSample = bCalDistSamples.iloc[randNum]
            perSample = np.random.normal(Bperc90, Bper90Err)

            logageSample = perSample * mSample + dSample * distVal + bSample

            trials.append(logageSample)

        trials = np.random.normal(trials,abs(np.random.normal(0.129,0.024,np.size(trials))))

        medA = np.median(trials)
        stdA = np.std(trials)
        medAges = [10**medA, np.asarray([(10**medA - 10**(medA-stdA)), (10**(medA+stdA) - 10**medA)])]

        a = 10**np.percentile(trials, [16, 50, 84])
        b = np.diff(a)
        c = 10**np.percentile(trials, 50)

        varBPage, varBPageErr = c, b
        
        
        # RP-band
        trials = []

        for i in range(iterations):
            randNum = random.randint(0,len(rCalDistSamples)-1)

            mSample, bSample, fSample, dSample = rCalDistSamples.iloc[randNum]
            perSample = np.random.normal(Rperc90, Rper90Err)

            logageSample = perSample * mSample + dSample * distVal + bSample

            trials.append(logageSample)

        trials = np.random.normal(trials,abs(np.random.normal(0.089,0.023,np.size(trials))))

        medA = np.median(trials)
        stdA = np.std(trials)
        medAges = [10**medA, np.asarray([(10**medA - 10**(medA-stdA)), (10**(medA+stdA) - 10**medA)])]

        a = 10**np.percentile(trials, [16, 50, 84])
        b = np.diff(a)
        c = 10**np.percentile(trials, 50)

        varRPage, varRPageErr = c, b
    
    
    
    return varGage, varGageErr, varBPage, varBPageErr, varRPage, varRPageErr



def combineAge(varGage, varGageErr, varBPage, varBPageErr, varRPage, varRPageErr):
     
    # using lower age error     
    wAvg0 = np.average( np.asarray([varGage, varBPage, varRPage]), weights = np.asarray([1/(varGageErr[0])**2, 1/(varBPageErr[0])**2, 1/(varRPageErr[0])**2]))
    wErr0 = np.sqrt(1/sum(np.asarray([1/(varGageErr[0])**2, 1/(varBPageErr[0])**2, 1/(varRPageErr[0])**2])))

    # using upper age error
    wAvg1 = np.average( np.asarray([varGage, varBPage, varRPage]), weights = np.asarray([1/(varGageErr[1])**2, 1/(varBPageErr[1])**2, 1/(varRPageErr[1])**2]),)
    wErr1 = np.sqrt(1/sum(np.asarray([1/(varGageErr[1])**2, 1/(varBPageErr[1])**2, 1/(varRPageErr[1])**2])))

    # combining
    wAvg = np.average([wAvg0, wAvg1])
    wErr = np.array(list(zip([wErr0], [wErr1]))).T
    
    return wAvg, [wErr[0][0], wErr[1][0]]


def plotXYZ(group, band='G'):
    
    Gllbb = bc.radec_to_lb(group['ra'][filters(group)] , group['dec'][filters(group)] , degree=True)
    Gxyz = bc.lbd_to_XYZ( Gllbb[:,0] , Gllbb[:,1] , group['dist'][filters(group)] , degree=True).T
    
    
    fig,axs = plt.subplots(2,2)
    fig.set_figheight(16)
    fig.set_figwidth(16)
    fig.subplots_adjust(hspace=0.03,wspace=0.03)


    cdict2 = {'red':  [(0.0,  1.0, 0.0),
                       (1.0,  1.0, 0.0),
                       (1.0,  1.0, 0.0)],
             'green': [(0.0,  0.0, 0.0),
                       (1.0, 0.0, 0.0),
                       (1.0,  0.0, 0.0)],
             'blue':  [(0.0,  0.0, 1.0),
                       (1.0,  0.0, 1.0),
                       (1.0,  0.0, 1.0)]} 
    my_cmap2 = matplotlib.colors.LinearSegmentedColormap('my_colormap2',cdict2,256)


    norm = mpl.colors.Normalize(vmin=min(group[f'var{band}'][filters(group)]), vmax=max(group[f'var{band}'][filters(group)]))
    scalarMap = cmx.ScalarMappable(norm=norm,cmap=my_cmap2)

    color = scalarMap.to_rgba(group[f'var{band}'][filters(group)])

    for i in range(len(color)):
        axs[0,0].plot( Gxyz[0][i] , Gxyz[1][i] , 'o' , markersize=18 , mew=3, color = color[i])
        axs[0,1].plot( Gxyz[2][i] , Gxyz[1][i] , 'o' , markersize=18 , mew=3, color = color[i])
        axs[1,0].plot( Gxyz[0][i] , Gxyz[2][i] , 'o' , markersize=18 , mew=3, color = color[i])

    axs[0,0].set_xlabel(r'$X$ (pc)',fontsize=20,labelpad=10)
    axs[0,0].set_ylabel(r'$Y$ (pc)',fontsize=20,labelpad=10)

    axs[1,0].set_xlabel(r'$X$ (pc)',fontsize=20,labelpad=10)
    axs[1,0].set_ylabel(r'$Z$ (pc)',fontsize=20,labelpad=10)

    axs[0,1].set_xlabel(r'$Z$ (pc)',fontsize=20,labelpad=10)
    axs[0,1].set_ylabel(r'$Y$ (pc)',fontsize=20,labelpad=10)

    axs[0,0].xaxis.set_ticks_position('top')
    axs[0,1].xaxis.set_ticks_position('top')
    axs[0,1].yaxis.set_ticks_position('right')

    axs[0,0].xaxis.set_label_position('top')
    axs[0,1].xaxis.set_label_position('top')
    axs[0,1].yaxis.set_label_position('right')

    for aa in [0,1]:
        for bb in [0,1]:
            axs[aa,bb].tick_params(top=True,bottom=True,left=True,right=True,direction='in',labelsize=18)

    fig.delaxes(axs[1][1])
    strsize = 26



    cbaxes = fig.add_axes([0.55,0.14,0.02,0.34])
    cb = plt.colorbar( scalarMap , cax=cbaxes )
    cb.set_label( label=f'var{band}' , fontsize=24 , labelpad=20 )
    cb.ax.tick_params(labelsize=18)

    plt.show()
    
    return


def plotHist(group, band='G'):
    plt.figure(figsize=(10,8))

    counts, bins = np.histogram(group[f'var{band}'][filters(group)])
    plt.stairs(counts, bins, fill=True, alpha = .75)
    
    Gperc90, Gper90Err, Bperc90, Bper90Err, Rperc90, Rper90Err = calcVar90(group)
    
    if band=='G': 
        plt.vlines(Gperc90, 0, max(counts), color='orange')
    elif band=='BP':
        plt.vlines(Bperc90, 0, max(counts), color='orange')
    elif band=='RP':
        plt.vlines(Rperc90, 0, max(counts), color='orange')

    plt.xlabel(f'var{band}', fontsize=18 )
    plt.ylabel('Count', fontsize=18 )

    plt.show()
