import numpy as np
from scipy.integrate import newton_cotes
from scipy import stats
#import graph plotting libraries
import matplotlib.pyplot as mplt
import matplotlib.colors as col
from matplotlib.colors import from_levels_and_colors
import seaborn as sns
#import astronomical data manipulations libraries
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#import libraries for sys/IO management
import multiprocessing as mp
import csv as csv

filterWave = dict()
filterTransmission = dict()
model = np.array([])
redshift = -1
lumDist = -1
raDec_extinction = -1
modelNumAgeMap = np.array([0e6,1e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,1e7,1.5e7,2e7,3e7,4e7,5e7,6e7,7e7,8e7,9e7,1e8,2e8,3.02e8,4.01e8,5.03e8,6.04e8,7.06e8,8.02e8,9.11e8,1.007e9,2014e6,3036e6,4029e6,5053e6,6074e6,8061e6,9028e6,11010e6,14000e6]) / 1e6 #rounding units to Myr
metalicityMap = np.array(['0.0004', '0.004','0.008', '0.02'])  #solar metalicity = 0.01

wave = np.array([15235.9,13645.3,12304.9,11029.6,10364.5,9008.9,7985.4,7651.8,6266.6,5810.1,4709.8,4332.7,4007.5,3351.9,2710.5,2359.6])

#wave = np.array([])
blambda = np.array([])
#waveErr = np.array([2874.18,3940.89,3005.20,4995.87,2917.03,1275.75,2069.25,1517.18,1389.62,2253.99,1398.97,879.58,947.87,550.32,435.65,469.48])/2
NUM_OF_AGE_MODELS = len(modelNumAgeMap)
#print NUM_OF_AGE_MODELS
NUM_OF_METAL_MODELS = len(metalicityMap)
filterName = ['f225w','f275w','f336w','f390w','F435W','F475W','F606W','F625W','F775W','F814W','F850LP','f105w','f110w','f125w','f140w','f160w']
#numOfFilters = len(filterName)
NUM_OF_FILTERS = len(filterName)
hydrogenFilters = []
rms_err = np.array([-1])

modelFluxMetal = np.array([-1])
region = np.array([-1])
separation = -1
rebinShape  = np.array([-1])
observedFluxReg = np.array([[-1],[-1]])
poiMap = np.array([[-1],[-1]])
fitsName = " "
snMap = np.array([[]])
photflam = np.array([])
exptime = np.array([])
fitParm = np.array([[]])


a174_16_arr = np.array([-1]) #an array correcting for inter-milky way extinction

# flags to ensure valid execution
initReg = 0
initSN = 0

homePath = " "
imgName = " "
rawName = " "
trial = " "

def init_config(z,err,radec,cosmology = None):
    global redshift
    global lumDist
    global rms_err
    global wave
    global blambda
    global raDec_extinction
    if cosmology is None :
        print "No cosmology model input. Assumes Lambda-CDM cosmology."
        cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
    redshift = z
    lumDist = (cosmo.luminosity_distance(z).to(u.cm)).value
    rms_err = np.flip(err,0)
    wave /= (1+z)
    wave = np.flip(wave,0)
    blambda = ext_ccm(wave)
    raDec_extinction = radec

def set_hydrogenFilteres(arr):
    global hydrogenFilters
    hydrogenFilters = arr

def set_workingDirectory(home):
    global homePath
    homePath = home

def set_trialID(t):
    global trial
    trial = t

def set_imageName(name):
    global imgName
    imgName = name

def set_rawName(name):
    global rawName
    rawName = name

def get_model():
    global modelFluxMetal
    modelFluxMetal = np.array([ importModel(met) for met in metalicityMap ]) # metalicity * age * filter
    return modelFluxMetal

def set_region(arr,aperture):
    global region
    global separation
    global rebinShape
    region = arr
    separation = aperture
    rebinShape = [ (region[3] - region[2])/separation +1, (region[1] - region[0])/separation +1 ] #shape of y, x

def set_fits(name):
    global fitsName
    fitsName = name

def set_photflam(arr):
    global photflam
    arr = np.flip(arr,0)
    photflam = arr

def set_exptime(arr):
    global exptime
    arr = np.flip(arr,0)
    exptime = arr

def save_config(config_file):
    saveData = dict()
    saveData['redshift'] = redshift
    saveData['lumDist'] = lumDist
    saveData['rms_err'] = rms_err
    saveData['wave'] = wave
    saveData['blambda']=blambda
    saveData['hydrogenFilters']=hydrogenFilters
    saveData['homePath'] = homePath
    saveData['raDec_extinction'] = raDec_extinction
    #saveData['trial'] = trial
    #saveData['imgName'] = imgName
    saveData['region'] = region
    saveData['separation'] = separation
    saveData['rebinShape'] = rebinShape
    saveData['fitsName'] = fitsName
    saveData['a174_16_arr'] = a174_16_arr
    saveData['modelFluxMetal'] = modelFluxMetal
    saveData['observedFluxReg'] = observedFluxReg
    saveData['rawName'] = rawName
    np.save(config_file, saveData)
    print "Configuration is saved"

def load_config(config_file=None):
    if config_file is None:
        print "Error: please enter the config file name. Program forced to stop."
        print 5./0
    else:
        conf = np.load(config_file).item()
        global redshift
        global lumDist
        global rms_err
        global wave
        global blambda
        global hydrogenFilters
        global homePath
        global region
        global separation
        global a174_16_arr
        global modelFluxMetal
        global observedFluxReg
        global rebinShape
        global rawName


        #print "Debug: ", conf

        redshift = conf['redshift']
        lumDist = conf['lumDist']
        rms_err = conf['rms_err']
        wave = conf['wave']
        blambda = conf['blambda']
        hydrogenFilters = conf['hydrogenFilters']
        homePath = conf['homePath']
        region = conf['region']
        separation = conf['separation']
        a174_16_arr = conf['a174_16_arr']
        modelFluxMetal = conf['modelFluxMetal']
        observedFluxReg = conf['observedFluxReg']
        rebinShape = conf['rebinShape']
        raDec_extinction = conf['raDec_extinction']
        rawName = conf['rawName']
        print "Configuration successfully loaded."


def save_fit_parm(parms,fname=None):
    if fname is None:
        parmsFileName = 'regFit_' + imgName +'_parms.npy'
        np.save(parmsFileName,parms)
    else:
        np.save(fname,parms)

def load_fit_parm(fName=None):
    if fName is None:
        print "Please enter file name. Program forced to stop."
        print 5/0
        return -1
    else:
        arg = np.load(fName)
        print "Fitting paramters successfully loaded. Size = ", arg.shape
        return arg

def importData(homePath,metalicity=0.008):
    #the argument homePath for specifying where (in the computer) are the data being stored
    def readFile(fileName):
        ifile = open(fileName,'rt')
        reader = csv.reader(ifile)
        # to detect how many columns are there
        firstLine = next(reader)[0].split()
        container = np.array([]).reshape(0,np.shape(firstLine)[0])
        firstLine = np.array([firstLine]).astype(np.float) #cast string type to float
        container = np.append(container,firstLine,axis=0)
        rowNum = 1
        for row in reader: #noted that this start from the second line, the first line is used to determine num of columns
            rowNum += 1
            if row:
                row = np.array([row[0].split()]).astype(np.float)
                container = np.append(container,row,axis=0)
        return container

    # start reading files with different naming styles and formats
    def readFilters():
        for i in filterName[0:4]:
            fileFullName = str(i + '.UVIS1_plain.tab')
            filePath = homePath  + fileFullName
            data = readFile(filePath)
            filterWave[i] = data[:,1]
            filterTransmission[i] = data[:,2]

        for i in filterName[4:11]:
            fileFullName = str('wfc_' + i + '.dat')
            filePath = homePath + fileFullName
            data = readFile(filePath)
            filterWave[i] = data[:,0]
            filterTransmission[i] = data[:,1]

        for i in filterName[11:16]:
            fileFullName = str(i+ '.IR_plain.tab')
            filePath = homePath + fileFullName
            data = readFile(filePath)
            filterWave[i] = data[:,1]
            filterTransmission[i] = data[:,2]

    #read model
    def readModel(metalicity):
        fileFullName = 'Z=' + metalicity + '_Flux_AgeHomogeneous'
        filePath = homePath+fileFullName
        global model
        model = readFile(filePath)

    readFilters()
    readModel(metalicity)

def importModel(metalicity):
    global a174_16
    print 'Importing model data...'
    importData(homePath,metalicity) #the global variable (type = dictionary of arrays) are now filled with the values in the file
    modelWave = model[:,0]* (1+redshift)
    #modelSpectrum = np.power(10, modelWave*.174/(-2.5))
    modelSpectrum = np.power(10, modelWave*raDec_extinction/(-2.5))
    modelFlux = np.array([]).reshape(0,16)
    alambda = ext_ccm(modelWave)
    a174 = np.power(10, alambda*raDec_extinction/-2.5)
    a174_16 , dum = correctionForFilters(lambdaRange=modelWave,sed=a174)
    #calculate the flux for 33 models with different age parameters
    print "Correcting models for filters properties..."
    for m in range(NUM_OF_AGE_MODELS-1): #for python, range32 is exclusive(not loop to 32), but idl is inclusive (include 32)
        #sed: Spectual Energy Distribution
        np_corr, np0_corr = correctionForFilters(lambdaRange=modelWave,sed=model[:,m]*1e20)
        flux_ = np.array([[np_corr[fil] for fil in filterName]])
        modelFlux = np.append(modelFlux, flux_, axis=0)
    return modelFlux


def correctionForFilters(lambdaRange,sed):
    #dont use 'lambda' as a variable name in python, that is a keyword reserved for lambda-calculus computation
    #to reproduce the idl_tabulate procedure, which was not ordinarily defined in python
    def idl_tabulate(x, f, p=5):
        def eval_newton_cotes(x, f):
            if x.shape[0] < 2 :
                return 0
            rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
            weights = newton_cotes(rn)[0]
            return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
        ret = 0
        for idx in xrange(0, x.shape[0], p - 1) :
            ret += eval_newton_cotes(x[idx:idx + p], f[idx:idx + p])
        return ret

    #now, initialize the required values
    filterTransmission_interpolate = dict()

    lambdaRangeUse = np.array([])
    sedUse = np.array([])
    for i in range(len(lambdaRange)):
        if 1000 < lambdaRange[i] < 20000:
            lambdaRangeUse = np.append(lambdaRangeUse,lambdaRange[i])
            sedUse = np.append(sedUse,sed[i])

    #print lambdaRangeUse[-100:]
    c = 3e10
    abzero = 3.631e-20 * c / np.power( (lambdaRangeUse*1e-8) , 2) *1e-8

    # do interpolation
    for key, value in filterTransmission.iteritems():
        filterTransmission_interpolate[key] = np.interp(lambdaRangeUse,filterWave[key],value)

    # integrate the observed with the filters
    np_ab = dict()
    np_0ab = dict()
    for key, value in filterTransmission.iteritems():
        divisor = idl_tabulate(lambdaRangeUse,lambdaRangeUse*filterTransmission_interpolate[key])
        np_ab[key] = idl_tabulate(lambdaRangeUse,sedUse*lambdaRangeUse*filterTransmission_interpolate[key]) / divisor
        np_0ab[key] = idl_tabulate(lambdaRangeUse,abzero*lambdaRangeUse*filterTransmission_interpolate[key]) / divisor
    return np_ab, np_0ab

def calFluxMag(np_ab, np_0ab):
    mag_ab = dict()
    mag_0ab = dict()
    flux = dict()
    for key,value in np_ab.iteritems():
        mag_ab[key] = -2.5*np.log10(np_ab[key]/np_0ab[key])
        mag_0ab[key] = -2.5*np.log10(np_0ab[key])
        flux[key] = np.power( 10, (-.4 * (mag_ab[key] + mag_0ab[key])) )
    return flux, mag_ab

def get_milky_way_ext():
    global a174_16_arr
    a174_16_arr =  np.array( [ a174_16[key] for key in filterName] )

def read_FitsFile(readRegion):
    dataCube = fits.getdata(homePath+'Avg_convolved_3sersic_cube.fits')
    #the 0 dimension contains the flux in 16 filters, 1 dim is Y-coord, 2 dim is X-coord
    dataCube = dataCube[:,readRegion[2]:readRegion[3]+1,readRegion[0]:readRegion[1]+1] #+1 for python is exclusive but IDL is inclusive
    observation = np.array( [ np.sum(dataCube[i,:,:])/16. for i in range(16) ] ) #the first 16 is the apparent area of the star, while the last 16 is for 16 filters
    def stackAndPlotLog():
        stackImg = np.zeros((dataCube.shape[1],dataCube.shape[2]))
        for i in range(dataCube.shape[1]):
            for j in range(dataCube.shape[2]):
                for k in range(dataCube.shape[0]):
                    stackImg[i,j] += dataCube[k,i,j]
        #plot
        from matplotlib.colors import LogNorm
        mplt.imshow(stackImg, cmap='gray' ,norm=LogNorm())
        mplt.colorbar()
        mplt.xlabel('x')
        mplt.ylabel('y')
        mplt.show()

    #stackAndPlotLog()
    return observation

def mcmcSample(obs,model,poi,ignore=[], sampleNum=30000):
    print "Performing sampling..."
    def whitening(data):
        mean = np.mean(obs)
        stddev = np.std(obs)
        return (data - mean) / stddev

    def chiSquare(p,model):
        sumTemp = 0
        for i in range(model.shape[0]):
            hostExt = np.power(10, (p[1] * blambda[i] / -2.5))
            #sumTemp += ( ((p[0]* model[i] * hostExt - obs[i]) /err_[i]) **2 )
            error = np.sqrt(err_[i]**2 + poi[i]**2)
            sumTemp += ( (( np.exp(p[0])* model[i] * hostExt - obs[i]) / error) **2 )
        return sumTemp/model.shape[0]

    #print chiSquare(np.array([.2323,.01]), model[9,:])
    #initialize the markov chain
    def stateTransition(state0):
        stateNew = np.zeros(state0.shape)
        for i in range(state0.shape[0]-1):
            stateNew[i] = state0[i] + np.random.normal(loc=0,scale=variation[i])
        stateNew[-1] = -1
        return stateNew

    #handle the filter to be discarded
    if (len(ignore) == NUM_OF_FILTERS) :
        print "It ignores all the filters"
        print ignore
        return np.ones((1,5))* -1
    else:
        obs = np.array( [ obs[j] for j in range(NUM_OF_FILTERS) if not(np.isin(j,ignore)) ] )
        model = np.array( [ model[:,:,j] for j in range(NUM_OF_FILTERS) if not(np.isin(j,ignore)) ] )
        print model.shape
        err_ = np.array( [rms_err[j] for j in range(NUM_OF_FILTERS) if not(np.isin(j,ignore)) ] )

    ageIdxHigh = NUM_OF_AGE_MODELS -1
    ageIdxLow = 0
    metalIdxHigh = NUM_OF_METAL_MODELS -1
    metalIdxLow = 0
    extIdxHigh = 1

    #variation = np.array([8e-3,6e-3,1,1]) #p0,p1,age,metalicity
    variation = np.array([5e-1,6e-3,2,1]) #p0,p1,age,metalicity , log normal
    #variation = np.array([3,0.001,2,1]) #p0,p1,age,metalicity , fix ext = 0
    #variation = np.array([1,6e-3,2,1]) #p0,p1,age,metalicity , fix ext = 0
    #variation = np.array([10,0,2,1]) #p0,p1,age,metalicity , fix ext = 0
    state = np.array([2,2e-1,25,2,-1]) #three paramters, p0,Av (extinction), ageIdx , metal, chi-sq; good init: age=15, log normal
    #state = np.array([6,0.1,25,2,-1]) #three paramters, p0,Av (extinction), ageIdx , metal, chi-sq; good init: age=15
    #state = np.array([60,0,25,2,-1]) #three paramters, p0,Av (extinction), ageIdx , metal, chi-sq; good init: age=15
    stateEnergy = chiSquare(state[0:2], model[:,int(state[3]),int(state[2])])
    state[-1] = stateEnergy
    #sampleNum = 10000
    sampleInfo = np.ones( (sampleNum,state.shape[0]) )
    sampleInfo *= -1

    #rejection counter
    rej = 0
    #start sampling
    for i in range(sampleNum):
        #propose a new state
        proposal = stateTransition(state)
        outRange = (not (ageIdxLow <= proposal[2] <= ageIdxHigh) ) or (proposal[0] <= 0) or (not (0 <= proposal[1] <= extIdxHigh) ) or (not( metalIdxLow <= proposal[3] <= metalIdxHigh))
        #print "The state; outside while: ", state
        while outRange :
            rej +=1
            #print "Sample rejected, it is the: ", rej , "-th time"
            #print "propose: ", proposal
            #print "state in while: ", state
            proposal = stateTransition(state)
            #print "state after proposal: ", state
            outRange = (not (ageIdxLow <= proposal[2] <= ageIdxHigh) ) or (proposal[0] <= 0) or (not (0 <= proposal[1] <= extIdxHigh) ) or (not( metalIdxLow <= proposal[3] <= metalIdxHigh))

        #proposalEnergy = chiSquare(proposal[0:2], model[int(proposal[3]),int(proposal[2]),:])
        proposalEnergy = chiSquare(proposal[0:2], model[:,int(proposal[3]),int(proposal[2])])
        proposal[-1] = proposalEnergy
        #calculate the acceptance ratio
        diffEnergy = (proposalEnergy - stateEnergy)
        acceptance = np.exp( -diffEnergy )
        if diffEnergy < 0:
            state = proposal
            stateEnergy = proposalEnergy
        else:
            u = np.random.rand()
            if u < acceptance:
                state = proposal
                stateEnergy = proposalEnergy
        sampleInfo[i,:] = state
        #print "Sample completed with state: ", state
    print "last sample: ", sampleInfo[-1,:]
    return sampleInfo

def ext_ccm(lam):
    Rv = 3.1

    xx = 1e4 / lam
    alam = np.zeros(len(lam))

    for i in range(len(xx)):
        if xx[i] >8:
            alam[i] =5.0

        if 3.3<xx[i]<=8:
            xt = xx[i]
            afac = 1.752 - .316 * xt - .0104 / ((xt-4.67)**2 + .341)
            bfac = -3.090 + 1.825 * xt + 1.206 / ((xt-4.62)**2 + .263)
            Fa = -.04473*(xt-5.9)**2 - .009779*(xt-5.9)**3
            Fb = .2130*(xt - 5.9)**2 + 0.1207*(xt-5.9)**3
            if 5.9 <= xx[i] <=8.0 :
                afac = afac + Fa
                bfac = bfac + Fb
            alam[i] = afac + bfac  /Rv

        if 1.1<=xx[i]<=3.3 :
            yy = xx[i] - 1.82
            afac = 1.0 + .17699*yy - .50447*yy**2 - 0.02427*yy**3+0.72085*yy**4 + 0.01979*yy**5 - .77530*yy**6 + .32999*yy**7
            bfac = 1.41338*yy + 2.28305*yy**2 + 1.07233 * yy **3 - 5.38434*yy**4 - 0.62251*yy**5 + 5.30260*yy**6 - 2.09002*yy**7
            alam[i] = afac + bfac / Rv

        if .3<xx[i]<1.1 :
            yy = xx[i]**1.61
            afac = .574 * yy
            bfrac = -.527*yy
            alam[i] = afac + bfac/Rv

        if xx[i] < .3 : #Abitary extrapolate
            yy = xx[i] ** 1.61
            afac = .574*yy
            bfac = -.527*yy
            alam[i] = afac + bfac / Rv

    return alam


def normalize(data):
    total = np.sum(data)
    return data/total


def contourPlot(x,y,ax,color='b'):
    x = x.astype(np.float)
    y = y.astype(np.float)
    nbins = 23
    #preprocessing: Smoothing, using Guassian kernal density estimation for the data
    data = np.vstack((x,y))
    k = stats.kde.gaussian_kde(data,)
    k.set_bandwidth(bw_method=k.factor*1.7)
    xi, yi = np.mgrid[ x.min():x.max():nbins*1j,y.min():y.max():nbins*1j ]
    zi = k(np.vstack([xi.flatten(), yi.flatten()])).reshape(xi.shape)

    #ax.pcolormesh(xi,yi,zi.reshape(xi.shape))

    #find the density
    levels = 1.0 - np.exp(-0.5 * np.arange(1., 1.6, 0.5) ** 2) #constructing the confidence interval of sigams
    rgba_color = col.colorConverter.to_rgba(color)
    contourMap = [list(rgba_color) for i in levels] + [rgba_color]
    for i,inten in enumerate(levels):
        contourMap[i][-1] *=  float(i) / ( len(levels)+1 )

    intensity, xEdge, yEdge = np.histogram2d(x,y,bins=[nbins,nbins])
    intensity = zi
    intensityFlat = intensity.flatten()
    idx = np.argsort(intensityFlat)[::-1] #decending order
    intensityFlat = intensityFlat[idx] #to prevent sorting again. same as np.sort(intensityFlat)
    normalization =  np.cumsum(intensityFlat)
    normalization /= normalization[-1]
    regionCut = np.empty(len(levels))
    for i, inten in enumerate(levels):
        try:
            regionCut[i] = intensityFlat[normalization <= inten][-1]
        except:
            regionCut[i] = intensityFlat[0]
    m = np.diff(regionCut) == 0
    while np.any(m):
        regionCut[np.where(m)[0][0]] *= 1.0 - 1e-4
        m = np.diff(regionCut) == 0
    regionCut.sort()

    # produce the contour lines
    H2 = intensity.min() + np.zeros((intensity.shape[0] + 4, intensity.shape[1] + 4))
    H2[2:-2, 2:-2] = intensity
    H2[2:-2, 1] = intensity[:, 0]
    H2[2:-2, -2] = intensity[:, -1]
    H2[1, 2:-2] = intensity[0]
    H2[-2, 2:-2] = intensity[-1]
    H2[1, 1] = intensity[0, 0]
    H2[1, -2] = intensity[0, -1]
    H2[-2, 1] = intensity[-1, 0]
    H2[-2, -2] = intensity[-1, -1]

    X1, Y1 = 0.5 * (xEdge[1:] + xEdge[:-1]), 0.5 * (yEdge[1:] + yEdge[:-1])
    X2 = np.concatenate([
        X1[0] + np.array([-2, -1]) * np.diff(X1[:2]),
        X1,
        X1[-1] + np.array([1, 2]) * np.diff(X1[-2:]),
    ])
    Y2 = np.concatenate([
        Y1[0] + np.array([-2, -1]) * np.diff(Y1[:2]),
        Y1,
        Y1[-1] + np.array([1, 2]) * np.diff(Y1[-2:]),
    ])


    #plot the contour lines
    #ax = mplt.subplot()
    #print np.concatenate([[0], regionCut, [intensity.max()*(1+1e-4)]])
    ax.contourf(X2,Y2,H2.T, np.concatenate([[0], regionCut, [intensity.max()*(1+1e-4)]]) , colors=contourMap)



def plotResults(obs,poi, data):
    def covariancePlot(data):
        print "Now plotting graph..."
        data = data.astype(np.float)

        #minIdx = np.argmin(data[1:,4])
        #minData = data[minIdx+1,:]

        burnIn = 8000
        numDim = 4

        #plot configuration
        #sns.set_style("dark")
        sns.set_style("ticks")

        #fig, plotArray = mplt.subplots(nrows=numDim,ncols=numDim)
        from matplotlib.gridspec import GridSpec
        fig = mplt.figure()
        gs = GridSpec(numDim,numDim)
        plotArray = []
        for i in range(numDim):
            tempList = []
            for j in range(numDim):
                tempList.append(mplt.subplot(gs[i,j]))
            plotArray.append(tempList)


        axisParameterMap = {
            0 : 'log Normalization' ,
            1 : '$A_v$ ' ,
            2 : 'Age / $Myr$' ,
            3 : 'Metalicity / $Z_{\odot}$' ,
            4 : '$\chi^2 $' ,
        }


        mean = np.array([ np.mean(data[:,i]) for i in range(numDim) ])
        cov = np.cov(data[:,0:numDim].T)
        print "mean: ", mean, "\n cov: " ,cov

        histSwitch = True
        kdeSwitch = True
        axisRange = []
        axisRange.append( np.linspace(0,1,5) )
        axisRange.append( np.linspace(0,.2,5) )
        axisRange.append( np.linspace(2.2,3.4,5) )
        axisRange.append( np.linspace(2.9,3.7,5) )
        axisRoom = np.array([ .01, .93, .1, .066 ])
        #tickPos = [ [0.1,0.3,.5,0.7,0.9,1.1], [0.1,0.25,0.4,0.55,0.7,0.85] , [5,10,15,20,25], [0,1,2,4] ]
        #tickPos = [ [.5,0.6,.7,.8,0.9], [0.1,0.2,0.3,0.4,0.5] , [5,7,9,11,13], [0,1,2,4] ]
        tickPos = [ [0.1,0.3,.5,0.7,0.9,1.1], [0,0.1,0.2,0.3,0.4,0.5] , [10,15,20,23], [0,1,2,4] ]
        tickLabels = [ tickPos[0], tickPos[1], modelNumAgeMap[tickPos[2]], ["0","20%","40%","100%"] ]



        minimumChiIdx = np.argmin(data[:,-1])
        print "minimumChi: ", data[minimumChiIdx,:]
        print "age: " , modelNumAgeMap[ int(data[minimumChiIdx,2]) ], " Myr"
        print "metalicity: ", metalicityMap[ int(data[minimumChiIdx,3]) ], "solar Z"

        for i in range(numDim):
            for j in range(numDim):
                if i != j and j<i:
                    print "Now plotting: ", i,j
                    #graphPlot = sns.kdeplot(data[burnIn:,j], data[burnIn:,i],shade=True, cmap='Blues' , ax=plotArray[i,j])
                    contourPlot(data[burnIn:,j],data[burnIn:,i],ax=plotArray[i][j])
                    '''
                    plotArray[i][j].set_xticks(tickPos[j])
                    plotArray[i][j].set_yticks(tickPos[i])
                    plotArray[i][j].set_xticklabels(tickLabels[j])
                    plotArray[i][j].set_yticklabels(tickLabels[i])
                    '''
                    #plotArray[i,j].scatter(data[minimumChiIdx,j],data[minimumChiIdx,i] , s=75 , marker='*' , c='yellow')
                    #plotArray[i,j].set_ylim(axisRange[i][0] - axisRoom[i], axisRange[i][-1]+axisRoom[i] )
                    #plotArray[i,j].set_xlim(axisRange[j][0] - axisRoom[j], axisRange[j][-1] + axisRoom[j])
                    if i == 0:
                        pass
                        #plotArray[i,j].set_yticks(axisRange[i] )
                        #plotArray[i,j].set_xticks(axisRange[j] )
                    elif j==0:
                        pass
                        #plotArray[i,j].set_yticks(axisRange[i] )
                        #plotArray[i,j].set_xticks(axisRange[j] )
                    else:
                        pass
                        #plotArray[i,j].set_yticks(axisRange[i] )
                        #plotArray[i,j].set_xticks(axisRange[j] )

                if i == j :
                    if i == 0 :
                        pass
                        #plotArray[i,j].set_xticks(axisRange[j] )
                    else:
                        pass
                        #plotArray[i,j].set_xticks(axisRange[j] )
                    plotArray[i][j].set_yticks([])
                    #plotArray[i,j].set_xlim(axisRange[j][0] - axisRoom[j] ,axisRange[j][-1] + axisRoom[j])
                    if i == 1 or i >= 4:
                        #sns.kdeplot(data[burnIn:,i],ax = plotArray[i][j], bw=(np.max(data[:,i])-np.min(data[:,i]) - 2*xGap)/30, kernel='gau' )
                        #sns.kdeplot(data[burnIn:,i],ax = plotArray[i][j], bw=0.33, kernel='gau' )
                        #plotArray[i,j].set_xlim(np.min(data[:,j])+xGap , np.max(data[:,j])-xGap )
                        sns.distplot(data[burnIn:,i], kde=kdeSwitch,  hist=histSwitch,color='b', ax =plotArray[i][j])
                        #sns.distplot(data[burnIn:1000,i], kde=True, hist=False, ax =plotArray[i][j])
                        '''
                        plotArray[i][j].set_xticks(tickPos[j])
                        plotArray[i][j].set_xticklabels(tickLabels[j])
                        '''
                    elif i==0:
                        #sns.kdeplot(data[burnIn:,i],ax = plotArray[i][j], bw="silverman" , kernel='gau')
                        #sns.distplot(data[burnIn:,i], kde=kdeSwitch, hist=histSwitch, fit=stats.norm, fit_kws={'color' : '#3a67ad'}, ax =plotArray[i][j])
                        sns.distplot(data[burnIn:,i], kde=kdeSwitch, hist=histSwitch,color='b', ax =plotArray[i][j])
                        '''
                        plotArray[i][j].set_xticks(tickPos[j])
                        plotArray[i][j].set_xticklabels(tickLabels[j])
                        '''
                    elif i == numDim-1:
                        sns.distplot(data[burnIn:,i], kde=kdeSwitch,  hist=histSwitch,color='b', ax =plotArray[i][j])
                        '''
                        plotArray[i][j].set_xticks(tickPos[j])
                        plotArray[i][j].set_xticklabels(tickLabels[j])
                        '''
                    else:
                        sns.distplot(data[burnIn:,i], kde=kdeSwitch,  hist=histSwitch,color='b', ax =plotArray[i][j])
                        #sns.distplot(data[burnIn:1000,i], kde=True, hist=False,color='b', ax =plotArray[i][j])
                        #sns.kdeplot(data[burnIn:,i],ax = plotArray[i][j], bw= (np.max(data[:,i])-np.min(data[:,i]) - 2*xGap)/3, kernel='gau')
                        '''
                        plotArray[i][j].set_xticks(tickPos[j])
                        plotArray[i][j].set_xticklabels(tickLabels[j])
                        '''

                if j != 0 and i!=numDim-1:
                    plotArray[i][j].set_xticklabels([])
                    plotArray[i][j].set_yticklabels([])
                elif j == 0 and i !=numDim-1:
                    plotArray[i][j].set_xticklabels([])
                elif i == numDim-1 and j !=0:
                    plotArray[i][j].set_yticklabels([])
                if j == 0:
                    plotArray[i][j].set_ylabel(axisParameterMap[i])

                if i == numDim-1:
                    plotArray[i][j].set_xlabel(axisParameterMap[j])
                    if (j == i):
                        pass
                        #plotArray[i,j].set_xlim(mean[i]-cov[i,i], mean[i]+cov[i,i])
                if j == numDim -1 and i!= numDim-1:
                    pass
                    #plotArray[i,j].set_ylim(mean[j]-cov[j,j], mean[j]+cov[j,j])

                if i == j and i == 0:
                    plotArray[i][j].tick_params(labelleft='off')
                if j>i:
                    fig.delaxes(plotArray[i][j])
                #plotArray[-1,-1].set_xticks( np.round( np.linspace( np.min(data[:,j])+xGap, np.max(data[:,j])-xGap , 5 ),1) )
                plotArray[0][0].set_yticklabels([])
                plotArray[0][0].set_ylabel('')


        #Start overplot
        sns.set_style("darkgrid")
        ax = mplt.subplot(gs[0:2,2:])
        ax.set_xlabel('Rest wavelength [Ang]')
        ax.set_xlim((1000,11500))
        obsPlot = ax.scatter(wave,obs, marker='o', s =50, alpha = 0.5 , label="Observed Data")
        error = np.sqrt(rms_err**2 + poi**2)
        print "background errors: ", rms_err
        print "Poi: ", poi
        print "Total: ", error
        #ax.errorbar(wave,obs,yerr=rms_err,barsabove=True)
        ax.errorbar(wave,obs,yerr=error,barsabove=True)
        #ax.set_ylabel("Flux")
        ax.set_yscale('log')
        m = int(data[minimumChiIdx,2])
        bestFitLABEL = "Best fit, model age = " + str(modelNumAgeMap[m]) + "Myr"
        hostExt = np.power(10, (data[minimumChiIdx,1] * blambda[:] / -2.5))
        bestFitPlot, = ax.plot(wave,np.exp(data[minimumChiIdx,0])*modelFluxMetal[ int(data[minimumChiIdx,3]), int(data[minimumChiIdx,2]) ] * hostExt,color='r' ,lw=3.0, label=bestFitLABEL)

        mplt.legend(handles=[obsPlot, bestFitPlot] , loc=1)
        #end overplot
        sns.despine()
        mplt.show()

    def plotBestFit():
        sns.set_style("darkgrid")
        fig = mplt.figure()
        ax = fig.add_subplot(1,1,1)
        ax3 = fig.add_subplot(414,sharex=ax)
        ax2 = ax.twiny()
        #plot the observation first
        obsPlot = ax.scatter(wave,observedReg, marker='o', s =50, alpha = 0.5 , label="Observed Data")
        #ax.errorbar(np.flip(wave,0),obs,yerr=np.flip(rms_err,0),barsabove=True)
        ax.errorbar(wave,observedReg,yerr=rms_err,barsabove=True)

        ax.set_xlabel('Rest wavelength [Ang]')
        ax.set_xlim((1000,11500))
        #ax.set_xscale('log')

        ax2.set_xlim(ax.get_xlim())
        #ax2.set_xticks(np.flip(wave,0))
        ax2.set_xticks(wave)
        ax2.set_xticklabels(filterName)
        ax2.set_xlabel('Filter')

        ax.set_ylabel("Flux")
        ax.set_yscale('log')

        ax3.set_ylabel('Resisduals')
        ax3.locator_params(nbins=5, axis='y')
        #regionUseString = "Region : "  + str(regionUse+1)

        #plot the models
        m = int(data[minimumChiIdx,2])
        bestFitLABEL = "Best fit, model age = " + str(modelNumAgeMap[m]) + "Myr"
        hostExt = np.power(10, (data[minimumChiIdx,1] * blambda[:] / -2.5))
        modelSED = data[minimumChiIdx,0]*modelFluxMetal[ int(data[minimumChiIdx,3]), int(data[minimumChiIdx,2]) ] * hostExt
        bestFitPlot, = ax.plot(wave, modelSED ,color='r' ,lw=3.0, label=bestFitLABEL)
        '''
        fixAgeParm = [1.24976381e+01, 4.84931894e-05 ,3.10000000e+01, 5.22576462e-01]
        hostExt = np.power(10, (fixAgeParm[1] * blambda[:] / -2.5))
        ageFixSED = fixAgeParm[0]*modelFluxMetal[ int(fixAgeParm[3]), int(fixAgeParm[2]) ] * hostExt
        ageOnlyPlot, = ax.plot(wave,ageFixSED, color='g', lw=3.0, label="Fix 3Gyr")
        '''
        mplt.legend(handles=[obsPlot, bestFitPlot] , loc=1)

        resisduals = ax3.plot(wave, (modelSED - observedReg)/rms_err, color='r' )
        #resisdualsAgeFix = ax3.plot(wave, (ageFixSED - observedReg)/rms_err, color='g')

        mplt.show()

    ageOnlyFitParm = np.zeros((12,2))
    #ageOnlyFitParm[5,:] = [0.49999947451, 17] #metallicity fix 40% solar
    ageOnlyFitParm[5,:] = [0.500021547017, 17] #metallicity fix 1 solar
    ageOnlyFitParm[2,:] = [0.501037624731, 20] #metallicity fix 1 solar
    ageOnlyFitParm[1,:] = [0.500004248224, 19] #metallicity fix 1 solar
    def plotResiduals():
        sns.set_style("darkgrid")
        fig = mplt.figure()
        ax = fig.add_subplot(111)
        ax2 = ax.twiny()

        ax.set_xlabel('Rest wavelength [Ang]')
        ax.set_xlim((1000,11500))
        #ax.set_xscale('log')
        ax.set_ylabel('Resisduals ($\sigma$)')

        ax2.set_xlim(ax.get_xlim())
        #ax2.set_xticks(np.flip(wave,0))
        ax2.set_xticks(wave)
        ax2.set_xticklabels(filterName)
        ax2.set_xlabel('Filter')


        #calculate the resisduals
        m = int(data[minimumChiIdx,2])
        resisduals = np.zeros(16)
        hostExt = np.power(10, (data[minimumChiIdx,1] * blambda[:] / -2.5))
        bestFitLABEL = "Metallicity fit, model age = " + str(modelNumAgeMap[m]) + "Myr " + ",Metalicity = " + metalicityMap[ int(data[minimumChiIdx,3]) ]
        resisduals = ( data[minimumChiIdx,0]*modelFluxMetal[ int(data[minimumChiIdx,3]), int(data[minimumChiIdx,2]) ] * hostExt - obs ) / rms_err
        bestFitPlot, = ax.plot(wave,resisduals,color='r' ,lw=3.0, label=bestFitLABEL)

        #resisduals = ( ageOnlyFitParm[5,0] * modelFluxMetal[ 3, int(ageOnlyFitParm[5,1]) ] - obs   ) / rms_err
        #ageFitLabel = "Age only, model age = " + str(modelNumAgeMap[ int(ageOnlyFitParm[5,1]) ]) + "Myr " + ",Metalicity = 40% $Z_{\odot}$"
        resisduals = ( ageOnlyFitParm[regionUse,0] * modelFluxMetal[ 4, int(ageOnlyFitParm[regionUse,1]) ] - obs   ) / rms_err
        ageFitLabel = "Age only, model age = " + str(modelNumAgeMap[ int(ageOnlyFitParm[regionUse,1]) ]) + "Myr " + ",Metalicity = 1 $Z_{\odot}$"
        ageFitPlot, = ax.plot(wave,resisduals,color='b', lw=2.0, label=ageFitLabel, linestyle='--')



        mplt.legend(handles=[bestFitPlot,ageFitPlot] , loc=1)
        mplt.show()

    minimumChiIdx = np.argmin(data[:,-1])
    #plotBestFit()
    #plotResiduals()
    covariancePlot(data)

def mcmcStepCheck(data):
    mplt.plot(data[:,0])
    mplt.plot(data[:,1])
    mplt.plot(data[:,2])
    mplt.plot(data[:,3])
    mplt.show()



def read_fits(showReg=False):
    #dataCube = fits.getdata(homePath+'Avg_convolved_3sersic_cube.fits')
    dataCube = fits.getdata(homePath+fitsName)
    #the 0 dimension contains the flux in 16 filters, 1 dim is Y-coord, 2 dim is X-coord
    global observedFluxReg
    observedFluxReg = np.array([]).reshape(0,NUM_OF_FILTERS)

    #showRegion = dataCube[:,regionUse[2]:regionUse[3],regionUse[0]:regionUse[1]]
    #np.save('showReg.npy', showRegion)

    for i in range(region[0],region[1]+1,separation):
        for j in range(region[2],region[3]+1,separation):
                observedFlux = dataCube[:,j:j+separation,i:i+separation]
                observedFlux = np.array( [ np.sum(observedFlux[k,:,:])/16. for k in range(NUM_OF_FILTERS) ] )
                observedFlux = np.flip(observedFlux,0)
                observedFlux *= photflam
                observedFlux_corr = np.array([observedFlux / a174_16_arr])
                #observedFlux_corr = np.array( [ observedFlux ] )
                observedFluxReg = np.append(observedFluxReg,observedFlux_corr,axis=0)
    print "Shape of observation region (after rebin): " ,observedFluxReg.shape
    if showReg:
        try:
            img_flip = np.flip(observedFluxReg,axis=1)
            img = np.array([ img_flip[:,k].reshape(rebinShape[1],rebinShape[0]).T for k in range(NUM_OF_FILTERS) ])
            hdu = fits.PrimaryHDU(img)
            hdulist = fits.HDUList([hdu])
            img_n = 'rebinImg_' + imgName +'.fits'
            hdulist.writeto(img_n)
            print "Rebin image drawn, name = ", img_n
        except:
            print "Error in draw rebin image"
    #np.save('regAllRebin_clusterMem2.npy', observedFluxReg)
    return observedFluxReg

def read_poi(showPoi=False):
    dataCube = fits.getdata(homePath+"/noise/"+rawName)
    #the 0 dimension contains the flux in 16 filters, 1 dim is Y-coord, 2 dim is X-coord
    global poiMap
    poiMap = np.array([]).reshape(0,NUM_OF_FILTERS)


    #showRegion = dataCube[:,regionUse[2]:regionUse[3],regionUse[0]:regionUse[1]]
    #np.save('showReg.npy', showRegion)

    for i in range(region[0],region[1]+1,separation):
        for j in range(region[2],region[3]+1,separation):
                poi = dataCube[:,j:j+separation,i:i+separation]
                poi = np.array( [ np.sum(poi[k,:,:])/16. for k in range(NUM_OF_FILTERS) ] )
                poi = np.flip(poi,0)
                poi = poi.reshape((1,NUM_OF_FILTERS))
                poi = np.sqrt(poi*exptime)*photflam/exptime
                poiMap = np.append(poiMap,poi,axis=0)
    print "Shape of poi (after rebin): " ,observedFluxReg.shape

    if showPoi:
        try:
            img_flip = np.flip(poiMap,axis=1)
            img = np.array([ img_flip[:,k].reshape(rebinShape[1],rebinShape[0]).T for k in range(NUM_OF_FILTERS) ])
            hdu = fits.PrimaryHDU(img)
            hdulist = fits.HDUList([hdu])
            img_n = 'rebinPoi_' + imgName +'.fits'
            hdulist.writeto(img_n)
            print "Rebin image drawn, name = ", img_n
        except:
            print "Error in draw rebin image"

    return poiMap


def fitAllReg():
    parms = np.array([]).reshape(0,5)
    nullArr = np.array([-1,-1,-1,-1,-1]).reshape(1,5)
    for i in range(observedFluxReg.shape[0]):
        #avgSig = np.sum( np.abs(observedFluxReg[i,5:] / rms_err[5:]) )/ (16.-5)
        avgSig = np.sum( np.abs(observedFluxReg[i,5:] / rms_err[5:]) ) / (NUM_OF_FILTERS - 5)
        #sig = observedFluxReg[i,:] / rms_err[:]
        #print sig
        if np.all(avgSig > 2):
            #discard shortest filters if signalNoise ratio <= 1 [too noisy]
            #discard = np.array( [ k for k in range(NUM_OF_FILTERS) if (observedFluxReg[i,k] - rms_err[k]) <= 0  ] )
            '''
            if np.any( ( observedFluxReg[i,:5] - rms_err[:5] ) <0 ):
                discard = np.array([0,1,2,3,4])
            else:
                discard = np.array([])
            '''
            samples = mcmcSample(observedFluxReg[i,:], modelFluxMetal,poiMap[i,:], sampleNum=30000)
            minimumChiIdx = np.argmin(samples[:,-1])
            parm = np.array( [ samples[minimumChiIdx,:] ] )
            parms = np.append(parms, parm, axis=0)
        else:
            parms = np.append(parms,nullArr,axis=0)
    return parms

def signalNoise():
    snMap = np.ones((observedFluxReg.shape[0],16))
    snMap *= -1
    for i in range(16):
        snMap[:,i] = observedFluxReg[:,i] / rms_err[i]
    #print snMap.shape
    #print snMap[]
    #np.save("signalNoiseRatio_allReg_ageHomo.npy", snMap)
    return snMap

def drawModelImg(parms):
    print 'Drawing model image'
    shape =  rebinShape
    #print rebinShape
    parms = np.array([ parms[:,i].reshape(shape[1],shape[0]).T for i in range(4) ])

    img = np.ones((16,parms.shape[1],parms.shape[2])) #shape y, shape x
    img *= -1

    for i in range(parms.shape[1]): #y
        for j in range(parms.shape[2]): #x
            if np.all(parms[:,i,j] != -1):
                hostExt = np.power(10, (parms[1,i,j] * blambda[:] / -2.5))
                img[:,i,j] = ( np.exp(parms[0,i,j])*modelFluxMetal[ int(parms[3,i,j]), int(parms[2,i,j]) ] * hostExt )
                img[:,i,j] = np.flip(img[:,i,j],0)
            else:
                img[:,i,j] = np.ones(16)
                img[:,i,j] *= -1
    print 'Storing the img to fits'
    hdu = fits.PrimaryHDU(img)
    hdulist = fits.HDUList([hdu])
    img_n = 'modelImg_' + imgName +'.fits'
    hdulist.writeto(img_n)
    return 0

def checkModel(observedFluxReg,parms, i,j):
    #i: along y axis; j: along x-axis
    '''
    #swap the y-axis, x-axis:
    temp = i
    i = j
    j = temp
    '''
    def plotFittings():
        sns.set_style("darkgrid")
        fig = mplt.figure()
        ax = fig.add_subplot(111)
        ax2 = ax.twiny()
        #plot the observation first
        #obsPlot = ax.scatter(np.flip(wave,0),observedFluxReg, marker='o', s =50, alpha = 0.5 , label="Observed Data")
        obsPlot = ax.scatter(wave,observedFluxReg[:,i,j], marker='o', s =50, alpha = 0.5 , label="Observed Data")
        #ax.errorbar(np.flip(wave,0),observedFluxReg,yerr=np.flip(rms_err,0),barsabove=True)
        ax.errorbar(wave,observedFluxReg[:,i,j],yerr=rms_err,barsabove=True)

        ax.set_xlabel('Rest wavelength [Amg]')
        ax.set_xlim((1000,11500))
        #ax.set_xscale('log')

        ax2.set_xlim(ax.get_xlim())
        #ax2.set_xticks(np.flip(wave,0))
        ax2.set_xticks(wave)
        ax2.set_xticklabels(filterName)
        ax2.set_xlabel('Filter')

        ax.set_ylabel("Intensity /$10^{-19}$")
        ax.set_yscale('log')
        #regionUseString = "Region : "  + str(regionUse+1)

        #plot the models
        m = int(parms[2,i,j])
        bestFitLABEL = "Best fit, model age = " + str(modelNumAgeMap[m]) + "Myr"
        bestFitPlot, = ax.plot(wave,img,color='r' ,lw=3.0, label=bestFitLABEL)

        mplt.legend(handles=[obsPlot, bestFitPlot] , loc=1)
        mplt.show()

    shape =  [ (region[3] - region[2])/separate +1, (region[1] - region[0])/separate +1 ] #shape of y, x
    parms = np.array([ parms[:,k].reshape(shape[1],shape[0]).T for k in range(5) ])
    if np.all(parms[:,i,j] != -1):
        hostExt = np.power(10, (parms[1,i,j] * blambda[:] / -2.5))
        img = ( parms[0,i,j]*modelFluxMetal[ int(parms[3,i,j]), int(parms[2,i,j]) ] * hostExt )
        #img = np.flip(img,0)
        observedFluxReg = np.array([ observedFluxReg[:,k].reshape(shape[1],shape[0]).T for k in range(NUM_OF_FILTERS) ])
        print "Chi-square: (calculated when sampling) " , parms[-1,i,j]
        print "model parameters: " , parms[:-1,i,j]
        print "model prediction: ", img
        print "observation: ", observedFluxReg[:,i,j]
        #print "Chi-square: (calculated when reconstruction) ", (observedFluxReg[:,i,j] - img)/rms_err
        plotFittings()
    else:
        print "Unfitted region: S/N ratio too low"

def drawResidualsImg(parms):
    shape =  rebinShape
    parms = np.array([ parms[:,k].reshape(shape[1],shape[0]).T for k in range(5) ])
    img = np.ones(16)
    img *= -1
    nullArr = img
    nullArr *= 100
    img_ = np.ones((16,parms.shape[1],parms.shape[2]))
    img_ *= -1
    for i in range(parms.shape[1]):
        for j in range(parms.shape[2]):
            if np.all (parms[:,i,j] !=  -1):
                hostExt = np.power(10, (parms[1,i,j] * blambda[:] / -2.5))
                img = ( np.exp(parms[0,i,j])*modelFluxMetal[ int(parms[3,i,j]), int(parms[2,i,j]) ] * hostExt )
                observedFlux = np.array([ observedFluxReg[:,k].reshape(shape[1],shape[0]).T for k in range(NUM_OF_FILTERS) ])
                poiLin = np.array([ poiMap[:,k].reshape(shape[1],shape[0]).T for k in range(NUM_OF_FILTERS) ])
                errorBar = np.sqrt( rms_err**2 + poiLin[:,i,j]**2 )
                resisduals = ( img - observedFlux[:,i,j] ) / errorBar
                resisduals = np.flip(resisduals,0)
                img_[:,i,j] =  resisduals
            else:
                img_[:,i,j] = nullArr
    #print img_
    hdu = fits.PrimaryHDU(img_)
    hdulist = fits.HDUList([hdu])
    img_n = 'resisdualsImg_' + imgName + '.fits'
    hdulist.writeto(img_n)
    return 0

def calibrateModel():

    for metal in metalicityMap[:]:
    #metal = metalicityMap[2]
        rowNum = 0
        fileName = homePath + 'kroupaSpectraList/'+'Z='+metal+"_kroupa_IMF_fcov_0_SFR_inst_Spectra_list_new"
        ifile = open(fileName,'rt')
        reader = csv.reader(ifile)
        container = np.array([]).reshape(0,1)
        for row in reader: #noted that this start from the second line, the first line is used to determine num of columns
            rowNum += 1
            if row:
                row = np.array([row])
                container = np.append(container,row,axis=0)
        #do
        agesModel=np.ones(rowNum)
        agesModel *= -1 #initialize to -1 to distinguish invalid values
        agesModel = np.array([ container[k,0].split('Age',1)[1] for k in range(rowNum) ])
        agesModel.astype(np.float)
        #fluxArr = np.ones(NUM_OF_FILTERS)
        #waveArr = np.array([]).reshape(0,NUM_OF_FILTERS)
        #flux = np.array([]).reshape(0,NUM_OF_FILTERS)
        #there are 47 ages models

        #initialize the fileData, need to make the first column contain the wavelength label
        fileData = np.array([]).reshape(0,2)
        fileName = homePath + 'kroupaSpectraList/'+ container[0,0]
        #print fileName
        ifile = open(fileName,'rt')
        reader = csv.reader(ifile)
        for row in reader: #noted that this start from the second line, the first line is used to determine num of columns
            if row:
                row = np.array([row[0].split()]).astype(np.float)
                fileData = np.append(fileData,row,axis=0)
        fileData[:,1] /= (4*np.pi*lumDist**2) #Luminosity distance in unit of cm
        #can now start loop through all other ages model
        for ageIdx in range(1,rowNum):
            rowNum = 0
            fileName = homePath + 'kroupaSpectraList/'+ container[ageIdx,0]
            #print fileName
            ifile = open(fileName,'rt')
            reader = csv.reader(ifile)
            tempStore = np.array([]).reshape(0,1)
            for row in reader: #noted that this start from the second line, the first line is used to determine num of columns
                rowNum += 1
                if row:
                    row = np.array([row[0].split()]).astype(np.float)
                    row = row[:,1]
                    tempStore = np.vstack((tempStore,row))
                    #tempStore = np.append(tempStore,row,axis=0)
            tempStore[:,0] /= (4*np.pi*lumDist**2)
            fileData = np.append(fileData,tempStore,axis=1)

        print fileData.shape
        #calibration finished, now write back the calibrated file into files
        fNameWrite = homePath + 'Z='+ metal+"_Flux_AgeHomogeneous"
        np.savetxt(fNameWrite, fileData , delimiter=' ')
        print "file \" "  + fNameWrite + "\" saved"

def plotResiduals(data):
    resisduals = data[:,-1].reshape(rebinShape[1],rebinShape[0]).T
    colorMap, colorNorm = from_levels_and_colors([-1,0,0.5,1,2,5,10,50,1000],[(0,0,0,1),(.05,.05,.5,0.3),(.05,.05,.65,0.5),(.05,.05,.8,0.7),(.05,.05,.95,0.9),(.3,.05,.5,1),(.8,.05,.2,1),(1,0,0,1)])
    #mplt.imshow(resisduals ,  interpolation='none', cmap = colorMap, norm=colorNorm)
    mplt.pcolormesh(resisduals ,  cmap = colorMap, norm=colorNorm)
    mplt.colorbar()
    mplt.show()

#mplt.imshow(ages, interpolation='none')
def plotAges(data):

    ages = data[:,2].reshape(rebinShape[1],rebinShape[0]).T
    colorMap, colorNorm = from_levels_and_colors([-1,0,10,15,20,25,30,35,38]
        ,[(0,0,0,1),(.05,.05,.5,0.3),(.05,.05,.65,0.5),(.05,.05,.8,0.7),(.4,.05,.05,0.8),(.7,.05,.05,0.9),(.9,.05,.05,0.9),(1,0,0,1)])
    #this color map is for age
    mplt.pcolormesh(ages, cmap=colorMap, norm=colorNorm)
    #mplt.pcolormesh(ages)
    cbar = mplt.colorbar()
    cbar.ax.set_yticklabels([-1,0,modelNumAgeMap[10],modelNumAgeMap[15],modelNumAgeMap[20], modelNumAgeMap[25],modelNumAgeMap[30],modelNumAgeMap[35],modelNumAgeMap[37]])

    mplt.show()

def plotMetals(data):
    metals = data[:,3].reshape(rebinShape[1],rebinShape[0]).T
    colorMap, colorNorm = from_levels_and_colors([-1,0,2,4],[(0,0,0,1),(.05,.05,.65,0.5),(1,0,0,1)])
    #this color map is for age
    mplt.pcolormesh(metals, cmap=colorMap, norm=colorNorm)
    cbar = mplt.colorbar()
    #cbar.ax.set_yticklabels([-1,0,metalicityMap[2],metalicityMap[3],metalicityMap[4]])
    cbar.ax.set_yticklabels([-1,0,metalicityMap[2],metalicityMap[3]])
    mplt.show()

def plotExtinctions(data):
    extinctions = data[:,1].reshape(rebinShape[1],rebinShape[0]).T
    colorMap, colorNorm = from_levels_and_colors([-1,0,.3,.8,1,10],[(0,0,0,1),(.05,.05,.5,0.3),(.05,.05,.65,0.5),(.05,.05,.8,.5),(1,0,0,1)])
    #this color map is for age
    mplt.pcolormesh(extinctions, cmap=colorMap, norm=colorNorm)
    cbar = mplt.colorbar()
    mplt.show()

def plotMass(data):
    mass = data[:,0].reshape(rebinShape[1],rebinShape[0]).T

    colorMap, colorNorm = from_levels_and_colors([-1,0,.1,.3,.5,1],[(0,0,0,1),(.05,.05,.5,0.3),(.05,.05,.65,0.5),(.05,.05,.8,.5),(1,0,0,1)])
    #this color map is for age
    mplt.pcolormesh(mass, cmap=colorMap, norm=colorNorm)
    cbar = mplt.colorbar()
    mplt.show()

'''
#set output specification
set_trialID(t='clusterMem')
set_imageName(name= 'MACS0329_bcg')

error_bar = np.array([.00708978*2,.00603376*2,.0101939*2,.00861733*2,.00964306*2,.0220976*2,.0153818*2,.0243875*2,.0304843,.0239322*2,.0903711,.157352,.125819,.264335,.467002,.761221])
#init_config(z=0.45, err=error_bar)
init_config(z=0.449, err=error_bar, radec=.1912)
set_workingDirectory(home = '/home/user/Documents/research/galaxyFormation/YGG_DATA/')
corrected_model = get_model()
get_milky_way_ext()

#read observation data
#set_region([685,785,581,669],aperture = 4)
#set_region([423,575, 526,634], aperture=4)
#set_region([1035,1335, 275,455], aperture=4) #cm1
#set_region([1332,1460,1433,1569], aperture=4) #cm3
#set_region([165,269,442,550], aperture=4) #cm8
#set_region([260,400,893,1045], aperture=4) #cm6
#set_region([1132,1284,507,615], aperture=4) #cm7
#set_region([380,580,1537,1754], aperture=4) #spiral
#set_region([1564,1704,487,631], aperture=4) #spiral
set_region([533,609,573,653], aperture=4) #RXJ1437,cm1

#set_fits(name='Avg_convolved_3sersic_CM_cube.fits')
set_fits(name='Avg_convolved_3sersic_CM_cube.fits')
read_fits(False)

#optional: save the configs
#save_config("config_RXJ1437_cm1.npy")
save_config("config_MACS0329_bcg.npy")
'''
#load_config("config_CM6.npy")
set_imageName("MACS0329")
set_rawName("rawDataCubeConvoluted_count.fits")
phflam_temp = np.array([1.91,1.46,2.21,1.52,3.03,15.072126,6.93,9.97,11.9,7.86,17.9,31.8,50.3,134,319,450])
exp_temp = np.array([7741.09,2311.736,2414.669,5026.421,2514.67,10560,14680,4096,3848,3878,10560,4068,4820,4781,7354,7243])
set_photflam(phflam_temp)
set_exptime(exp_temp) #exposure time
read_poi(False)

print poi_map

#do fitting
#parms = fitAllReg()
#save_fit_parm(parms)

#load parms
#parms = load_fit_parm("regFit_clusterMem2_ext0_largeStep2_parms_try.npy")
#parms = load_fit_parm("regFit_bcg_CM-subtract-2_parms_try.npy")
#parms = load_fit_parm("regFit_cm1_parms_try.npy")
#parms = load_fit_parm("regFit_cm3_parms_try.npy")
#parms = load_fit_parm("regFit_CM6_poi_parms.npy")

'''
observedReg = np.array([ observedFluxReg[:,k].reshape(rebinShape[1],rebinShape[0]).T for k in range(NUM_OF_FILTERS) ])
poiMap = np.array([ poiMap[:,k].reshape(rebinShape[1],rebinShape[0]).T for k in range(NUM_OF_FILTERS) ])
observedReg = observedReg[:,20 ,17]
poi = poiMap[:,20,17]
#samples = mcmcSample(observedReg,modelFluxMetal,[0,1,2,3,4],sampleNum=30000)
samples = mcmcSample(observedReg,modelFluxMetal,poi,[],sampleNum=30000)
plotResults(observedReg,poi,samples)
'''
#drawModelImg(parms)
#drawResidualsImg(parms)
#plotAges(parms)
#plotMetals(parms)
#plotResiduals(parms)
#print modelNumAgeMap[35]
