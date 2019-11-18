# Program to read input data and split it into training and testing datasets. -> New approach
import ROOT
import array
import numpy as np

ENR=[20,50,80,100,120,200,250,300]
# ENR=[20]
INPATH1 ="IN_ROOT/Pion_14_11_19_Sim_"
INPATH2 ="_ForML.root"
# branches = ["beamEnergy","EE_Sum","FH_Sum","AH_Sum"]
branches = ["EE_Sum","FH_Sum","AH_Sum"]

def GetData(sp_val):      # train = sp_val*data , test = (1-sp_val)*data
    print ("%.1f percent of data taken for training" % (100*sp_val))
    rawdat = []
    xtrain = []
    ytrain = []
    xtest = []
    ytest = []
    energy=[]
    
    for i in ENR:
        fname = "%s%s%s" % (INPATH1,str(i),INPATH2)
        infile = ROOT.TFile(fname, "READ")
        tree = infile.Get("For_ML")

        temp=[]
        for jj in range(len(branches)):
            _temp = array.array( 'd', [0.])
            tree.SetBranchAddress(branches[jj],_temp)
            temp.append(_temp)
        
        temp_en = array.array('f',[0.])
        tree.SetBranchAddress("beamEnergy",temp_en)

        # for ii in range(300):
        for ii in range(tree.GetEntries()):
            tree.GetEntry(ii)
            temp2=[]
            for kk in range(len(branches)):
                temp2.append(temp[kk][0])
            rawdat.append(temp2)
            energy.append(temp_en)
        print("%s loaded" % fname)
    rawdat = np.array(rawdat)
    energy = np.array(energy)
    # rawdat=rawdat[:,:,0]
    
    ######## Shuffling the data
    NN = len(energy)
    id = np.arange(NN)
    np.random.shuffle(id)

    rawdat = rawdat[id]
    energy = energy[id]
    
    SPLIT = int(sp_val*NN)
    xtrain = rawdat[0:SPLIT]
    ytrain = energy[0:SPLIT]
    xtest = rawdat[SPLIT:NN]
    ytest = energy[SPLIT:NN]
    
    return [xtrain, ytrain, xtest, ytest];
