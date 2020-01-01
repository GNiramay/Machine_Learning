# Program to read input data and split it into training and testing datasets. -> New approach
import ROOT
import array
import numpy as np

# ENR=[20,50,80,100,120,200,250,300]
# ENR=[20]

# ENR = tuple storing which energies to use
# train = sp_val*data , test = (1-sp_val)*data
def GetData(sp_val,ENR):
    # INPATH1 ="IN_ROOT/Pion_22_11_19_Sim_"
    INPATH1 ="IN_ROOT/Pion_03_12_19_Sim_"

    INPATH2 ="_ForML.root"
    
    # branches = ["beamEnergy","EE_Sum","FH_Sum","AH_Sum"]
    branches = ["EE_Sum","FH_Sum","AH_Sum"]
    # branches = ["En_"+str(i) for i in range(79)]
    # branches = ["En_"+str(i) for i in range(1,41)]

    print ("%.1f percent of data taken for training" % (100*sp_val))

    # rawdat = []
    # energy =[]

    rawdat1 = []
    xtrain1 = []
    ytrain1 = []
    xtest1 = []
    ytest1 = []
    energy1 =[]
    
    rawdat2 = []
    xtrain2 = []
    ytrain2 = []
    xtest2 = []
    ytest2 = []
    energy2 =[]

    for i in ENR:
        Q=[0.0,0.0,0.0]
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
            ######## Classifying evets
            if GetDet(temp2)==1:
                rawdat1.append(temp2)
                energy1.append(temp_en)
                Q[1]=Q[1]+1
            elif GetDet(temp2)==2:
                rawdat2.append(temp2)
                energy2.append(temp_en)
                Q[2]=Q[2]+1
            else:
                Q[0]=Q[0]+1
            ######## End
        print("%s loaded" % fname)
        print("Beam Energy: %.f GeV" % i)
        print("mips fraction: %.3f" % (Q[0]/(tree.GetEntries())))
        print("EE shower start fraction: %.3f" % (Q[1]/tree.GetEntries()))
        print("FH shower start fraction: %.3f" % (Q[2]/tree.GetEntries()))

    rawdat1 = np.array(rawdat1)
    energy1 = np.array(energy1)

    rawdat2 = np.array(rawdat2)
    energy2 = np.array(energy2)

    # rawdat=rawdat[:,:,0]
    
    ######## Shuffling the data
    NN1 = len(energy1)
    id1 = np.arange(NN1)
    np.random.shuffle(id1)

    rawdat1 = rawdat1[id1]
    energy1 = energy1[id1]
    
    SPLIT1 = int(sp_val*NN1)
    xtrain1 = rawdat1[0:SPLIT1]
    ytrain1 = energy1[0:SPLIT1]
    xtest1 = rawdat1[SPLIT1:NN1]
    ytest1 = energy1[SPLIT1:NN1]
    ############### For FH
    NN2 = len(energy2)
    id2 = np.arange(NN2)
    np.random.shuffle(id2)

    rawdat2 = rawdat2[id2]
    energy2 = energy2[id2]
    
    SPLIT2 = int(sp_val*NN2)
    xtrain2 = rawdat2[0:SPLIT2]
    ytrain2 = energy2[0:SPLIT2]
    xtest2 = rawdat2[SPLIT2:NN2]
    ytest2 = energy2[SPLIT2:NN2]
    
    return [[xtrain1,xtrain2], [ytrain1,ytrain2], [xtest1,xtest2], [ytest1,ytest2]];

# EData : array with 3 elements: [EE,FH,AH]
# This will return 0,1,2.
# 0 -> MIPs in EE and FH -> should throw this event
# 1 -> SS in EE
# 2 -> SS in FH
def GetDet(EData):
    EnEE = EData[0]
    EnFH = EData[1]
    if EnEE>100:
        return 1
    elif EnFH>60:
        return 2
    else:
        return 0
