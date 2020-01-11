import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt

import LoadData_3_Inputs as ld # To get the input data points

import ROOT
import keras
from keras.models import Sequential
from keras.layers import Dense, Activation
from ROOT import TFile,TTree
import array
import sys
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' # to supress unwanted errors
os.environ['CUDA_VISIBLE_DEVICES'] = '-1' # To force GPU usage

def ML_Regression(ENR):
    DATE = "aaaa_16_12_19"
    if len(ENR)==1:
        EName = "%s_E%s" % (DATE,ENR[0])
    else:
        EName = "%s_E_All" % DATE
    [xtrain,ytrain,xtest,ytest] = ld.GetData(0.75,ENR,'QGSP_FTFP_BERT_EML')
    [ntrain,nnode] = np.shape(xtrain)

    ############ ML model setup
    model = Sequential()
    model.add(Dense(nnode, input_dim=nnode, activation='linear'))
    model.add(Dense(1, input_dim=nnode, activation='linear'))
    model.summary()
    
    ############ ML Training
    CFNAME = "mse"                  # cost function-> min. standard error
    model.compile(optimizer='adam', loss='mse', metrics=['mse'])
    model.fit(xtrain,ytrain, batch_size=100, epochs=20, verbose=0, shuffle=True)
    
    # ############ Get weights
    # WEIGHTS = model.layers[0].get_weights()
    # weight2 = np.array(WEIGHTS)
    # # with open('OUT_TXT/Weights_%s_%s.txt' % (CFNAME,EName),'w') as ff:
    # #     np.savetxt(ff,weight2[0],fmt="%.4f")
    # # print(weight2)
            
    ############ ML Testing
    predict = model.predict(xtest)
    
    # EvalOut = TFile("OUT_ROOT/CF_%s_%s.root" % (CFNAME,EName),"RECREATE")
    EvalOut = TFile("OUT_ROOT/QGSP_FTFP_BERT_EML_3_Inputs_EAll.root","RECREATE")
    ttr  = TTree("MainTree","to check performance of ML training")
        
    BeamEn = np.array(ytest[0])
    RecoEn = np.array(predict[0])
        
    Y1 = ttr.Branch("BeamEn",BeamEn,"BeamEn/F")
    Y2 = ttr.Branch("RecoEn",RecoEn,"RecoEn/F")
        
    for ii in range(len(ytest)):
        BeamEn[0] = ytest[ii]
        RecoEn[0] = predict[ii]
        ttr.Fill()
    ttr.Write()
    EvalOut.Close()

# ML_Regression([20])
# ML_Regression([50])
# ML_Regression([80])
# ML_Regression([100])
# ML_Regression([120])
# ML_Regression([200])
# ML_Regression([250])
# ML_Regression([300])
ML_Regression([20,50,80,100,120,200,250,300])
