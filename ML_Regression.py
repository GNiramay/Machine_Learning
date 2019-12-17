import tensorflow as tf
import numpy as np
# import matplotlib.pyplot as plt
import LoadData as ld
import ROOT
import keras
from keras.models import Sequential
from keras.layers import Dense, Activation
from ROOT import TFile,TTree
import array
import sys
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' # to supress unwanted errors

def ML_Regression(ENR):
    DATE = "16_12_19"
    if len(ENR)==1:
        EName = "%s_E%s" % (DATE,ENR[0])
    else:
        EName = "%s_E_All" % DATE
    [xtrain,ytrain,xtest,ytest] = ld.GetData(0.75,ENR)
    [ntrain,nnode] = np.shape(xtrain)

    # with tf.device('/gpu:0'):
    ############ ML model setup
    model = Sequential()
    model.add(Dense(1, input_dim=nnode, activation='linear', use_bias=False))
    model.summary()
    
    ############ ML Training
    CFNAME = "mse"                  # cost function-> min. standard error
    model.compile(optimizer='adam', loss='mse', metrics=['mse'])
    model.fit(xtrain,ytrain, batch_size=100, epochs=100, shuffle=False, verbose=0)
    
    ############ Get weights
    WEIGHTS = model.layers[0].get_weights()
    weight2 = np.array(WEIGHTS)
    with open('OUT_TXT/Weights_%s_%s.txt' % (CFNAME,EName),'w') as ff:
        np.savetxt(ff,weight2[0],fmt="%.4f")
            
    ############ ML Testing
    predict = model.predict(xtest)
    
    EvalOut = TFile("OUT_ROOT/CF_%s_%s.root" % (CFNAME,EName),"RECREATE")
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


ML_Regression([20])
# ML_Regression([50])
# ML_Regression([80])
# ML_Regression([100])
# ML_Regression([120])
# ML_Regression([200])
# ML_Regression([250])
# ML_Regression([300])
