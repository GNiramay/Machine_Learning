# Program for linear regression using ML.
# model: RecoEn = [W1 W2 W3]*[EE FH AH]' with W1, W2, W3 chosen dynamically

# import LoadData as ld
import LoadData_Nocut as ld
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from ROOT import TFile,TTree
import array
import os
import sys
import tensorflow.compat.v1 as tf1

def Do_Regression(ENR,PhList):
    N_iter = 10                  # no. of iterations
    mu = 0.1                        # mean val. of initial weights
    sigma=0.01                       # std.dev of initial weights
    LR = 0.01                       # Learning rate
    DATE = "08_01_20"               # Date of running

    if len(ENR)==1:
        EName = "%s_E%s_NEW" % (DATE,ENR[0])
    else:
        EName = "%s_E_All_NEW" % DATE

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    
    [xtr,ytr,xtst,ytst] = ld.GetData(0.75,ENR,PhList)

    xtrain = xtr
    ytrain = ytr
    xtest = xtst
    ytest = ytst
    [ntrain,nnode] = np.shape(xtrain)
    print(ytrain[1:4])
    print(xtrain[1:4])

    with tf.device('/gpu:0'):               # For GPU
        InX = tf.placeholder(tf.float32,[None,nnode],name="InX")
        InY = tf.placeholder(tf.float32,[None,1],name="InY")
        Wt0 = tf.Variable(tf.random.normal([nnode,3],mean=mu,stddev=sigma),name="Wt0")
        Y1  = tf.matmul(InX,Wt0,name="Y1")

        # y_temp = tf.constant([1.0,1.0,1.0,1.0,1.0],dtype=tf.float32,shape=[1,5])
        y_temp = tf.constant([1.0,1.0,1.0],dtype=tf.float32,shape=[1,3])
        TRU = tf.matmul(InY,y_temp,name="TRU")
        z1 = tf.abs(TRU-Y1)
        II = tf.argmin(z1,axis=1)
        Y_out = Y1[:,II]
    
    cf = tf.reduce_sum(tf.pow((InY-Y_out),2)/InY)
    CFNAME = "Chi2"

    train_step = tf1.train.AdamOptimizer(LR).minimize(cf)
    init = tf1.global_variables_initializer()
    sess = tf1.Session()
    sess.run(init)
        
    CFVar=[]
    xdat=[]
    
    for i in range(N_iter):
        sess.run(train_step,feed_dict={InX: xtrain, InY: ytrain})
        xdat.append(i)
        CFVar.append(sess.run(cf,feed_dict={InX: xtrain, InY: ytrain}))
    plt.scatter(xdat,CFVar,s=1)
    plt.show()
    plt.savefig('bbb.png')

    WEIGHTS = sess.run(Wt0)
    print("\n",WEIGHTS,"\n")

    EvalOut = TFile("OUT_ROOT/AAA.root","RECREATE")
    ttr2  = TTree("Closure","testing the model on training data")
    y_clos = sess.run(Y_out,feed_dict={InX: xtrain,InY: ytrain})
    E_clos_tru = np.array(ytrain[0])
    E_clos_pred = np.array(y_clos[0])

    Y3 = ttr2.Branch("E_clos_tru",E_clos_tru,"E_clos_tru/F")
    Y4 = ttr2.Branch("E_clos_pred",E_clos_pred,"E_clos_pred/F")
    for ii in range(len(ytrain)):
        E_clos_tru[0] = ytrain[ii]
        E_clos_pred[0] = y_clos[ii]
        ttr2.Fill()
    ttr2.Write()
    EvalOut.Close()

Do_Regression([20],"QGSP_FTFP_BERT_EML")
