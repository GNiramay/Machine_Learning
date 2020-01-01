# Program for linear regression using ML.
import LoadData as ld
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from ROOT import TFile,TTree
import array
import os
import sys
import tensorflow.compat.v1 as tf1

def Do_Regression(ENR):
    N_iter = 10000                  # no. of iterations
    mu = 0.1                        # mean val. of initial weights
    sigma=0.2                       # std.dev of initial weights
    LR = 0.01                       # Learning rate
    DATE = "26_12_19"               # Date of running

    if len(ENR)==1:
        EName = "%s_E%s" % (DATE,ENR[0])
    else:
        EName = "%s_E_All" % DATE

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    
    [xtr,ytr,xtst,ytst] = ld.GetData(0.75,ENR)
    Detector=["EE","FH"]

    for det in [0,1]:
        print("Shower start in %s" % Detector[det])
        xtrain = xtr[det]
        ytrain = ytr[det]
        xtest = xtst[det]
        ytest = ytst[det]
        [ntrain,nnode] = np.shape(xtrain)
        # with tf.device('/cpu:0'):               # For CPU
        with tf.device('/gpu:0'):               # For GPU
            x = tf1.placeholder(tf.float32, [None, nnode],name="detwise_energy")
            Wt = tf.Variable(tf.random.normal([nnode, 1], mean=mu, stddev=sigma),name="MyWeights")
            y_pred = tf.matmul(x, Wt)
            y_true = tf1.placeholder(tf.float32, [ntrain, 1])

        # cost_fun = tf.reduce_sum(tf.pow((y_true-y_pred),2)/y_true)       # Cost function -> simple chisquare
        # cost_fun = tf1.math.reduce_std(y_pred)/tf1.math.reduce_mean(y_pred) # Cost function -> RMS/Mean

        # CFNAME = "Chi2"
        # CFNAME = "RMSMean"

        cf1 = tf1.math.reduce_std(y_pred)
        cf2 = tf.reduce_sum(tf.pow((y_true-y_pred),2)/y_true)
        CFNAME = "RMS_Chi2"

        train_step_1 = tf1.train.AdamOptimizer(LR).minimize(cf1)
        train_step_2 = tf1.train.AdamOptimizer(LR).minimize(cf2)

        # saver = tf1.train.Saver([Wt])
        init = tf.initialize_all_variables()
        sess = tf1.Session()
        sess.run(init)
        
        CFVar=[]
        xdat=[]
        cccc2=[]
        xx2=[]
        for i in range(N_iter+1):
            sess.run(train_step_1,feed_dict={x: xtrain, y_true: ytrain})
            sess.run(train_step_2,feed_dict={x: xtrain, y_true: ytrain})

            # TempCF = sess.run(cost_fun,feed_dict={x: xtrain, y_true: ytrain})/ntrain
            # CFVar.append(TempCF)
            xdat.append(i)
            if i>3*N_iter/4:
                cccc2.append(sess.run(cf2,feed_dict={x: xtrain, y_true: ytrain}))
                xx2.append(i)
            if (10*i)%N_iter==0:
                print(100*i/N_iter,"%")
                tt = "%.3f\t%.3f" % (sess.run(cf1,feed_dict={x: xtrain, y_true: ytrain}),sess.run(cf2,feed_dict={x: xtrain, y_true: ytrain}))
                print(tt)

        WEIGHTS = sess.run(Wt)

        with open('OUT_TXT/Weights_%s_%s_%s.txt' % (CFNAME,EName,Detector[det]),'w') as ff:
            np.savetxt(ff,WEIGHTS, fmt="%.4f")
        # plt.scatter(xdat,CFVar,s=1)
        plt.scatter(xx2,cccc2)
        plt.show()
        plt.savefig('PNG/LossFunctionVar_%s_%s_%s.png' % (CFNAME,EName,Detector[det]))
        
        y_eval = sess.run(y_pred,feed_dict={x: xtest})
        y_clos = sess.run(y_pred,feed_dict={x: xtrain})

        # WEE = tf.math.scalar_mul(WEIGHTS[0],xtest[:,0])
        # WFH = tf.math.scalar_mul(WEIGHTS[1],xtest[:,1])
        # WAH = tf.math.scalar_mul(WEIGHTS[2],xtest[:,2])

        temp_W1EE = np.array(WEIGHTS[0]*xtest[0,0]) # W1*EE
        temp_W2FH = np.array(WEIGHTS[1]*xtest[0,1]) # W2*FH
        temp_W3AH = np.array(WEIGHTS[2]*xtest[0,2]) # W3*AH
        
        ###### Store the output of the test as root tree.
        EvalOut = TFile("OUT_ROOT/CF_%s_%s_%s_TwoStep.root" % (CFNAME,EName,Detector[det]),"RECREATE")
        ttr1  = TTree("Test","testing the model on unknown test data")
        ttr2  = TTree("Closure","testing the model on training data")
        
        E_test_tru = np.array(ytest[0])
        E_test_pred = np.array(y_eval[0])
        E_clos_tru = np.array(ytrain[0])
        E_clos_pred = np.array(y_clos[0])
        
        Y1 = ttr1.Branch("E_test_tru",E_test_tru,"E_test_tru/F")
        Y2 = ttr1.Branch("E_test_pred",E_test_pred,"E_test_pred/F")
        En1= ttr1.Branch("W1EE",temp_W1EE,"W1EE/F")
        En2= ttr1.Branch("W2FH",temp_W2FH,"W2FH/F")
        En3= ttr1.Branch("W3AH",temp_W3AH,"W3AH/F")

        Y3 = ttr2.Branch("E_clos_tru",E_clos_tru,"E_clos_tru/F")
        Y4 = ttr2.Branch("E_clos_pred",E_clos_pred,"E_clos_pred/F")
    
        for ii in range(len(ytest)):
            E_test_tru[0] = ytest[ii]
            E_test_pred[0] = y_eval[ii]
            temp_W1EE[0] = WEIGHTS[0]*xtest[ii,0]
            temp_W2FH[0] = WEIGHTS[1]*xtest[ii,1]
            temp_W3AH[0] = WEIGHTS[2]*xtest[ii,2]
            ttr1.Fill()
        ttr1.Write()

        for ii in range(len(ytrain)):
            E_clos_tru[0] = ytrain[ii]
            E_clos_pred[0] = y_clos[ii]
            ttr2.Fill()
        ttr2.Write()

        EvalOut.Close()

Do_Regression([20])
Do_Regression([50])
Do_Regression([80])
Do_Regression([100])
Do_Regression([120])
Do_Regression([200])
Do_Regression([250])
Do_Regression([300])
# # # Do_Regression([20,50,80,100,120,200,250,300])
