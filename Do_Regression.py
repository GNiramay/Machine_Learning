# Program to do actual regression on the available data.
# import LoadData as ld
import LoadData as ld
# import Test_Model as tm
import tensorflow as tf
# import tensorflow_probability as tfp
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from ROOT import TFile,TTree
import array
# Just disables the warning, doesn't enable AVX/FMA
import os
import sys
import tensorflow.compat.v1 as tf1

def Do_Regression(ENR):
    N_iter = 10000                  # no. of iterations
    mu = 0.1                        # mean val. of initial weights
    sigma=0.2                       # std.dev of initial weights
    LR = 0.01                       # Learning rate
    DATE = "06_12_19"               # Date of running

    if len(ENR)==1:
        EName = "%s_E%s" % (DATE,ENR[0])
    else:
        EName = "%s_E_All" % DATE

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    
    [xtrain,ytrain,xtest,ytest] = ld.GetData(0.75,ENR)
    [ntrain,nnode] = np.shape(xtrain)
    
    with tf.device('/cpu:0'):               #if one has the gpu version of tensorflow installed, one could indicate tf.device('/gpu:0') instead
        # x = tf1.placeholder(tf.float32, [None, nnode])
        x = tf1.placeholder(tf.float32, [None, nnode],name="detwise_energy")
        Wt = tf.Variable(tf.random.normal([nnode, 1], mean=mu, stddev=sigma),name="MyWeights")
        y_pred = tf.matmul(x, Wt)
        y_true = tf1.placeholder(tf.float32, [ntrain, 1])

    cost_fun = 1/tf1.math.reduce_std(y_pred)                        # RMS inverse -> just to check performance
    # cost_fun = tf.reduce_sum(tf.pow((y_true-y_pred),2)/y_true)       # Cost function -> simple chisquare
    # cost_fun = tf.reduce_sum(tf.pow(y_true-y_pred,2))/ntrain      # Cost function -> Avg. error
    # cost_fun = tf1.math.reduce_std(y_pred)/tf1.math.reduce_mean(y_pred) # Cost function -> RMS/Mean
    # # cost_fun = tf.nn.moments(y_pred,axes=0)[1]                         # Cost function -> std. dev

    
    # lmd = 0.5                  # Lagrange's multiplier lambda
    # mean1 = tf1.reduce_mean(y_pred)
    # mean2 = tf1.math.reduce_mean(tf.pow(y_pred,2))
    # cost_fun = (1+lmd)*mean2-2*mean1-lmd*tf.pow(mean1,2) # My cost function 1

    CFNAME = "InvRMS"
    # CFNAME = "Chi2"
    # CFNAME = "AvgErr"
    # CFNAME = "RMSMean"
    # CFNAME = "MyCF1"

    train_step = tf1.train.AdamOptimizer(LR).minimize(cost_fun)
    
    saver = tf1.train.Saver([Wt])
    init = tf.initialize_all_variables()
    # init = tf1.global_variables_initializer
    sess = tf1.Session()
    sess.run(init)
    
    All_chi=[]
    for i in range(N_iter+1):
        sess.run(train_step,feed_dict={x: xtrain, y_true: ytrain})
        # MyChi = sess.run(cost_fun,feed_dict={x: xtrain, y_true: ytrain})/ntrain
        # print sess.run(Wt)
        # All_chi.append(MyChi)
        if (10*i)%N_iter==0:
           print(100*i/N_iter,"%")
        # print sess.run(Wt)
    WEIGHTS = sess.run(Wt)
    # print WEIGHTS
    # with open('Weights_Chi2_%s_GeV_21_11_19.txt' % ENR[0],'w') as ff:
    # with open('OUT_TXT/Weights_Chi2_%s.txt' % EName,'w') as ff:
    with open('OUT_TXT/Weights_%s_%s.txt' % (CFNAME,EName),'w') as ff:
        # ff.write(WEIGHTS[0])
        # print >> ff , WEIGHTS
        np.savetxt(ff,WEIGHTS, fmt="%.4f")
        orig_stdout = sys.stdout
        sys.stdout = ff
        # print(WEIGHTS)
        sys.stdout = orig_stdout
    # plt.plot(All_chi)
    # plt.show()

    # saver = tf1.train.Saver()
    # save_path = saver.save(sess, "OUT_MODEL/TempModel")
    # print("Model saved in path: %s" % save_path)
    
    # for v in tf1.get_default_graph().as_graph_def().node:
    #       print v.name
    
    #### Testing the model
    #### Following step is under construction. Till then
    #### use the normal way.
    # tm.TestIt(xtest,ytest,"OUT_MODEL/TempModel.ckpt")
    
    y_eval = sess.run(y_pred,feed_dict={x: xtest})
    # ytest=np.array(ytest)
    # y_eval=np.array(y_eval)
    
    ###### Store the output of the test as root tree.
    # EvalOut = TFile("OUT_ROOT/CF_Chi2_%s_GeV_21_11_19.root" % ENR[0],"RECREATE")
    # EvalOut = TFile("OUT_ROOT/CF_Chi2_%s.root" % EName,"RECREATE")
    EvalOut = TFile("OUT_ROOT/CF_%s_%s.root" % (CFNAME,EName),"RECREATE")
    ttr  = TTree("MainTree","to check performance of ML training")
    
    temp1 = np.array(ytest[0])
    temp2 = np.array(y_eval[0])
    
    Y1 = ttr.Branch("temp1",temp1,"temp1/D")
    Y2 = ttr.Branch("temp2",temp2,"temp2/F")
    
    for ii in range(len(ytest)):
        temp1 = np.array([ytest[ii][0]])
        temp2 = np.array(y_eval[ii])
        ttr.Fill()
    ttr.Write()
    EvalOut.Close()

# Do_Regression([20])
# Do_Regression([50])
# Do_Regression([80])
# Do_Regression([100])
# Do_Regression([120])
Do_Regression([200])
# Do_Regression([250])
# Do_Regression([300])
# Do_Regression([20,50,80,100,120,200,250,300])

# os.system("root -l -b -q Driver_Script.C"); # To store histograms
