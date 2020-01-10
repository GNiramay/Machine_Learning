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
    N_iter = 100                  # no. of iterations
    # N_iter = 10000                  # no. of iterations
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

    # print("Shower start in %s" % Detector[det])
    xtrain = xtr
    ytrain = ytr
    xtest = xtst
    ytest = ytst
    [ntrain,nnode] = np.shape(xtrain)
    print(ytrain[1:4])
    print(xtrain[1:4])
        
    tempwt=[[0.01,0.05,0.02],[0,0.05,0.02],[0,0,0.01]]



    with tf.device('/gpu:0'):               # For GPU
        InX = tf.placeholder(tf.float32,[None,nnode],name="InX")
        InY = tf.placeholder(tf.float32,[None,1],name="InY")
        Wt0 = tf.Variable(tf.random.normal([nnode,3],mean=mu,stddev=sigma),name="Wt0")
        # Wt0 = tf.Variable(tempwt,shape=[nnode,3],dtype=tf.float32,name="Wt0")
        Y1  = tf.matmul(InX,Wt0,name="Y1")

        # y_temp = tf.constant([1.0,1.0,1.0,1.0,1.0],dtype=tf.float32,shape=[1,5])
        y_temp = tf.constant([1.0,1.0,1.0],dtype=tf.float32,shape=[1,3])
        TRU = tf.matmul(InY,y_temp,name="TRU")

        z1 = tf.abs(tf.math.subtract(Y1,TRU))
        z2 = tf.negative(z1)
        Wt1 = tf.exp(z2,name="Wt1")
        # Wtemp = tf.exp(z2,name="Wt1")
        # Wt1 = tf.nn.relu(Wtemp)

        Y2  = tf.multiply(Y1,Wt1,name="Y2")
        Y_out= tf.matmul(Y2,y_temp,transpose_b=True,name="Y_out")

    print("\n",InX,"\n",InY,"\n",Wt0,"\n",Y1,"\n",TRU,"\n",Wt1,"\n",Y2,"\n",Y_out,"\n")
    
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

        # YO = sess.run(Y_out,feed_dict={InX: xtrain, InY: ytrain})
        # YOne = sess.run(Y1,feed_dict={InX: xtrain, InY: ytrain})
        # TRU_ = sess.run(TRU,feed_dict={InX: xtrain, InY: ytrain})
        # Zone = sess.run(z1,feed_dict={InX: xtrain, InY: ytrain})
        # Ztwo = sess.run(z2,feed_dict={InX: xtrain, InY: ytrain})
        # Wtone = sess.run(Wt1,feed_dict={InX: xtrain, InY: ytrain})
        # Ytwo = sess.run(Y2,feed_dict={InX: xtrain, InY: ytrain})

        # print("Y_out\n",YO)
        # print("YOne\n",YOne)
        # print("TRU\n",TRU_)
        # print("z1\n",Zone)
        # print("z2\n",Ztwo)
        # print("Wtone\n",Wtone)
        # print("Ytwo\n",Ytwo)

        # # print("Y_out\n",Y_out.eval())
        # # print("YOne\n",Y1.eval())
        # # print("TRU\n",TRU.eval())
        # # print("z1\n",z1.eval())
        # # print("z2\n",z2.eval())
        # # print("Wtone\n",Wt1.eval())
        # # print("Ytwo\n",Y2.eval())

    plt.scatter(xdat,CFVar,s=1)
    plt.show()
    plt.savefig('bbb.png')

    WEIGHTS = sess.run(Wt0)
    print("\n",WEIGHTS,"\n")
    # y_eval = sess.run(Y_out,feed_dict={InX: xtest,  InY: ytest})
    # print(np.shape(y_eval))
    # print(y_eval[1:4])
    # with open("OUT_TXT/aaa.txt",'w') as ff:
        # np.savetxt(ff,,fmt="%f")

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
        # if ii<10:
        #     print(y_clos[ii])
    ttr2.Write()
    EvalOut.Close()

Do_Regression([20],"QGSP_FTFP_BERT_EML")
