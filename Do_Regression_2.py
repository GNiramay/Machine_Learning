# Program to do actual regression on the available data.
# import LoadData as ld
import LoadData_2 as ld
import Test_Model as tm
import tensorflow as tf
# import tensorflow_probability as tfp
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from ROOT import TFile,TTree
import array
# Just disables the warning, doesn't enable AVX/FMA
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

[xtrain,ytrain,xtest,ytest] = ld.GetData(0.75)
[ntrain,nnode] = np.shape(xtrain)

N_iter = 100                    # no. of iterations
mu = 0.1                        # mean val. of initial weights
sigma=0.2                       # std.dev of initial weights
LR = 0.01                       # Learning rate

with tf.device('/cpu:0'):               #if one has the gpu version of tensorflow installed, one could indicate tf.device('/gpu:0') instead
    # x = tf.compat.v1.placeholder(tf.float32, [None, nnode])
    x = tf.compat.v1.placeholder(tf.float32, [None, nnode],name="detwise_energy")
    Wt = tf.Variable(tf.random.normal([nnode, 1], mean=mu, stddev=sigma),name="MyWeights")
    y_pred = tf.matmul(x, Wt)
    y_true = tf.compat.v1.placeholder(tf.float32, [ntrain, 1])

cost_fun = tf.reduce_sum(tf.pow(y_true-y_pred,2))/ntrain      # Cost function -> simple chisquare
# cost_fun = tf.nn.moments(y_pred,axes=0)[1]                         # Cost function -> std. dev
# cost_fun = tf.compat.v1.math.reduce_sum(y_pred)

train_step = tf.compat.v1.train.AdamOptimizer(LR).minimize(cost_fun)

saver = tf.compat.v1.train.Saver([Wt])
init = tf.initialize_all_variables()
# init = tf.compat.v1.global_variables_initializer
sess = tf.compat.v1.Session()
sess.run(init)

All_chi=[]
for i in range(N_iter+1):
    sess.run(train_step,feed_dict={x: xtrain, y_true: ytrain})
    MyChi = sess.run(cost_fun,feed_dict={x: xtrain, y_true: ytrain})/ntrain
    # print sess.run(Wt)
    All_chi.append(MyChi)
print sess.run(Wt)
# plt.plot(All_chi)
# plt.show()

# saver = tf.compat.v1.train.Saver()
# save_path = saver.save(sess, "OUT_MODEL/TempModel")
# print("Model saved in path: %s" % save_path)

# for v in tf.compat.v1.get_default_graph().as_graph_def().node:
#       print v.name

#### Testing the model
#### Following step is under construction. Till then
#### use the normal way.
# tm.TestIt(xtest,ytest,"OUT_MODEL/TempModel.ckpt")

y_eval = sess.run(y_pred,feed_dict={x: xtest})
# ytest=np.array(ytest)
# y_eval=np.array(y_eval)

###### Store the output of the test as root tree.
EvalOut = TFile("OUT_ROOT/CF_Chi2.root","RECREATE")
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
