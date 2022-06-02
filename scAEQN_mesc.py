import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from numpy import unique
import time
from itertools import chain
from model.model import AutoEncoder,datasets,quant
from sklearn.cluster import KMeans
from sklearn import metrics

start = time.time()
#load data
data_path = './data/mesc/count_mesc.csv'  # 38616*469
batch_path = './data/mesc/batch_mesc.csv'
nrows,ncols,data= datasets(data_path)

cell = pd.read_csv('./data/mesc/cell_type_mesc.csv',index_col=0)
quant_data  = quant(data_path,batch_path,method="row/column")
#quant_data = pd.read_csv('./data/mesc/quant_mesc.csv',header=None)

# Hyper Parameters
EPOCH= 10000
#BATCH_SIZE= 256
LR= 8e-5         # learning rate
eps = 2e-3

# 将dataframe数据转换为torch格式
x_data = torch.from_numpy(data.values)
y_data = torch.from_numpy(quant_data.values)
x_data = x_data.float()
y_data = y_data.float()

n_outputs = np.shape(quant_data)[0]

autoencoder = AutoEncoder(ncols,ace1=nn.ReLU(),ace2=nn.ReLU(),acd1=nn.ReLU(),acd2=nn.ReLU(),acd3=nn.ReLU())
optimizer = torch.optim.Adam(autoencoder.parameters(),lr=LR)
loss_func = nn.MSELoss()

for epoch in range(EPOCH):
    encoded, decoded,corr= autoencoder(x_data)
    loss = loss_func(corr, y_data)  # mean square error
    optimizer.zero_grad()  # clear gradients for this training step
    loss.backward()    # backpropagation, compute gradients
    optimizer.step()    # apply gradients
#xunliann
    if loss < eps:
        print('##############################')
        print('STOP EPOCH:', epoch, '| train loss: %.8f' % loss.data.numpy())
        break
    elif epoch % 5 == 0:
        print('Epoch: ', epoch, '| train loss: %.8f' % loss.data.numpy())

#得到潜在层的特征
latent = encoded.detach().numpy()
output = decoded.detach().numpy()

#聚类
ARI = 0
for i in range(10):
    cluster = KMeans(n_clusters=3)
    cluster.fit(latent)
    cluster_result = cluster.labels_
    unique(cluster_result)
    ARI = ARI+metrics.adjusted_rand_score(cell['label'], cluster_result)

ARI = ARI/10
print(ARI)
end = time.time()
print('Running time: %s Seconds'%(end-start))

##保存数据
latent=pd.DataFrame(latent.T,columns=data.index)
output=pd.DataFrame(output.T,index=data.columns,columns=data.index)
latent.to_csv('E:/mymodel/methodcompare/mesc/mesc_scaeqn_latent.csv')
output.to_csv('E:/mymodel/methodcompare/mesc/mesc_scaeqn_output.csv')



