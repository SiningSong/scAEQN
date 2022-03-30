import pandas as pd
import torch
import torch.nn as nn
from model.model import AutoEncoder,datasets,quant
from model.model2 import AutoEncoder,datasets,quant
from sklearn.cluster import KMeans
from sklearn import metrics
from numpy import unique
import time
from itertools import chain

start = time.time()

#load data
data_path = './data/human/count_human.csv'
batch_path = './data/human/batch_human.csv'
quant_data_path = './data/human/quant_human.csv'
nrows,ncols,data= datasets(data_path)
quant_data = quant(data_path,batch_path,method="row/column")

cell = pd.read_csv('./data/human/cell_type_human.csv',index_col=0)

# Hyper Parameters
maxite = 10000
LR= 5e-4       # learning rate
eps=7e-2
# 将dataframe数据转换为torch格式
x_data = torch.from_numpy(data.values)
y_data = torch.from_numpy(quant_data.values)
x_data = x_data.float()
y_data = y_data.float()
#n_outputs = np.shape(quant_data)[0]

#autoencoder = AutoEncoder(ncols,n_outputs,ace1=nn.Tanh(),ace2=nn.Softplus(),ace3=nn.Softplus(),acd1=nn.Softplus(),acd2=nn.Softplus(),acd3=nn.Sigmoid())
autoencoder = AutoEncoder(ncols,ace1=nn.ReLU(),ace2=nn.ReLU(),acd1=nn.ReLU(),acd2=nn.ReLU(),acd3=nn.ReLU())
optimizer = torch.optim.Adam(autoencoder.parameters(),lr=LR)
loss_func = nn.MSELoss()

for epoch in range(maxite):
    encoded, decoded,corr= autoencoder(x_data)
    loss = loss_func(corr, y_data)  # mean square error
    optimizer.zero_grad()  # clear gradients for this training step
    loss.backward()    # backpropagation, compute gradients
    optimizer.step()    # apply gradients

    if loss<eps:
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
    cluster = KMeans(n_clusters=4)
    cluster.fit(latent)
    cluster_result = cluster.labels_
    unique(cluster_result)
    ARI = ARI+metrics.adjusted_rand_score(cell['index'], cluster_result)
ARI = ARI/10
print(ARI)
end = time.time()
print('Running time: %s Seconds'%(end-start))

##保存数据
latent=pd.DataFrame(latent.T,columns=data.index)
output=pd.DataFrame(output.T,index=data.columns,columns=data.index)
latent.to_csv('E:/mymodel/methodcompare/human pancreas/humanpanc_scaeqn_latent.csv')
output.to_csv('E:/mymodel/methodcompare/human pancreas/humanpanc_scaeqn_output.csv')

#绘制PCA图
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
newData=pca.fit_transform(latent)
print(newData[:4])
x1=[n[0] for n in newData]
x2=[n[1] for n in newData]
batch = pd.read_csv(batch_path)

import matplotlib.pyplot as plt
plt.figure(figsize=(10,10))
fig1 = plt.subplot(221)
plt.title("celltype")
plt.xlabel("PC1")
plt.ylabel("PC2")
fig1.scatter(x1,x2,c=cluster_result,marker='o',label=["beta","alpha","delta","gamma"])
#fig1.legend(["beta","alpha","delta","gamma"])

fig2=plt.subplot(222)
plt.title("batch")
plt.xlabel("PC1")
plt.ylabel("PC2")
cm=plt.cm.get_cmap('RdYlBu')
fig2.scatter(x1,x2,c=batch['index'],marker='o',cmap=cm,label=['Non T2D 1', 'Non T2D 10', 'Non T2D 11', 'Non T2D 12', 'Non T2D 2',
       'Non T2D 3', 'Non T2D 4', 'Non T2D 5', 'Non T2D 6', 'Non T2D 7','Non T2D 8', 'Non T2D 9'])
#fig2.legend(['Non T2D 1', 'Non T2D 10', 'Non T2D 11', 'Non T2D 12', 'Non T2D 2',
       #'Non T2D 3', 'Non T2D 4', 'Non T2D 5', 'Non T2D 6', 'Non T2D 7','Non T2D 8', 'Non T2D 9'])
#plt.savefig("E:/mymodel/figure/mymodel_mouse.png",dpi=400)
plt.show()
