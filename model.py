import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()



def datasets(data_path):
    data = pd.read_csv(data_path,index_col=0)
    data = data.T
    nrows = np.shape(data)[0]
    ncols = np.shape(data)[1]

    return nrows,ncols,data

def quant(data_path,batch_path):    #要求，data需要没有行列标题，batch处理为数值型
    robjects.r("""
        quant <- function(data_path,batch_path){
        library(QuantNorm)
        data<-read.csv(data_path,row.names = 1)
        batch<-read.csv(batch_path)
        batch<-batch$index
        quantdata <- QuantNorm(data,batch,logdat=F,method="vectorize",cor_method='pearson',max=50)
        return(quantdata)
    }
        """)
    quantdata = robjects.r['quant'](data_path, batch_path)
    data = pd.read_csv(data_path,index_col=0)
    quant_data = pd.DataFrame(quantdata)
    quant_data = pd.DataFrame(quant_data.values.reshape(data.shape[1], data.shape[1]))
    return quant_data



class AutoEncoder(nn.Module):                    #固定套路
    def __init__(self,
                 n_inputs,
                 n_outputs,
                 ace1,
                 ace2,
                 ace3,
                 acd1,
                 acd2,
                 acd3):
        super(AutoEncoder, self).__init__()
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs
        # 定义编码器层数神经元个数及激活函数
        self.encoder = nn.Sequential(
            nn.Linear(n_inputs,512),
            ace1,
            nn.Linear(512, 128),
            ace2,
            nn.Linear(128, 64),
            ace3,
            nn.Linear(64, 12),
            )
        # 定义解码器层数神经元个数及激活函数
        self.decoder = nn.Sequential(
            nn.Linear(12, 64),
            acd1,
            nn.Linear(64, 128),
            acd2,
            nn.Linear(128, n_outputs),
            acd3,  # compress to a range (0, 1)
        )

    # 将编码和解码连接起来
    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return encoded,decoded
































