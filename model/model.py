import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import rpy2.robjects as robjects
import time

from rpy2.robjects.packages import importr
import torch.nn.functional as F


def datasets(data_path):
    data = pd.read_csv(data_path, index_col=0)
    data = data.T
    nrows = np.shape(data)[0]
    ncols = np.shape(data)[1]
    return nrows, ncols, data


def quant(data_path, batch_path,method):
    robjects.r("""
        quant <- function(data_path,batch_path,method){
        library(QuantNorm)
        data<-read.csv(data_path,row.names = 1)
        batch<-read.csv(batch_path)
        batch<-batch$index
        quantdata <- QuantNorm(data,batch,logdat=F,method=method,cor_method='pearson',max=50)
        return(quantdata)
    }
        """)
    quantdata = robjects.r['quant'](data_path, batch_path,method)
    data = pd.read_csv(data_path, index_col=0)
    quant_data = pd.DataFrame(quantdata)
    quant_data = pd.DataFrame(quant_data.values.reshape(data.shape[1], data.shape[1]))
    return quant_data


def corrcoef(x):
    """
    定义tensor中的Pearson相关系数，按行求
    :param x: x为tensor格式的数据集
    :return: x的pearson相关系数矩阵，大小为x的行数*行数
    """
    x_mean=torch.mean(x,axis=1)
    x_reducemean=x-x_mean.reshape(x.shape[0],1)
    frac = x.shape[1] - 1
    numerator = torch.matmul(x_reducemean, x_reducemean.T)
    x_var = torch.var(x, axis=1, unbiased=True).reshape(x.shape[0], 1)
    denominator = torch.sqrt(torch.matmul(x_var, x_var.T)) * frac
    corr = numerator / denominator
    return corr


# 定义相关系数层
class corrlayer(nn.Module):
    def __init__(self):
        # 定义参数
        super(corrlayer, self).__init__()

    # 定义前向传播
    def forward(self, input):
        result = 1 - corrcoef(input)
        return result


# 定义模型
class AutoEncoder(nn.Module):  # 固定套路
    def __init__(self,
                 n_inputs,
                 ace1,
                 ace2,
                 acd1,
                 acd2,
                 acd3
                 ):
        super(AutoEncoder, self).__init__()
        self.n_inputs = n_inputs
        # 定义编码器层数神经元个数及激活函数
        self.encoder = nn.Sequential(
            nn.Linear(n_inputs, 128),
            ace1,
            nn.Linear(128,64),
            ace2,
            nn.Linear(64, 12),
        )
        # 定义解码器层数神经元个数及激活函数
        self.decoder = nn.Sequential(
            nn.Linear(12, 64),
            acd1,
            nn.Linear(64,128),
            acd2,
            nn.Linear(128, n_inputs),
            acd3
        )
        # 模型中加入相关系数
        self.corrlayer = corrlayer()

    # 将编码和解码连接起来
    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        corr = self.corrlayer(decoded)
        return encoded, decoded, corr



