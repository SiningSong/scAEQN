import pandas as pd # 用于处理表格数据
import numpy as np # 用于科学计算
import matplotlib.pyplot as plt  # 绘图的核心库
from matplotlib.font_manager import FontProperties # 字体属性管理器，知道就好

data_path = 'E:/mymodel/methodcompare/results/'
data = 'mesc/'
metric = 'ARISampled_CT'
data = pd.read_csv(data_path+data+metric+'ARI_Sampledmesc_allmethod.csv')
