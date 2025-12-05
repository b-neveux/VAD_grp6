import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

df = pd.read_excel("Sources/data/Aircraft_Data.xlsx")
var_x = "MTOW"
var_y = "ZFW"

X = df[var_x]
Y = df[var_y]

sns.pairplot(df.iloc[:,0:9])

plt.scatter(X, Y)
plt.xlabel(var_x)
plt.ylabel(var_y)
plt.show()