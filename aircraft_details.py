import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

db = pd.read_excel("Sources/data/Aircraft_Data.xlsx")
var_x = "MTOW"
var_y = "ZFW"

X = db[var_x]
Y = db[var_y]

plt.scatter(X, Y)
plt.xlabel(var_x)
plt.ylabel(var_y)
plt.show()