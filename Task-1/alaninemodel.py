#   alaninemodel.py
#   Input: CVS file containing amino acid counts of individual proteins
#   Output: simple linear regression model of Alanine

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

data = pd.read_csv('RAP-DB_protein.csv')

id = data["GenID"]
length = data['Length'].to_numpy()
aCount = data['A'].to_numpy()

# linear regression model
x = length.reshape(-1, 1)
y = aCount
model = LinearRegression(fit_intercept=False)   # length = 0 implies count = 0
model.fit(x, y)

r_sq = model.score(x, y)
slope = model.coef_[0]
xfit = np.linspace(0, np.amax(length), int(np.amax(length) / 10))
yfit = model.predict(xfit[:, np.newaxis])

print("Coefficient of Determination:", r_sq)
print("Slope:", slope)
print("Equation:  AlanineCount = %lf(Length)\n" % slope)

# matplotlib Scatter plot
plt.rcParams["figure.figsize"]=(15, 10)
font = {'family': 'sans serif',
        'weight': 'bold',
        'size': 30}
plt.rc('font', **font)
plt.xlabel("Length of Protein")
plt.ylabel("Alanine Count")

plt.scatter(length, aCount, marker="x")
plt.plot(xfit, yfit, color="black")

plt.show()

# find protein with maximum deviation
diff = []
for i in range(len(id)):
    predictedCount = model.predict([[length[i]]])
    diff.append(abs(predictedCount - aCount[i]))

i_max=diff.index(max(diff))
print("Protein with maximum deviation:\nId:", id[i_max], " Length:", length[i_max], " Actual Count:", aCount[i_max], " Predicted Count:", model.predict([[length[i_max]]])[0])
