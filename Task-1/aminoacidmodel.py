#   aminoacidmodel.py
#   Input: CVS file containing amino acid counts of individual proteins
#   Output: coefficient of determination of linear models of each amino acids

from typing import List, Tuple, Any
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

r_sq: list[tuple[Any, float]] = []

data = pd.read_csv('RAP-DB_protein.csv')
length = data['Length'].to_numpy().reshape(-1, 1)
model = LinearRegression(fit_intercept=False)   # length = 0 implies count = 0

cols = ["GenID", "Length"]
for column in data.drop(cols, axis=1):
	count = data[column]
	model.fit(length, count)
	r_sq.append((column, model.score(length, count)))
r_sq.sort(key=lambda x: x[1])
r_sq.reverse()
print(r_sq)
