"""
@author: Radosław Pławecki
"""

import pandas as pd

filename = "./genome-data-processing/data/human_genome_filtered.csv"
data = pd.read_csv(filename, delimiter=';')
df = pd.DataFrame(data)

print(df.head)
