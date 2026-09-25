import pandas as pd
import matplotlib.pyplot as plt

mut21 = pd.read_csv("results/mutations/barcode21.vcf", comment="#", sep="\t", header=None)
mut22 = pd.read_csv("results/mutations/barcode22.vcf", comment="#", sep="\t", header=None)

mut21['pos'] = mut21[1]
mut22['pos'] = mut22[1]

shared = set(mut21['pos']).intersection(set(mut22['pos']))
unique21 = len(set(mut21['pos']) - set(mut22['pos']))
unique22 = len(set(mut22['pos']) - set(mut21['pos']))

plt.bar(["Unique 21", "Shared", "Unique 22"], [unique21, len(shared), unique22],
        color=["skyblue", "orange", "lightgreen"])
plt.title("Mitochondrial mutation comparison: G1 vs G2")
plt.ylabel("Variant count")
plt.show()
