# Importing libraries

import numpy as np
import pandas as pd
import random


# Defining composition

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")


def composition(seq):
    return {aa: seq.count(aa) / (len(seq) - seq.count("U")) for aa in amino_acids}


# Loading data and computing compositions

alpha = open("alpha_40.txt", mode="r").readlines()
alpha = [x.strip() for x in alpha[1::2]]
alpha_comp = composition("".join(alpha))

beta = open("beta_40.txt", mode="r").readlines()
beta = [x.strip() for x in beta[1::2]]
beta_comp = composition("".join(beta))

overall_comp = pd.DataFrame([alpha_comp, beta_comp], index=["TMH (alpha)", "TMB (beta)"]).transpose()

alpha = pd.DataFrame(alpha, columns=["Sequence"])
beta = pd.DataFrame(beta, columns=["Sequence"])

for amino_acid in amino_acids:
    alpha[amino_acid] = alpha["Sequence"].apply(lambda x: x.count(amino_acid) / (len(x) - x.count("U")))
    beta[amino_acid] = beta["Sequence"].apply(lambda x: x.count(amino_acid) / (len(x) - x.count("U")))

# Computing Fisher discriminant ratio

alpha_stats = alpha.drop("Sequence", axis=1).describe().loc[["mean", "std"]]
beta_stats = beta.drop("Sequence", axis=1).describe().loc[["mean", "std"]]
fisher_ratios = {aa: (alpha_stats.loc["mean"][aa] - beta_stats.loc["mean"][aa]) ** 2 / (alpha_stats.loc["std"][aa] ** 2 + beta_stats.loc["std"][aa] ** 2) for aa in amino_acids}
fisher_ratios = pd.DataFrame(fisher_ratios, index=["FDR"]).transpose().sort_values("FDR", ascending=False)

# Performing the discrimination


def predict_class(dev1, dev2):
    if dev1 > dev2:
        return "Beta"
    else:
        return "Alpha"


predict_class = np.vectorize(predict_class)

alpha_alpha_devs = []
alpha_beta_devs = []
for i in range(len(alpha)):
    alpha_dev = 0
    beta_dev = 0
    for aa in amino_acids:
        alpha_dev += abs(alpha.iloc[i][aa] - alpha_comp[aa])
        beta_dev += abs(alpha.iloc[i][aa] - beta_comp[aa])
    alpha_alpha_devs.append(alpha_dev)
    alpha_beta_devs.append(beta_dev)
alpha_devs = pd.DataFrame([alpha_alpha_devs, alpha_beta_devs]).transpose()
alpha_devs.columns = ["Alpha dev", "Beta dev"]
alpha = pd.concat([alpha, alpha_devs], axis=1)
alpha["Prediction"] = predict_class(alpha["Alpha dev"], alpha["Beta dev"])
tp = alpha["Prediction"].value_counts()["Alpha"]
fn = alpha["Prediction"].value_counts()["Beta"]

beta_alpha_devs = []
beta_beta_devs = []
for i in range(len(beta)):
    alpha_dev = 0
    beta_dev = 0
    for aa in amino_acids:
        alpha_dev += abs(beta.iloc[i][aa] - alpha_comp[aa])
        beta_dev += abs(beta.iloc[i][aa] - beta_comp[aa])
    beta_alpha_devs.append(alpha_dev)
    beta_beta_devs.append(beta_dev)
beta_devs = pd.DataFrame([beta_alpha_devs, beta_beta_devs]).transpose()
beta_devs.columns = ["Alpha dev", "Beta dev"]
beta = pd.concat([beta, beta_devs], axis=1)
beta["Prediction"] = predict_class(beta["Alpha dev"], beta["Beta dev"])
tn = beta["Prediction"].value_counts()["Beta"]
fp = beta["Prediction"].value_counts()["Alpha"]

# Computing the performance

sensitivity = tp / (tp + fn)
specificity = tn / (tn + fp)
accuracy = (tp + tn) / (tp + tn + fp + fn)

print("\n100% TRAIN TEST")
print(f"TP: {tp}")
print(f"TN: {tn}")
print(f"FP: {fp}")
print(f"FN: {fn}")
print(f"Sensitivity: {sensitivity}")
print(f"Specificity: {specificity}")
print(f"Accuracy: {accuracy}")

# Repeating for a 50% train-test split

alpha = open("alpha_40.txt", mode="r").readlines()
alpha = [x.strip() for x in alpha[1::2]]
beta = open("beta_40.txt", mode="r").readlines()
beta = [x.strip() for x in beta[1::2]]

random.seed(4110)
train_alpha = random.choices(alpha, k=int(len(alpha) * 0.5))
test_alpha = list(set(alpha) - set(train_alpha))
train_beta = random.choices(beta, k=int(len(beta) * 0.5))
test_beta = list(set(beta) - set(train_beta))

alpha_comp = composition("".join(train_alpha))
beta_comp = composition("".join(train_beta))

test_alpha = pd.DataFrame(test_alpha, columns=["Sequence"])
test_beta = pd.DataFrame(test_beta, columns=["Sequence"])

for amino_acid in amino_acids:
    test_alpha[amino_acid] = test_alpha["Sequence"].apply(lambda x: x.count(amino_acid) / (len(x) - x.count("U")))
    test_beta[amino_acid] = test_beta["Sequence"].apply(lambda x: x.count(amino_acid) / (len(x) - x.count("U")))

alpha_alpha_devs = []
alpha_beta_devs = []
for i in range(len(test_alpha)):
    alpha_dev = 0
    beta_dev = 0
    for aa in amino_acids:
        alpha_dev += abs(test_alpha.iloc[i][aa] - alpha_comp[aa])
        beta_dev += abs(test_alpha.iloc[i][aa] - beta_comp[aa])
    alpha_alpha_devs.append(alpha_dev)
    alpha_beta_devs.append(beta_dev)
alpha_devs = pd.DataFrame([alpha_alpha_devs, alpha_beta_devs]).transpose()
alpha_devs.columns = ["Alpha dev", "Beta dev"]
test_alpha = pd.concat([test_alpha, alpha_devs], axis=1)
test_alpha["Prediction"] = predict_class(test_alpha["Alpha dev"], test_alpha["Beta dev"])
tp = test_alpha["Prediction"].value_counts()["Alpha"]
fn = test_alpha["Prediction"].value_counts()["Beta"]

beta_alpha_devs = []
beta_beta_devs = []
for i in range(len(test_beta)):
    alpha_dev = 0
    beta_dev = 0
    for aa in amino_acids:
        alpha_dev += abs(test_beta.iloc[i][aa] - alpha_comp[aa])
        beta_dev += abs(test_beta.iloc[i][aa] - beta_comp[aa])
    beta_alpha_devs.append(alpha_dev)
    beta_beta_devs.append(beta_dev)
beta_devs = pd.DataFrame([beta_alpha_devs, beta_beta_devs]).transpose()
beta_devs.columns = ["Alpha dev", "Beta dev"]
test_beta = pd.concat([test_beta, beta_devs], axis=1)
test_beta["Prediction"] = predict_class(test_beta["Alpha dev"], test_beta["Beta dev"])
tn = test_beta["Prediction"].value_counts()["Beta"]
fp = test_beta["Prediction"].value_counts()["Alpha"]

sensitivity = tp / (tp + fn)
specificity = tn / (tn + fp)
accuracy = (tp + tn) / (tp + tn + fp + fn)

print("\n50% TRAIN TEST")
print(f"TP: {tp}")
print(f"TN: {tn}")
print(f"FP: {fp}")
print(f"FN: {fn}")
print(f"Sensitivity: {sensitivity}")
print(f"Specificity: {specificity}")
print(f"Accuracy: {accuracy}")
