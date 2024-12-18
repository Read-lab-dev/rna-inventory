from pyscenic import binarization
import pandas as pd
auc_mtx = pd.read_csv("reg.csv", index_col="Regulon")
# time consuming
auc_binary = binarization.binarize(auc_mtx, num_workers=8)
auc_binary[0]
auc_binary[1]
