import json
import pandas as pd
from rpy2.robjects.packages import importr
import rpy2.robjects as ro


biocPkgTools = importr('BiocPkgTools')
biocPkgList = biocPkgTools.biocPkgList()
# save bioctools package list to a csv file
# convert to pandas dataframe
biocPkgList = pd.DataFrame(biocPkgList)
biocPkgList.to_csv('biocPkgList.csv')
