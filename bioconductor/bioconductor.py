import pandas as pd
from rpy2.robjects.packages import importr

biocPkgTools = importr('BiocPkgTools')
biocPkgList = biocPkgTools.biocPkgList()
biocPkgList = pd.DataFrame(biocPkgList)
# save each row as a separate dictionary with column names as keys
biocPkgList = biocPkgList.to_dict('records')
print(biocPkgList)
