import json
import pandas as pd
from rpy2.robjects.packages import importr

biocPkgTools = importr('BiocPkgTools')
biocPkgList = biocPkgTools.biocPkgList()
df = pd.DataFrame(biocPkgList)

print(df.get_lo)
