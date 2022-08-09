from rpy2.robjects.vectors import StrVector
import rpy2.robjects.packages as rpackages
utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)
packages = ('Biostrings', 'dyplr', 'tidyverse')
utils.install_packages(StrVector(packages))


# def String_to_DNA(x):
#     return Biostrings.DNAString(x)


# omicron['dna'] = omicron.apply(lambda row: String_to_DNA(row['spike']), axis=1)
