# Hemolysis_Correlation
This script was used by the paper Multi-omics reveals early hemolysis of red blood cells during storage altered by heritable differences in GPx4 abundance.

# Installation
This python script is known to run in python 3.7.3 with scipy installed. If you would like to run the script, please install python along with SciPy following the instructions from their respective repositories:

https://www.python.org/<br>
https://www.scipy.org/install.html

# Running
This script needs to be run in the same folder as the 2 files that it is correlating (Hemolysis.csv and Proteins.csv). I have included these two files in the repository for examination. When finished, a new folder will appear (Outputs). Within this folder, there are 3 files labeled Data Count, Pearson’s, and Spearman’s.

## Data Count.csv
This file contains the number of points shared between the two datasets being correlated. Each row is a molecule where the columns are different measurements of hemolysis. The value shown in the intersection is the number of patients for which hemolysis at that day was measured (column) as well as that molecule (row).
## Pearson’s.csv
Same as Data Count, just with Pearson’s correlation coefficients reported instead of the number of data points used for the correlation.
## Spearman’s.csv
Same as Data Count, just with Spearman’s correlation coefficients reported instead of the number of data points used for the correlation.
