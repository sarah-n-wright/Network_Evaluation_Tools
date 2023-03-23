import pandas as pd
import csv


if __name__ == '__main__':
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/BioPlex_293T_Network_10K_Dec_2019.tsv"
    d1 = pd.read_csv(datafile, sep="\t", header=0, engine='python', quoting=csv.QUOTE_NONE)
    print("No quoting")
    print(d1.head())
    d2 = pd.read_csv(datafile, sep="\t", header=0, quotechar='"', quoting=csv.QUOTE_ALL)
    print("QUOTE ALL")
    print(d2.head())