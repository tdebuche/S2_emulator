from collections import defaultdict
import numpy as np
import awkward as ak
import uproot
import math
import yaml

f = open('output_BeforeChange' + '.txt', 'r')
EMP = f.readlines()
f.close()

print(EMP[0])
print(EMP[1])
print(EMP[2])
print(EMP[3])
print(EMP[4])
print(EMP[5])
