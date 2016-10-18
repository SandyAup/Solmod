import os
import glob
import sys
files = glob.glob('data*.dat')

i = 0
for file_name in files:
	i += 1
	os.rename(file_name, ("data_PAMELA%i.dat" % i))
