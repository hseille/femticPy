# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 13:56:30 2023

@author: sei029
"""


# read file and change a column

import pandas as pd
import numpy as np

file_dir = '.'
file_name = 'output.1.ele'

ele1 = pd.read_csv('%s/%s'%(file_dir,file_name), 
                   names = ['n','a','b','c','d','e'],
                   skiprows=1,
                   skipfooter=1,
                   delim_whitespace=True,
			engine='python')

header = pd.read_csv('%s/%s'%(file_dir,file_name), nrows=1, header=None)

ele1['e'] = np.where(ele1['e'] > 30, 30, ele1['e'])

with open('%s/output.1.ele'%file_dir, 'w') as file:
    file.write('%s\n'%header[0].values[0])
    ele1.to_csv(file, header=False, index=False,sep=' ',lineterminator='\n')
    #ele1.to_csv(file, header=False, index=False,sep=' ',line_terminator='\n') # for pandas < 1.5