# -*- coding: utf-8 -*-
"""
@author: Hoel Seille (CSIRO)
"""

"""
functions that read and write EDI files

Z in linear scale; format: ohm

Pandas DF format:        
# ['FREQ',
# 'ZXXR','ZXXI','ZXX.VAR',
# 'ZXYR','ZXYI','ZXY.VAR',
# 'ZYXR','ZYXI','ZYX.VAR',
# 'ZYYR','ZYYI','ZYY.VAR',
# 'TXR.EXP','TXI.EXP','TXVAR.EXP',
# 'TYR.EXP','TYI.EXP','TYVAR.EXP']
    
"""


def read(file_name):
    
    """
    Input:
    - complete path to the EDI file
    
    Output:
    - data (EDI pandas DF format)
    - site_id (string)
    - coord (panda DF of 7 values (['Lat_deg']['Lat_min']['Lat_sec']['Long_deg']['Long_min']['Long_sec']['Elev'])    
    """   
    
    import numpy as np
    import re
    import pandas as pd

    #file_name = 'ST2,1262-27.edi'
    #file_name = 'ST_topo7-27.edi'

    coord = pd.DataFrame({'Lat_deg':0,'Lat_min':0,'Lat_sec':0,'Long_deg':0,'Long_min':0,'Long_sec':0,'Elev':0},index=[0])
    tipper = False
    #READ SITE NAME        
    with open(file_name, 'r') as f:
        data = f.readlines()
        for i,line in enumerate(data):
            line = data[i]
            words = line.split()
            
            if any("DATAID" in s for s in words):
                words = ''.join(words)
                #print words
                site_id = re.search('\"([^"]+)', words).group(1)

            #READ NUMBER OF FREQUENCIES
            if any("NFREQ" in s for s in words):
                words = ''.join(words)
                nfreq_str = (re.findall('\d+', words))
                nfreq = int(nfreq_str[0])

             #READ LATITUDE
            if any("REFLAT" in s for s in words):
                words = ''.join(words)
                reflat_str = (re.findall('\-?\d+', words))
                coord['Lat_deg'] = str(reflat_str[0])
                coord['Lat_min'] = str(reflat_str[1])
                coord['Lat_sec'] = str('.'.join(reflat_str[2:]))				

            #READ LONGITUDE
            if any("REFLONG" in s for s in words):
                words = ''.join(words)
                reflong_str = (re.findall('\-?\d+', words))
                coord['Long_deg'] = str(reflong_str[0])
                coord['Long_min'] = str(reflong_str[1])
                coord['Long_sec'] = str('.'.join(reflong_str[2:]))	

            #READ ELEVATION
            if any("REFELEV" in s for s in words):
                words = ''.join(words)
                refelev_str = (re.findall('\-?\d+', words))
                coord['Elev'] = str('.'.join(refelev_str[0:]))				

            
            #READ ELEVATION
            
            if any(">TXR.EXP" in s for s in words):
                tipper = True


    #READ MT DATA
    if tipper:
        param = ['FREQ',
                 'ZXXR','ZXXI','ZXX.VAR',
                 'ZXYR','ZXYI','ZXY.VAR',
                 'ZYXR','ZYXI','ZYX.VAR',
                 'ZYYR','ZYYI','ZYY.VAR',
                 'TXR.EXP','TXI.EXP','TXVAR.EXP',
                 'TYR.EXP','TYI.EXP','TYVAR.EXP',]
    else:
        param = ['FREQ',
                 'ZXXR','ZXXI','ZXX.VAR',
                 'ZXYR','ZXYI','ZXY.VAR',
                 'ZYXR','ZYXI','ZYX.VAR',
                 'ZYYR','ZYYI','ZYY.VAR']            
    
    edi_data = np.empty((nfreq, len(param)))
    
    with open(file_name, 'r') as f:
        data = f.readlines()
        for i,line in enumerate(data):
            line = data[i]
            words = line.split()
            
            for col, data_type in enumerate(param):
                aa=[]            
                if ('>%s' %data_type) in words:
                    for k in range (1,100):                   
                        if any(">" in s for s in data[i+k].split()):
                                break
                        else:
                            a = data[i+k].split()
                            aa += a
                    edi_data[:,col] = aa
    
    
    # write to Pandas format
    edi_pd = pd.DataFrame(edi_data)
    edi_pd.columns = param

    return (edi_pd, site_id, coord)





def write(MT_Id, coord, edi_files_path, data, write_tz=True):
    
    """
    Input:
    - MT_Id (string)
    - coord (panda DF of 7 values (['Lat_deg']['Lat_min']['Lat_sec']['Long_deg']['Long_min']['Long_sec']['Elev'])
    - path to write the file
    - data (EDI pandas DF format)
    
    Output:
    - none, writes the file only 
    
    """    
    
    

    import os
    import numpy as np
    from math import isnan
    from datetime import datetime
    import csv

    params = list(data.columns)
    
    if not write_tz:params = params[:13]
    
    datestring = datetime.strftime(datetime.now(), '%Y/%m/%d_%H:%M:%S')
    
    file_loc = os.path.join(edi_files_path, '%s.edi' % (MT_Id)) 
    file = open(file_loc,'w') 
    file.write('>HEAD\n') 
    file.write('DATAID="%s"\nACQBY=" "\nFILEBY=""\nACQDATE=""\nFILEDATE=""\n'%(MT_Id))  #ADD LEMI DALOGGER NUMBER
    file.write('PROSPECT="Area Name"\nLOC="Area Name"\n') 
    file.write('STDVERS="SEG 1.0"\nPROGVERS="DEI_MT"\nPROGDATE=%s\nMAXSECT=999\nEMPTY=1.0e+32\n'%(datestring)) 
    file.write('\n')
    file.write('>INFO\n')
    file.write('\n') 
    file.write('>=DEFINEMEAS\n')
    file.write('UNITS=M\n')
    file.write('REFLAT=%s:%s:%s\n'%(str(coord['Lat_deg'].values)[2:-2],str(coord['Lat_min'].values)[2:-2],str(coord['Lat_sec'].values)[2:-2]))
    file.write('REFLONG=%s:%s:%s\n'%(str(coord['Long_deg'].values)[2:-2],str(coord['Long_min'].values)[2:-2],str(coord['Long_sec'].values)[2:-2]))
    file.write('REFELEV=%4.2f\n'%(coord['Elev'].values))
    file.write('\n')
    file.write('>=MTSECT\n')
    file.write('SECTID="%s"\n'%MT_Id)
    file.write('NFREQ=%d\n'%len(data['FREQ'].values))
    file.write('\n') 
    file.write('>!****FREQUENCIES****!\n')
     
    for k, para in enumerate(params):
        file.write('>%s //%d\n'%(para, len(data['FREQ'].values)))
        for i in range(0, len(data['FREQ'].values)):
            if isnan(data['%s'%para].values[i]):
                file.write('%e '%1.0e+32)
            else:
                file.write('%e  '%data['%s'%para].values[i])
                if (i+1)%6==0 and i!=0: file.write('\n')
					
        file.write('\n')
        
    file.write('>END')
    
    file.close() 