# -*- coding: utf-8 -*-
"""
@author: Hoël Seillé (CSIRO)
"""



# =============================================================================
# # TODO  28/02/2024
#
# - create observed data object attributes for the different datatypes: (have independent data blocks)
        # - with the apprpriate number of sites for each one:
        #     - self.mt_ids_Z   = None
        #     - self.mt_ids_Tz  = None
        #     - self.mt_data_Z  = None
        #     - self.mt_data_Tz = None
# -  plot induction arrows
# -  plot fit Tz
# 
#  When reading back the impedance, reapply the ocnjugate
# 
# Simpllfy the result class: do not read input data, just read coordinates
# =============================================================================





import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np
import pandas as pd
import os
from copy import deepcopy

import EDItools as edi

"""
This class prepares the model and data inputs required for the inversion. 
It consists of 2 parts:

    - dataGen creates the input data files for the inversion, using .edi files   
    
    - meshGen prepares the input data required to run the sequential shell 
    script that creates the input mesh and .vtk visualization files

Files required: 

    - dataGen:
        - .edi files
        - data coordinates
        
    - meshGen:
        - data coordinates
        - coast_line 
        - topography 
        - bathymetry 
    
"""

class DataGen():

    def __init__(self, survey, outdir):
        self.survey = survey

        
        ## Size of computational domain (in km)
        self.analysis_domain = [[-500.0, 500.0],
                                [-500.0, 500.0],
                                [-500.0, 500.0]]
        
        ## Mesh center (in km)
        self.center = [0.0, 0.0, 0.0]
        
        ## Data to be inverted for
        self.invert_Z = True
        self.invert_VTF = False
        self.invert_PT = False
        self.error_floor_Z = [0.05, 0.05, 0.05, 0.05] # [Zxx, Zxy, Zyx, Zyy]
        self.error_floor_Tz = 0.05 #(absolute value)
        self.error_floor_PT = None 
        
        
        # MT sites data
        self.mt_coords_utm = None
        # self.mt_coords = None
        #self.mt_ids = []
        self.ids_Z   = []
        self.ids_VTF  = []
        self.ids_PT  = []
        self.data_Z  = []
        self.data_VTF = []
        self.data_PT = []
        #self.mt_data = None
        #self.mt_data_inversion = None
        self.nRx = 0
        self.nFreq = 0   
        self.freqs = []

        self.topography = False
        self.bathymetry = False
        self.coast_line = False
        
        self.outdir = outdir
        
    
    def read_MTdata(self, edi_path):

        site_ids=[]
        freq_list=[]
        
        for file in os.listdir(edi_path):
            if file.endswith(".edi"):
                site_ids.append(file[:-4])
        site_ids = np.sort(site_ids)
        
        for site_id in site_ids: 
                    self.nRx += 1
                    edi_file = '%s.edi'%site_id
                    edi_file_path = '%s/%s'%(edi_path, edi_file)
                    print ('    reading %s'%edi_file)
                    #read edi file
                    edi_data, site_id, coord = edi.read(edi_file_path)
                    

                    
                    if 'ZXXR' in edi_data.columns:

                        data_Z = edi_data[['FREQ',
                                 'ZXXR','ZXXI','ZXX.VAR',
                                 'ZXYR','ZXYI','ZXY.VAR',
                                 'ZYXR','ZYXI','ZYX.VAR',
                                 'ZYYR','ZYYI','ZYY.VAR']]
                        
                        # remove masked data
                        masked=np.unique(np.where(data_Z[['ZXXR','ZXXI','ZXX.VAR',
                                                      'ZXYR','ZXYI','ZXY.VAR',
                                                      'ZYXR','ZYXI','ZYX.VAR',
                                                      'ZYYR','ZYYI','ZYY.VAR']]>= 1e10)[0])
                        data_Z.drop(masked, inplace=True)
                        data_Z = data_Z.reset_index(drop=True)
                        
                        self.data_Z.append(data_Z)
                        self.ids_Z.append(site_id)
                    
                    
                    if 'TXR.EXP' in edi_data.columns:

                        data_VTF = edi_data[['FREQ',
                                 'TXR.EXP','TXI.EXP','TXVAR.EXP',
                                 'TYR.EXP','TYI.EXP','TYVAR.EXP',]]
                        
                        masked=np.unique(np.where(data_VTF[[
                                 'TXR.EXP','TXI.EXP','TXVAR.EXP',
                                 'TYR.EXP','TYI.EXP','TYVAR.EXP',]]>= 1e10)[0])
                        
                        data_VTF.drop(masked, inplace=True)
                        data_VTF = data_VTF.reset_index(drop=True)
                        
                        self.data_VTF.append(data_VTF)
                        self.ids_VTF.append(site_id)
                    
                    
                    # data.append([site_id, edi_data])
                    # site_ids.append(site_id)
                    # coords.append(coord)
                    # list all the existing frequencies in the dataset
                    for freq in range(len(edi_data['FREQ'].values)):
                        if not edi_data['FREQ'].values[freq] in freq_list:
                            if freq == 0:
                                freq_list.append(edi_data['FREQ'].values[freq])
                            else:
                                if min(np.abs(np.log10(freq_list) - np.log10(edi_data['FREQ'].values[freq]))) < 0.01:
                                    idx = (np.abs(freq_list - edi_data['FREQ'].values[freq])).argmin()
                                    edi_data['FREQ'][freq] = freq_list[idx]
                                else:
                                    freq_list.append(edi_data['FREQ'].values[freq])

        freq_list=np.sort(freq_list)[::-1]

        # self.mt_ids = site_ids
        # self.mt_data = data
        self.nRx_Z = len(self.data_Z)
        self.nRx_VTF = len(self.data_VTF)
        self.nFreq = len(freq_list)
        self.freqs = freq_list
        
        print('Read %d .edi files'%self.nRx)


    def read_MTdata_coordinates(self, coordinates_path, coordinates_file):
        
        # read coordinates_file
        coords = pd.read_csv('%s/%s'%(coordinates_path,coordinates_file), 
                             header=None, delim_whitespace=(True),
                             names = ['id','east','north','z'],
                             dtype={'id': 'string',
                                    'east': 'float64',
                                    'north': 'float64',
                                    'z': 'float64'})
        
        coords = coords.sort_values('id')
        coords = coords.reset_index(drop=True) 
        
        self.mt_coords_utm = coords
        self.mt_coords_utm['z'] = - self.mt_coords_utm['z']
        self.mt_coords = self.mt_coords_utm.copy(deep=True)
        
        

    
    def center_data(self):
        """
        This function defines the coordinates of the MT stations relative to
        the mesh center
        """    
        east_min = self.mt_coords_utm['east'].min()
        east_max = self.mt_coords_utm['east'].max()
        north_min = self.mt_coords_utm['north'].min()
        north_max = self.mt_coords_utm['north'].max()
        east_center = east_min + (east_max - east_min) / 2
        north_center = north_min + (north_max - north_min) / 2
        
        self.mt_coords['east'] -= east_center
        self.mt_coords['north'] -= north_center
        self.anchor = [north_center, east_center]



    def convert_units(self, data_orig):
        # Conversion to appropriate units: FEMTIC uses ohms
        #   1 ohm = 10000(4*pi) [mV/km/nT]
        C = 10000/(4*np.pi)
        data = data_orig.copy(deep=True)
        data.loc[:,['ZXXR','ZXXI','ZXYR','ZXYI','ZYXR','ZYXI','ZYYR','ZYYI']] = data_orig.loc[:,['ZXXR','ZXXI','ZXYR','ZXYI','ZYXR','ZYXI','ZYYR','ZYYI']] / C
        data.loc[:,['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']] = data_orig.loc[:,['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']] / (C**2)
        return data 
    
    
    # def apply_error_floor_Z(self, data):
    #     #calculate error floor
    #     std = data[['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']]**0.5  
    #     for freq in range(len(data)):
    #         error_floor = ((((data['ZXYR'][freq]**2+data['ZXYI'][freq]**2)**0.5) * 
    #                         ((data['ZYXR'][freq]**2+data['ZYXI'][freq]**2) ** 0.5))**0.5)  * np.array(self.error_floor_Z)
            
    #         for i, comp in enumerate(['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']):
    #             if std[comp][freq] < error_floor[i]:
    #                 data[comp][freq] =  (error_floor[i])**2

    #     return data
    
    
    # def apply_error_floor_Tz(self, data):
    # #calculate error floor
    # std = Tz[['TXVAR.EXP','TYVAR.EXP']]**0.5   
    # for freq in range(len(data)):
    #     # error_floor = ((((Z['ZXYR'][freq]**2+Z['ZXYI'][freq]**2)**0.5) * 
    #     #                 ((Z['ZYXR'][freq]**2+Z['ZYXI'][freq]**2) ** 0.5))**0.5)  * np.array(self.error_floor_Z)
        
    #     for comp in enumerate(['TXVAR.EXP','TYVAR.EXP']):
    #         if std[comp][freq] < self.error_floor_Tz:
    #             data[comp][freq] =  (self.error_floor_Tz)**2

    # return data
  

    def apply_error_floor(self, data, data_type = 'Z'):
        
        if data_type == 'Z':
            #calculate error floor
            std = data[['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']]**0.5  
            for freq in range(len(data)):
                
                error_floor = ((((data['ZXYR'][freq]**2+data['ZXYI'][freq]**2)**0.5) * 
                                ((data['ZYXR'][freq]**2+data['ZYXI'][freq]**2) ** 0.5))**0.5)  * np.array(self.error_floor_Z)
                
                for i, comp in enumerate(['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']):
                    if std[comp][freq] < error_floor[i]:
                        data[comp][freq] =  (error_floor[i])**2
            
        elif data_type == 'VTF':
            std = data[['TXVAR.EXP','TYVAR.EXP']]**0.5   
            for freq in range(len(data)):
                for i, comp in enumerate(['TXVAR.EXP','TYVAR.EXP']):
                    if std[comp][freq] < self.error_floor_Tz:
                        data[comp][freq] =  (self.error_floor_Tz)**2

        return data
  
    
    def conjugate(self, data, data_type = 'Z'):
        # FEMTIC convention takes the complex conjugate of Z (exp(+iwt) assumed)
        
        if data_type == 'Z':
            data['ZXXI'] = -data['ZXXI']
            data['ZXYI'] = -data['ZXYI']
            data['ZYXI'] = -data['ZYXI']
            data['ZYYI'] = -data['ZYYI']
            
        elif data_type == 'VTF':
            data['TXI.EXP'] = -data['TXI.EXP']
            data['TYI.EXP'] = -data['TYI.EXP']
            
        return data
    
    
            
    def prep_data(self, freq_bandwidth = None, subsampling = 1):
        
        # subsample the original frequency list:
        if isinstance(subsampling, int):
            # subsample every th frequency
            self.freqs_inversion = self.freqs[::subsampling]
        else:
            # subsambple at the indexes specified in the list provided as subsampling
            self.freqs_inversion = self.freqs[subsampling]
        self.freqs_inversion = self.freqs_inversion[np.where((self.freqs_inversion<=freq_bandwidth[1]) & (self.freqs_inversion>=freq_bandwidth[0]))]
        
        #self.mt_data_inversion = deepcopy(self.mt_data)
        # flag for specifying that data has already been prepared ?
        
        if self.invert_Z:
        
            for rx in range(self.nRx_Z):
                
                rx_data = self.data_Z[rx]
                rx_data = rx_data[rx_data['FREQ'].isin(self.freqs_inversion)]
                #rx_data = rx_data.iloc[::subsampling, :]
                rx_data.reset_index(inplace=True, drop=True)
                
                # we convert the units 
                rx_data_ = self.convert_units(rx_data)
                # we take the complex conjugate 
                rx_data_ = self.conjugate(rx_data_,data_type = 'Z')   
                # we apply the error floors
                rx_data_ = self.apply_error_floor(rx_data_,data_type = 'Z')
                # if self.invert_Z:
                #     rx_data_ = self.apply_error_floor_Z(rx_data_)
                # if self.invert_VTF:
                #     rx_data_ = self.apply_error_floor_Tz(rx_data_)                
    
                self.data_Z[rx] = rx_data_


        if self.invert_VTF:
        
            for rx in range(self.nRx_VTF):
                                
                rx_data = self.data_VTF[rx]
                rx_data = rx_data[rx_data['FREQ'].isin(self.freqs_inversion)]
                #rx_data = rx_data.iloc[::subsampling, :]
                rx_data.reset_index(inplace=True, drop=True)

                # we take the complex conjugate
                rx_data_ = self.conjugate(rx_data,data_type = 'VTF')   
                # we apply the error floors
                rx_data_ = self.apply_error_floor(rx_data_,data_type = 'VTF')
                # if self.invert_Z:
                #     rx_data_ = self.apply_error_floor_Z(rx_data_)
                # if self.invert_VTF:
                #     rx_data_ = self.apply_error_floor_Tz(rx_data_)                
    
                self.data_VTF[rx] = rx_data_

        
    
    def write_observe(self, write = True, freq_bandwidth = None, subsampling = 1, plot=True):
        
        self.prep_data(freq_bandwidth = freq_bandwidth, subsampling = subsampling)
        self.plot_inversion_periods()
        
        file_loc = os.path.join(self.outdir, 'observe.dat') 
        file = open(file_loc,'w')

        if self.invert_Z:
        
            file.write('MT    %d\n'%(self.nRx_Z))
            
            for rx in range(self.nRx_Z):
    
                ind = np.where(self.ids_Z[rx] == self.mt_coords['id'].values)[0][0]
                
                data = self.data_Z[rx]
                
                file.write('%d  %d  %2.3f  %2.3f  \n'%(rx+1, rx+1, 
                                                     self.mt_coords['north'][ind],
                                                     self.mt_coords['east'][ind]))                
                nFreq = len(data)
                file.write('%d\n'%nFreq)
                
                for freq in range(nFreq):
                    file.write('%.5f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e \n'%(
                               data['FREQ'][freq],
                               data['ZXXR'][freq],data['ZXXI'][freq],data['ZXYR'][freq],data['ZXYI'][freq],
                               data['ZYXR'][freq],data['ZYXI'][freq],data['ZYYR'][freq],data['ZYYI'][freq],
                               data['ZXX.VAR'][freq]**0.5,data['ZXX.VAR'][freq]**0.5,data['ZXY.VAR'][freq]**0.5,data['ZXY.VAR'][freq]**0.5,
                               data['ZYX.VAR'][freq]**0.5,data['ZYX.VAR'][freq]**0.5,data['ZYY.VAR'][freq]**0.5,data['ZYY.VAR'][freq]**0.5))
        
        
        if self.invert_VTF:
        
            file.write('VTF    %d\n'%(self.nRx_VTF))
            
            for rx in range(self.nRx_VTF):
    
                ind = np.where(self.ids_VTF[rx] == self.mt_coords['id'].values)[0][0]
                
                file.write('%d  %d  %2.3f  %2.3f  \n'%(rx+1001, rx+1001, 
                                                     self.mt_coords['north'][ind],
                                                     self.mt_coords['east'][ind]))
                data_tz = self.data_VTF[rx]
                
                # remove masked tipper points (>1+8) or inexistant points (<1e-8)
                
                for freq in range(len(data_tz)):
                    if (data_tz['TXR.EXP'][freq] > 1e+8 or
                        data_tz['TXR.EXP'][freq] > 1e+8 or 
                        data_tz['TXI.EXP'][freq] > 1e+8 or 
                        data_tz['TXI.EXP'][freq] > 1e+8 or
                        data_tz['TXR.EXP'][freq] < 1e+8 or
                        data_tz['TXR.EXP'][freq] < 1e+8 or 
                        data_tz['TXI.EXP'][freq] < 1e+8 or 
                        data_tz['TXI.EXP'][freq] < 1e+8) :
                            data_tz.drop(data_tz.index[freq])
                
                data_tz = data_tz.reset_index(drop=True)  
                
                nFreq = len(data_tz)
                file.write('%d\n'%nFreq)
                
                for freq in range(nFreq):
                    file.write('%.5f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e \n'%(
                               data_tz['FREQ'][freq],
                               data_tz['TXR.EXP'][freq],data_tz['TXI.EXP'][freq],data_tz['TYR.EXP'][freq],data_tz['TYI.EXP'][freq],
                               data_tz['TXVAR.EXP'][freq]**0.5,data_tz['TXVAR.EXP'][freq]**0.5,
                               data_tz['TYVAR.EXP'][freq]**0.5,data_tz['TYVAR.EXP'][freq]**0.5))
               
        
        
        file.write('END')
        file.close()
        

        
        print('Inversion data summary: \n')  
        print(f'  Inverting for {len(self.freqs_inversion)} frequencies\n\n')
        if self.invert_Z:        
            data_points_Z = np.hstack([self.data_Z[i]['FREQ'].values for i in range(self.nRx_Z)])
            unique, counts = np.unique(data_points_Z, return_counts=True)
            print('  Impedance Z: \n')
            print(f'  {self.nRx_Z} MT sites\n')
            print(f'  {len(data_points_Z)*8} Impedance tensor data points (full Z tensor, Real and Imag data) \n')
            print('  Data points per frequencies: \n')
            for i in range(len(unique)):
                #print(f'    T = {1/unique[i]} s  nData = {counts[i]} \n')
                print('    T = %.3Es  nData = %3d \n'%(1/unique[i],counts[i]))


        if self.invert_VTF:
            data_points_VTF = np.hstack([self.data_VTF[i]['FREQ'].values for i in range(self.nRx_VTF)])
            unique, counts = np.unique(data_points_VTF, return_counts=True)
            print('  Tipper VTF: \n')
            print(f'  {self.nRx_VTF} MT sites\n')
            print(f'  {len(data_points_VTF)*4} VTF full tensor (Real and Imag data) \n')
            print('  Data points per frequencies: \n')
            for i in range(len(unique)):
                #print(f'    T = {1/unique[i]} s  nData = {counts[i]} \n')
                print('    T = %.3Es  nData = %3d \n'%(1/unique[i],counts[i]))
            
                
    
            
    def plot_inversion_periods(self):
        
        fig, ax = plt.subplots(1,1, figsize=(10,3))
        
        if self.invert_Z:        
            data_points_Z = np.hstack([self.data_Z[i]['FREQ'].values for i in range(self.nRx_Z)])
            unique, counts = np.unique(data_points_Z, return_counts=True)
            ax.bar(np.log10(1/unique), counts, 0.1, color='k', alpha=0.5, label='Z') 

        if self.invert_VTF:
            data_points_VTF = np.hstack([self.data_VTF[i]['FREQ'].values for i in range(self.nRx_VTF)])
            unique, counts = np.unique(data_points_VTF, return_counts=True)
            ax.bar(np.log10(1/unique), counts, 0.1, color='r', alpha=0.5, label='VTF') 

                
        
        for i in range(self.nFreq):
            ax.axvline(x=np.log10(1/self.freqs[i]))
        for i in range(len(self.freqs_inversion)):
            ax.axvline(x=np.log10(1/self.freqs_inversion[i]), c= 'r', ls='--')   

        # ax.set_xscale('log')
        ax.set_title(f'{len(self.freqs_inversion)} inversion periods, {len(self.freqs)} original periods')
        
        ax.axvline(x=np.log10(1/self.freqs[0]), label='original periods')
        ax.axvline(x=np.log10(1/self.freqs_inversion[0]), c= 'r', ls='--', label='inversion periods')
        ax.legend()
        ax.set_xlabel('Log10 Periods (s)')
        ax.set_ylabel('nData')
            
            
    
    def write_inversion_control(self,
                                INV_METHOD = 1,
                                NUM_THREADS = 1,
                                DISTORTION = 0,
                                TRADE_OFF_PARAM = [3],
                                ITERATION = 10,
                                CONVERGE = 1.0,
                                ALPHA_WEIGHT = [1,1,1]
                                ):
        file_loc = os.path.join(self.outdir, 'control.dat') 
        file = open(file_loc,'w')
        file.write('NUM_THREADS\n')
        file.write('%d\n'%NUM_THREADS)
        file.write('MESH_TYPE\n')
        file.write('1\n')
        file.write('OUTPUT_PARAM_VTK\n')
        file.write('2\n')
        file.write('0 4\n')
        file.write('OFILE_TYPE\n')
        file.write('0\n')
        file.write('DISTORTION\n')
        file.write('%d\n'%DISTORTION)
        file.write('ALPHA_WEIGHT\n')
        for i in range(len(ALPHA_WEIGHT)):
            file.write('%.2f '%ALPHA_WEIGHT[i]) 
        file.write('\n')
        file.write('INV_METHOD\n')
        file.write('%d\n'%INV_METHOD)
        file.write('BOTTOM_RESISTIVITY\n')
        file.write('100.0\n')
        file.write('BOTTOM_ROUGHNING_FACTOR\n')
        file.write('1.0\n')
        file.write('OUTPUT_OPTION\n')
        file.write('0 0\n')
        file.write('TRADE_OFF_PARAM\n')
        for i in range(len(TRADE_OFF_PARAM)):
            file.write('%.2f '%TRADE_OFF_PARAM[i])      
        file.write('\n')    
        file.write('ITERATION\n')
        file.write('0 %d\n'%ITERATION)        
        file.write('RETRIAL\n')
        file.write('2\n') 
        file.write('CONVERGE\n')       
        file.write('%.1f\n'%CONVERGE)  
        file.write('STEP_LENGTH\n')
        file.write('0.5 0.1 1.0\n')           
        file.write('3\n')
        file.write('0.5 1.25\n')           
        file.write('END')
        file.close()    
                

    
    def plot_data_loc(self, plot_ids = False, zoom_core=False):
        fig, ax = plt.subplots(1,1,figsize=(8, 8))
        ax.plot(self.mt_coords['east'], self.mt_coords['north'], 'kv', ms=5, label = 'MT')
        if plot_ids:
            for mt_id in range(self.nRx_Z):
                ax.text(self.mt_coords['east'][mt_id], self.mt_coords['north'][mt_id]+0.1,
                         self.mt_coords['id'][mt_id])
            
        domain_boudary = np.array([[self.analysis_domain[0][0],
                          self.analysis_domain[0][0],
                          self.analysis_domain[0][1],
                          self.analysis_domain[0][1],
                          self.analysis_domain[0][0]],
                          [self.analysis_domain[1][0],
                           self.analysis_domain[1][1],
                           self.analysis_domain[1][1],
                           self.analysis_domain[1][0],
                           self.analysis_domain[1][0]]])
        
        ax.plot(domain_boudary[0], domain_boudary[1]  , 'k-', 
                 label = 'analysis_domain')
        
        if zoom_core:
            ax.set_xlim(self.mt_coords['east'].min()-10, self.mt_coords['east'].max()+10)
            ax.set_ylim(self.mt_coords['north'].min()-10, self.mt_coords['north'].max()+10)
        
        ax.legend()




    def plot_coast_line(self):
            nPolys = pd.read_csv(self.coast_line, nrows=1, names=['n'])['n'][0]
            coast_line_pd = pd.read_csv(self.coast_line, skiprows=1, delim_whitespace= True, names = ['north','east','i','j'])
            cb = np.where((coast_line_pd['i']==1) | (coast_line_pd['i']==-1))[0]
            
            for i in range(nPolys):
                if i == 0:
                    plt.plot(coast_line_pd['east'][:cb[i]+1],coast_line_pd['north'][:cb[i]+1], 'k-',lw=1)
                else:
                    plt.plot(coast_line_pd['east'][cb[i-1]+1:cb[i]+1],coast_line_pd['north'][cb[i-1]+1:cb[i]+1], 'k-',lw=1)
            

    def plot_topo_bathy(self):

            topo_pd = pd.read_csv(self.topography, delim_whitespace= True, names = ['north','east','z'])
            topo_pd = topo_pd[topo_pd['z'] > 0]
            plt.scatter(topo_pd['east'],topo_pd['north'],c=topo_pd['z'],cmap = 'Greens', s=2)  
  
            bathy_pd = pd.read_csv(self.bathymetry, delim_whitespace= True, names = ['north','east','z'])
            bathy_pd = bathy_pd[bathy_pd['z'] > 0]
            plt.scatter(bathy_pd['east'],bathy_pd['north'],c=bathy_pd['z'] ,cmap = 'Blues', s=2)                   
        

    #def plot_input_data(self):
        



class MeshGen():
    """
    MeshGen prepares the input data required to run the sequential shell 
    script that creates the input model and .vtk visualization files
    
    Files required: 
        - data coordinates
        - coast_line 
        - topography 
        - bathymetry 
    
    Writes out:
        - analysis_domain.dat
        - control.dat
        - makeMtr.param
        - obs_site.dat
        - observing_site.dat
        - resistivity_attr.dat

    """
    
    def __init__(self, survey, analysis_domain, center, mt_coords, nRx, outdir):
        self.survey = survey
        
        self.outdir = outdir
                
        ## Domains to include
        self.land = True
        self.sea = False
        #self.coast_line = False
        
        
        ## Size of computational domain (in km)
        self.analysis_domain = analysis_domain
        
        ## Mesh 
        self.center = center  #center (in km)
        self.rotation = 0.0
        
        # MT sites data
        self.mt_coords = mt_coords
        self.median_elevation = - np.median(self.mt_coords['z'])
        self.center = [self.center[0], self.center[1], self.median_elevation]
        # self.mt_data = mt_data
        self.nRx = nRx

        
        # Information about the ellipsoids to  control edge lengths
        self.ellipsoids_control = [6,
                            [40.0,  1.0, 0.0, 0.5, 0.7],
                            [60.0,  5.0, 0.0, 0.3, 0.5],
                            [100.0, 10.0, 0.0, 0.1, 0.3],
                            [200.0, 20.0, 0.0, 0.0, 0.0],
                            [300.0, 30.0, 0.0, 0.0, 0.0],
                            [500.0, 50.0, 0.0, 0.0, 0.0]]

        
        # topography interpolation parameters
        self.interpolate_sr = 100
        self.interpolate_npts = 3
        self.interpolate_minval = 1e-6
        # altitude parameters:
        self.altitude_file = 'topography.dat'
        self.altitude_min = 0.0
        self.altitude_max = 1e20
        # sea_depth parameters:
        self.sea_depth_file = 'bathymetry.dat'
        self.sea_depth_min = 0.01
        self.sea_depth_max = 1e20     
        
        self.ellipsoids_mtr = [10,
                            [40.0,    1.0,  0.0,  0.7,  0.9],
                            [45.0,   1.5,  0.0,  0.5,  0.7],
                            [50.0,    3.0,  0.0,  0.4,  0.7],
                            [60.0,    5.0,  0.0,  0.3,  0.5],
                            [80.0,    8.0,  0.0,  0.1,  0.3],
                            [100.0,  10.0,  0.0,  0.0,  0.0],
                            [200.0,  20.0,  0.0,  0.0,  0.0],
                            [300.0,  30.0,  0.0,  0.0,  0.0],
                            [400.0,  40.0,  0.0,  0.0,  0.0],
                            [500.0,  50.0,  0.0,  0.0,  0.0]]
        
        self.ellipsoids_observing_sites =  [5, 
                                            0.1, 0.02,  
                                            0.3, 0.05,  
                                            1.0, 0.10,  
                                            3.0, 0.30,  
                                            5.0, 0.50]        

        self.ellipsoids_obs_sites = [6,
                            [0.5, 0.10, 0.3],
                            [1.0, 0.20, 0.3],
                            [1.5, 0.30, 0.3],
                            [2.0, 0.50, 0.3],
                            [3.0, 1.00, 0.3],
                            [5.0, 2.00, 0.3]]

        
        self.ellipsoids_resistivity_attr = [9,
                            [40.0,       2.0,  0.0,  0.7],
                            [45.0  ,     3.0,  0.0,  0.7],
                            [50.0,       5.0,  0.0,  0.7],
                            [60.0,      10.0,  0.0,  0.6],
                            [100.0,    100.0,  0.0,  0.5],
                            [200.0,    200.0,  0.0,  0.3],
                            [300.0,    300.0,  0.0,  0.2],
                            [500.0,    500.0,  0.0,  0.1],
                            [1000.0,  1000.0,  0.0,  0.0]] 
        
        
        self.ellipsoids_resistivity_attr_sites = [2,
                            [3.0, 2.0],
                            [5.0, 3.0]]
        
        self.resistivity_sea = 0.25
        self.resistivity_starting_model = 100

    def write_analysis_domain(self):
        file_loc = os.path.join(self.outdir, 'analysis_domain.dat') 
        file = open(file_loc,'w') 
        for i in range(3):
            file.write(' %4.1f %4.1f\n'%(self.analysis_domain[i][0], self.analysis_domain[i][1]))
        file.close() 
        
    
    def write_control(self):
        file_loc = os.path.join(self.outdir, 'control.dat') 
        file = open(file_loc,'w')      
        file.write('CENTER\n')
        file.write('%4.1f %4.1f %4.1f\n'%(self.center[0],self.center[1],self.center[2]))
        file.write('ROTATION\n')
        file.write('%4.1f\n'%(self.rotation))
        file.write('NUM_THREADS\n1\nSURF_MESH\nELLIPSOIDS\n')
        file.write('%d'%(self.ellipsoids_control[0]))
        for tt in self.ellipsoids_control[1:]:
            file.write('\n')
            for ttt in tt:
                file.write('%4.1f '%(ttt))
        file.write('\n')
        file.write('INTERPOLATE\n')
        file.write('%4.1f\n'%self.interpolate_sr)
        file.write('%4.1f\n'%self.interpolate_npts)
        file.write('%4.7f\n'%self.interpolate_minval)

        file.write('ALTITUDE\n')
        file.write('topography.dat\n')
        file.write('%4.1f\n'%self.altitude_min)
        file.write('%4.1f\n'%self.altitude_max)          

        file.write('SEA_DEPTH\n')
        file.write('bathymetry.dat\n')
        file.write('%4.1f\n'%self.sea_depth_min)
        file.write('%4.1f\n'%self.sea_depth_max) 
        file.write('END')
        file.close() 


    def write_makeMtr(self):
        file_loc = os.path.join(self.outdir, 'makeMtr.param') 
        file = open(file_loc,'w')      
        file.write('%4.1f %4.1f %4.1f\n'%(self.center[0],self.center[1],self.center[2]))
        file.write('%4.1f\n'%(self.rotation))
        file.write('%d'%(self.ellipsoids_mtr[0]))
        for tt in self.ellipsoids_mtr[1:]:
            file.write('\n')
            for ttt in tt:
                file.write('%4.1f '%(ttt))
        file.close() 
        

    def write_obs_site(self):
        file_loc = os.path.join(self.outdir, 'obs_site.dat') 
        file = open(file_loc,'w')   
        file.write('%d\n'%self.nRx)
        for rx in range(self.nRx):
            file.write('%2.3f %2.3f %2.3f\n'%(self.mt_coords['north'][rx],
                                              self.mt_coords['east'][rx],
                                              self.mt_coords['z'][rx]))
            file.write('%d'%(self.ellipsoids_obs_sites[0]))
            for tt in self.ellipsoids_obs_sites[1:]:
                file.write('\n')
                for ttt in tt:
                    file.write('%4.3f '%(ttt))
            file.write('\n')
        file.close() 



    def write_observing_site(self):
        file_loc = os.path.join(self.outdir, 'observing_site.dat') 
        file = open(file_loc,'w')    
        file.write('%d\n'%self.nRx)
        for rx in range(self.nRx):
            file.write('%2.3f  %2.3f  '%(self.mt_coords['north'][rx],
                                         self.mt_coords['east'][rx]))
            file.write('%d  '%(self.ellipsoids_observing_sites[0]))
            for tt in self.ellipsoids_observing_sites[1:]:
                file.write('%4.3f '%(tt))
            file.write('\n')
        file.close() 
        
        
    def write_dummy_coast_line(self):
        file_loc = os.path.join(self.outdir, 'coast_line.dat') 
        file = open(file_loc,'w')   
        # define a boundary which covers whole of the computational domain.
        domain_boudary = np.array([[self.analysis_domain[0][0]-1,
                          self.analysis_domain[0][0]-1,
                          self.analysis_domain[0][1]+1,
                          self.analysis_domain[0][1]+1,
                          self.analysis_domain[0][0]-1],
                          [self.analysis_domain[1][0]-1,
                            self.analysis_domain[1][1]+1,
                            self.analysis_domain[1][1]+1,
                            self.analysis_domain[1][0]-1,
                            self.analysis_domain[1][0]-1]])
        
        file.write('1\n')
        for corner in range(3):    
            file.write('%4.3f %4.3f 0 0\n'%(domain_boudary[0][corner],
                                            domain_boudary[1][corner]))
        file.write('%4.3f %4.3f 1 0\n'%(domain_boudary[0][3],
                                        domain_boudary[1][3]))
        file.close() 

        
        
    def write_dummy_topo_file(self):
        file_loc = os.path.join(self.outdir, 'topography.dat') 
        file = open(file_loc,'w') 
        for north in np.arange(self.analysis_domain[0][0]-11, self.analysis_domain[0][1]+11,10):
            for east in np.arange(self.analysis_domain[1][0]-11, self.analysis_domain[1][1]+11,10):
                file.write('%4.3f %4.3f %4.3f\n'%(north, east, 0.000))
        file.close() 



    def write_dummy_bathy_file(self):
        file_loc = os.path.join(self.outdir, 'bathymetry.dat') 
        file = open(file_loc,'w') 
        for north in np.arange(self.analysis_domain[0][0]-11, self.analysis_domain[0][1]+11,10):
            for east in np.arange(self.analysis_domain[1][0]-11, self.analysis_domain[1][1]+11,10):
                file.write('%4.3f %4.3f %4.3f\n'%(north, east, -0.01))
        file.close() 

        

    def write_resistivity_attr(self):
        file_loc = os.path.join(self.outdir, 'resistivity_attr.dat') 
        file = open(file_loc,'w') 
        self.region_attributes = np.sum([self.land, self.sea]) + 1
        file.write('%d\n'%self.region_attributes)
        file.write('10 1.0e+9 -1 1\n')
        if self.sea:
            file.write('20 %.3f  -1 1\n'%self.resistivity_sea)
            file.write('30 %.1f   9 0\n'%self.resistivity_starting_model)
        else:
            file.write('20 %.1f   9 0\n'%self.resistivity_starting_model)
            
        file.write('%4.1f %4.1f %4.1f\n'%(self.center[0],self.center[1],self.center[2]))
        file.write('%4.1f\n'%(self.rotation))
        
        file.write('%d'%(self.ellipsoids_resistivity_attr[0]))
        for tt in self.ellipsoids_resistivity_attr[1:]:
            file.write('\n')
            for ttt in tt:
                file.write('%4.1f '%(ttt))
        file.write('\n')
        file.write('%d\n'%self.nRx)
        for rx in range(self.nRx):
            file.write('%2.3f  %2.3f  %2.3f\n'%(self.mt_coords['north'][rx],
                                         self.mt_coords['east'][rx],
                                         self.mt_coords['z'][rx]))
            file.write('%d  '%(self.ellipsoids_resistivity_attr_sites[0]))
            for tt in self.ellipsoids_resistivity_attr_sites[1:]:
                file.write('\n')
                for ttt in tt:
                    file.write('%4.3f '%(ttt))
            file.write('\n')
        file.close() 
        
    
    def write_inputs(self):
        self.write_analysis_domain()
        self.write_makeMtr()
        self.write_obs_site()
        self.write_observing_site()
        self.write_resistivity_attr()
        self.write_control()
        





     
class InvResults():
    """
    InvResults plot the statistics of the inversion
    No Files required
    """
    
    # def __init__(self, survey, mt_coords, mt_data, mt_ids, results_directory):
    def __init__(self, survey, mt_coords, ids_Z, ids_VTF, results_directory):
        self.mt_coords = mt_coords
        self.ids_Z = ids_Z
        self.ids_VTF = ids_VTF
        #self.nRx = len(mt_data)
        self.dir = results_directory 
        self.plot_dir = None
 
        self.nIter = None
        
        self.rms_total = None
        self.rms_Z = None
        self.rms_VTF = None
        self.rms_breakDown_Z = None
        self.rms_breakDown_VTF = None
        self.log_inversion = None
        
        ## Data to be inverted for
        self.invert_Z = False
        self.invert_VTF = False
        self.invert_PT = False
        
        #self.obs = mt_data
        self.resp_Z = None  
        self.resp_VTF = None 
        self.resp_PT = None 
    
    
    
    def create_res_dir(self):
        if self.plot_dir is None:
            self.plot_dir = './%s/plots/it%02d/'%(self.dir,self.nIter)
            if not os.path.isdir('./%s/plots/'%self.dir):
                os.mkdir('./%s/plots'%self.dir)
            if not os.path.isdir(self.plot_dir):
                os.mkdir(self.plot_dir)  
                

    def open_single_csv(self, csv_path):

        nMT  = 0
        nVTF = 0
        nPT  = 0
    
        with open(csv_path, 'r') as file:
            
            lines = file.readlines()
            
            for i, line in enumerate(lines):
                line = line.strip()
                if 'MT' in line:
                    data_type = 'MT'
                    line_MT = i
                    self.invert_Z = True            
                elif 'VTF' in line:
                    data_type = 'VTF'
                    line_VTF = i
                    self.invert_VTF = True
                elif 'PT' in line:
                    data_type = 'PT'
                    line_PT = i
                    self.invert_PT = True
                elif 'StaID' in line:
                    pass
                elif data_type:  
                    if data_type == 'MT':
                        nMT += 1
                    elif data_type == 'VTF':
                        nVTF += 1
                    elif data_type == 'PT':
                        nPT += 1
    
        #nDatasets = invert_Z + invert_VTF + invert_PT
        
        class output():
            pass
    
        if self.invert_Z:
            headers = [header.strip() for header in lines[line_MT+1].split(",")]
            data = [line.split(",") for line in lines[line_MT+2:line_MT+2+nMT]]
            df_Z = pd.DataFrame(data, columns=headers)
            df_Z.pop(df_Z.columns[-1])
            # cols = df_Z.columns.difference(['StaID'])
            # df_Z[cols] = df_Z[cols].astype(float).astype(int)
            # df_Z['StaID'] = df_Z['StaID'].apply(lambda x: x.strip())
            df_Z['StaID'] = df_Z['StaID'].astype(int) 
            df_Z = df_Z.astype({col: float for col in df_Z.columns[1:]})
            output.Z = df_Z
    
        if self.invert_VTF:
            headers = [header.strip() for header in lines[line_VTF+1].split(",")]
            data = [line.split(",") for line in lines[line_VTF+2:line_VTF+2+nVTF]]
            df_VTF = pd.DataFrame(data, columns=headers)
            df_VTF.pop(df_VTF.columns[-1])
            # cols = df_VTF.columns.difference(['StaID'])
            # df_VTF[cols] = df_VTF[cols].astype(float).astype(int)
            # df_VTF['StaID'] = df_VTF['StaID'].apply(lambda x: x.strip())
            df_VTF['StaID'] = df_VTF['StaID'].astype(int) 
            df_VTF = df_VTF.astype({col: float for col in df_VTF.columns[1:]})
            output.VTF = df_VTF
    
        if self.invert_PT:
            headers = [header.strip() for header in lines[line_VTF+1].split(",")]
            data = [line.split(",") for line in lines[line_PT+2:line_PT+2+nPT]]
            df_PT = pd.DataFrame(data, columns=headers)
            df_PT.pop(df_PT.columns[-1])
            # cols = df_PT.columns.difference(['StaID'])
            # df_PT[cols] = df_PT[cols].astype(float).astype(str)
            # df_PT['StaID'] = df_PT['StaID'].apply(lambda x: x.strip())
            df_PT['StaID'] = df_PT['StaID'].astype(int) 
            df_PT.astype({col: float for col in df_PT.columns[1:]})
            output.PT = df_PT
    
        return output
    
    
    def read_result_csv(self):
        
        self.create_res_dir()
    
        for root, dirs, files in os.walk(self.dir):
            for file in (files):
                if file.endswith('%d.csv'%self.nIter):
                    csv_path = '%s/%s'%(self.dir, file)
                    data_csv = self.open_single_csv(csv_path)
                    
                    if self.invert_Z:
                        try:
                            resp_Z_all
                        except NameError:
                            resp_Z_all = data_csv.Z
                        else:
                            resp_Z_all = pd.concat([resp_Z_all, data_csv.Z], ignore_index=True)
        
                            
                    if self.invert_VTF:
                        try:
                            resp_VTF_all
                        except NameError:
                            resp_VTF_all = data_csv.VTF
                        else:
                            resp_VTF_all = pd.concat([resp_VTF_all, data_csv.VTF], ignore_index=True)
                            
                            
                    if self.invert_PT:
                        try:
                            resp_PT_all
                        except NameError:
                            resp_PT_all = data_csv.PT
                        else:
                            resp_PT_all = pd.concat([resp_PT_all, data_csv.PT], ignore_index=True)
        
        
        if self.invert_Z:
            resp_Z_all.columns = resp_Z_all.columns.str.lstrip()
            resp_Z_all = resp_Z_all.sort_values(['StaID', 'Freq[Hz]'],
                          ascending = [True, False])     
        
            resp_Z_all = resp_Z_all.reset_index(drop=True)   
        
            for site in range(len(self.ids_Z)):
                idx = np.where(resp_Z_all['StaID'] == site+1) 
                # print(idx[0], '----', self.mt_ids[site])
                resp_Z_all['StaID'][idx[0]] = self.ids_Z[site]
        
            self.resp_Z = resp_Z_all        
    

        if self.invert_VTF:
            resp_VTF_all.columns = resp_VTF_all.columns.str.lstrip()
            resp_VTF_all = resp_VTF_all.sort_values(['StaID', 'Freq[Hz]'],
                          ascending = [True, False])     
        
            resp_VTF_all = resp_VTF_all.reset_index(drop=True)   
        
            for site in range(len(self.ids_VTF)):
                idx = np.where(resp_VTF_all['StaID'] == site+1001) 
                resp_VTF_all['StaID'][idx[0]] = self.ids_VTF[site]
                
            resp_VTF_all['Im(Tzx)_Cal'] = -resp_VTF_all['Im(Tzx)_Cal']
            resp_VTF_all['Im(Tzy)_Cal'] = -resp_VTF_all['Im(Tzy)_Cal']
            resp_VTF_all['Im(Tzx)_Obs'] = -resp_VTF_all['Im(Tzx)_Obs']
            resp_VTF_all['Im(Tzy)_Obs'] = -resp_VTF_all['Im(Tzy)_Obs']
        
            self.resp_VTF = resp_VTF_all    




        if self.invert_PT:
            resp_PT_all.columns = resp_PT_all.columns.str.lstrip()
            resp_PT_all = resp_PT_all.sort_values(['StaID', 'Freq[Hz]'],
                          ascending = [True, False])     
        
            resp_PT_all = resp_PT_all.reset_index(drop=True)   
        
            for site in range(len(self.mt_ids)):
                idx = np.where(resp_PT_all['StaID'] == site+1) 
                resp_PT_all['StaID'][idx[0]] = self.mt_ids[site]
        
            self.resp_PT = resp_PT_all
    
    
    
    # def read_result_csv(self):

    #     for root, dirs, files in os.walk(self.dir):
    #         for file in (files):
    #             if file.endswith('%d.csv'%self.nIter):
    #                 csv_path = '%s/%s'%(self.dir, file)
    #                 resp = pd.read_csv(csv_path, skiprows=1)
    #                 try:
    #                     resp_all
    #                 except NameError:
    #                     resp_all = resp
    #                 else:
    #                     #resp_all = resp_all.append(resp, ignore_index=True)
    #                     resp_all = pd.concat([resp_all, resp], ignore_index=True)
                        
    #     resp_all.columns = resp_all.columns.str.lstrip()
    #     resp_all = resp_all.sort_values(['StaID', 'Freq[Hz]'],
    #                   ascending = [True, False])     
    
    #     resp_all = resp_all.reset_index(drop=True)   
        
    #     for site in range(len(self.mt_ids)):
    #         idx = np.where(resp_all['StaID'] == site+1) 
    #         resp_all['StaID'][idx[0]] = self.mt_ids[site]
        
    #     self.resp = resp_all


    def plot_cnv(self, save_plot=True):
        # log_inv = pd.read_csv('%s/femtic.cnv'%self.dir)
        log_inv = pd.read_csv('%s/femtic.cnv'%self.dir, delim_whitespace=True)
        
        fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('Iter#')
        ax1.set_ylabel('RMS', color=color)
        ax1.plot(log_inv['Iter#'], log_inv['RMS'], 'o-',color=color ,ms=3)
        ax1.tick_params(axis='y', labelcolor=color)        
        
        ax1.set_xlim(0,log_inv['Iter#'].max() +1)
        ax1.set_yticks(np.arange(0,np.ceil(log_inv['RMS'].max()) +1 ))
        # ax.set_yticks(np.arange(ylim[0], ylim[1]))
        ax1.grid(lw=0.3)

        ax2 = ax1.twinx() 

        color = 'tab:blue'
        ax2.set_ylabel('Roughness', color=color)  
        ax2.plot(log_inv['Iter#'], log_inv['Roughness'], 'o-',color=color ,ms=3)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout() 
        # plt.show
        
        if save_plot:
            plt.savefig('%s/log_rms.png'%(self.plot_dir),dpi=300,bbox_inches='tight')
            # plt.close('all')


    
        
        
    def z2rhophy(self, FREQ,ZR,ZI,dZ):
        
        ZR *= 10000/(4*np.pi)
        ZI *= - 10000/(4*np.pi)
        dZ *= 10000/(4*np.pi)
        
        # # calcul of apparent resistivity and phases
        rho = ((ZR**2+ZI**2)*0.2/(FREQ))
        phy = np.degrees(np.arctan2(ZI,ZR))
        # # calcul of errors
        drho = 2*rho*dZ / (((ZR**2+ZI**2)**0.5))
        dphy = np.degrees(0.5 * (drho/rho))
        log10_drho = (0.3772 * (dZ**2/(ZR**2+ZI**2)))**0.5

        return (rho, phy, drho, dphy,log10_drho)  
    
    
    
    def compute_rms(self, residuals):
        lists_res = [np.array(residual.values.tolist(), dtype=float).flatten() for residual in residuals]
        combined = np.array([item for sublist in lists_res for item in sublist], dtype=float)
        rms = (sum(combined**2) / (len(combined)))**0.5
        return rms




    # def compute_rms_breakDown(self):
        
    #     # total RMS
    #     filter_col = [col for col in self.resp if col.endswith(('Res'))]
    #     residuals = self.resp[filter_col]
    #     self.rms_total  = self.compute_rms(residuals)
        
    #     self.rms_breakDown = pd.DataFrame(columns=['StaID','Total','Zxx','Zxy','Zyx','Zyy'])
        
    #     # RMS / site
    #     rms_sites = []
    #     rms_sites_id = []
    #     filter_col = [col for col in self.resp if col.endswith(('StaID','Res'))]
    #     residuals = self.resp[filter_col]   
        
    #     for site in range(len(self.mt_ids)):
    #         residuals_site = residuals.loc[np.where(residuals['StaID'] == self.mt_ids[site])[0]]
    #         residuals_site = residuals_site.drop(['StaID'], axis=1)
    #         rms_site = self.compute_rms(residuals_site)
    #         rms_sites_id.append(self.mt_ids[site])
    #         rms_sites.append(rms_site)
    #     rms_sites = np.array(rms_sites)
        
        
    #     # RMS / site / component
    #     rms_sites_comps = []
    #     for site in range(len(self.mt_ids)):
    #         rms_comps = []
    #         residuals_site = residuals.loc[np.where(residuals['StaID'] == self.mt_ids[site])[0]]
    #         residuals_site = residuals_site.drop(['StaID'], axis=1)
    #         comp = ['Zxx','Zxy','Zyx','Zyy']
    #         for i in range(4):
    #             filter_col = [col for col in residuals_site if comp[i] in col]
    #             residuals_site_comp = residuals_site[filter_col]
    #             rms_comp = self.compute_rms(residuals_site_comp)
    #             rms_comps.append(rms_comp)
    #         rms_sites_comps.append(rms_comps)
    #     rms_sites_comps  = np.array(rms_sites_comps)

    #     self.rms_breakDown['StaID'] = rms_sites_id
    #     self.rms_breakDown['Total'] = rms_sites
    #     self.rms_breakDown[['Zxx','Zxy','Zyx','Zyy']] = rms_sites_comps



    def compute_rms_breakDown(self):
        
        # total RMS
        residuals = []
        if self.invert_Z:
            filter_col = [col for col in self.resp_Z if col.endswith(('Res'))]
            residuals_Z = self.resp_Z[filter_col]
            self.rms_Z  = self.compute_rms([residuals_Z])
            residuals.append(residuals_Z)
        if self.invert_VTF:
            filter_col = [col for col in self.resp_VTF if col.endswith(('Res'))]
            residuals_VTF = self.resp_VTF[filter_col]
            self.rms_VTF  = self.compute_rms([residuals_VTF])
            residuals.append(residuals_VTF)            
        if self.invert_PT:
           filter_col = [col for col in self.resp_PT if col.endswith(('Res'))]
           residuals_PT = self.resp_PT[filter_col]
           residuals.append(residuals_PT)

        self.rms_total  = self.compute_rms(residuals)
        
        
        if self.invert_Z:
        
            self.rms_breakDown_Z = pd.DataFrame(columns=['StaID','Total','Zxx','Zxy','Zyx','Zyy'])
        
            # RMS / site
            rms_sites = []
            rms_sites_id = []
            filter_col = [col for col in self.resp_Z if col.endswith(('StaID','Res'))]
            residuals = self.resp_Z[filter_col]   
            
            for site in range(len(self.ids_Z)):
                residuals_site = residuals.loc[np.where(residuals['StaID'] == self.ids_Z[site])[0]]
                residuals_site = residuals_site.drop(['StaID'], axis=1)
                rms_site = self.compute_rms([residuals_site])
                rms_sites_id.append(self.ids_Z[site])
                rms_sites.append(rms_site)
            rms_sites = np.array(rms_sites)
            
            
            # RMS / site / component
            rms_sites_comps = []
            for site in range(len(self.ids_Z)):
                rms_comps = []
                residuals_site = residuals.loc[np.where(residuals['StaID'] == self.ids_Z[site])[0]]
                residuals_site = residuals_site.drop(['StaID'], axis=1)
                comp = ['Zxx','Zxy','Zyx','Zyy']
                for i in range(4):
                    filter_col = [col for col in residuals_site if comp[i] in col]
                    residuals_site_comp = residuals_site[filter_col]
                    rms_comp = self.compute_rms([residuals_site_comp])
                    rms_comps.append(rms_comp)
                rms_sites_comps.append(rms_comps)
            rms_sites_comps  = np.array(rms_sites_comps)
    
            self.rms_breakDown_Z['StaID'] = rms_sites_id
            self.rms_breakDown_Z['Total'] = rms_sites
            self.rms_breakDown_Z[['Zxx','Zxy','Zyx','Zyy']] = rms_sites_comps


        if self.invert_VTF:
        
            self.rms_breakDown_VTF = pd.DataFrame(columns=['StaID','Total','Tzx','Tzy'])
        
            # RMS / site
            rms_sites = []
            rms_sites_id = []
            filter_col = [col for col in self.resp_VTF if col.endswith(('StaID','Res'))]
            residuals = self.resp_VTF[filter_col]   
            
            for site in range(len(self.ids_VTF)):
                residuals_site = residuals.loc[np.where(residuals['StaID'] == self.ids_VTF[site])[0]]
                residuals_site = residuals_site.drop(['StaID'], axis=1)
                rms_site = self.compute_rms([residuals_site])
                rms_sites_id.append(self.ids_VTF[site])
                rms_sites.append(rms_site)
            rms_sites = np.array(rms_sites)
            
            
            # RMS / site / component
            rms_sites_comps = []
            for site in range(len(self.ids_VTF)):
                rms_comps = []
                residuals_site = residuals.loc[np.where(residuals['StaID'] == self.ids_VTF[site])[0]]
                residuals_site = residuals_site.drop(['StaID'], axis=1)
                comp = ['Tzx','Tzy']
                for i in range(2):
                    filter_col = [col for col in residuals_site if comp[i] in col]
                    residuals_site_comp = residuals_site[filter_col]
                    rms_comp = self.compute_rms([residuals_site_comp])
                    rms_comps.append(rms_comp)
                rms_sites_comps.append(rms_comps)
            rms_sites_comps  = np.array(rms_sites_comps)
    
            self.rms_breakDown_VTF['StaID'] = rms_sites_id
            self.rms_breakDown_VTF['Total'] = rms_sites
            self.rms_breakDown_VTF[['Tzx','Tzy']] = rms_sites_comps



                # ind = np.where(self.ids_VTF[rx] == self.mt_coords['id'].values)[0][0]
                
                # file.write('%d  %d  %2.3f  %2.3f  \n'%(rx+1001, rx+1001, 
                #                                      self.mt_coords['north'][ind],
                #                                      self.mt_coords['east'][ind]))



    def read_distorsion(self):
        for root, dirs, files in os.walk(self.dir):
            for file in (files):
                if file.endswith('distortion_iter%d.dat'%self.nIter):
                    file_path = '%s/%s'%(self.dir, file)

        self.gd_data = pd.read_csv(file_path, skiprows=1 ,delim_whitespace=True,names=['Cxx','Cxy','Cyx','Cyy','act'])


    def calc_gd_strength(self):
        # following Adveeda et al., 2015
        self.gd_data['gd'] = 0
        for i in range(len(self.gd_data)):
            gd_matrix = np.array([[self.gd_data['Cxx'].iloc[i], self.gd_data['Cxy'].iloc[i]], [self.gd_data['Cyx'].iloc[i], self.gd_data['Cyy'].iloc[i]]])
            self.gd_data['gd'].iloc[i] = np.linalg.norm(gd_matrix - np.identity(2))


    def plot_distorsion_map(self,
                     ax, 
                     xlim = [-100, 100],ylim = [-100, 100], 
                     vmin=1,vmax=3,
                     save_plot = False):

        self.read_distorsion()
        self.calc_gd_strength()

        print('Plotting nRMSE map for the impedance Z...')
        pc=ax.scatter(self.mt_coords['east'],self.mt_coords['north'],c='w',
                #vmin=vmin,vmax=vmax,cmap=cmap,
                marker='o', s=1+self.gd_data['gd'], linewidths=0.5,edgecolors='k',)
        ax.scatter(10,-10,s=1, linewidths=0.5,edgecolors='k',)
        ax.scatter(10,-10,s=10, linewidths=0.5,edgecolors='k',)
        ax.set_title('Galvanic Distortion strength')
    
        #ax.axis('equal')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.grid(lw=0.1)
        ax.set_xlabel('East (km)')
        ax.set_ylabel('North (km)')
        #plt.tight_layout()
        
        if save_plot:
                plt.savefig('%s/gd_map.png'%(self.plot_dir),dpi=300,bbox_inches='tight')


    



    def plot_rms_map(self,
                     ax, 
                     Z  = True,
                     VTF = False,
                     PT = False,
                     xlim = [-100, 100],ylim = [-100, 100], 
                     vmin=1,vmax=3,
                     ms_size = 40,
                     save_plot = True):

        
        #plt.figure(1,figsize=figsize)
        cmap = cm.get_cmap('bone_r', 12)
        
        if Z:
            print('Plotting nRMSE map for the impedance Z...')
            pc=ax.scatter(self.mt_coords['east'],self.mt_coords['north'],c=self.rms_breakDown_Z['Total'],
                    vmin=vmin,vmax=vmax,cmap=cmap,
                    marker='o', s=ms_size, linewidths=0.5,edgecolors='k',)
            ax.set_title('nRMSE Map    Total nRMSE Z = %.2f'%self.rms_Z)
            
        if VTF:
            print('Plotting nRMSE map for the vetical transfer function VTF...')
            
            self.mt_coords_VTF = pd.DataFrame(columns=self.mt_coords.columns)
            print(self.mt_coords_VTF)
            for rx in range(len(self.ids_VTF)):
                ind = np.where(self.ids_VTF[rx] == self.mt_coords['id'].values)[0][0]
                #self.mt_coords_VTF = self.mt_coords_VTF.append(self.mt_coords.loc[ind])
                self.mt_coords_VTF  = pd.concat([self.mt_coords_VTF, pd.DataFrame([self.mt_coords.loc[ind]])], ignore_index=True)
            
            pc=ax.scatter(self.mt_coords_VTF['east'],self.mt_coords_VTF['north'],c=self.rms_breakDown_VTF['Total'],
                    vmin=vmin,vmax=vmax,cmap=cmap,
                    marker='o', s=ms_size, linewidths=0.5,edgecolors='k',)
            ax.set_title('nRMSE Map    Total nRMSE VTF = %.2f'%self.rms_VTF)
            
        if PT:
            pc=ax.scatter(self.mt_coords['east'],self.mt_coords['north'],c=self.rms_breakDown['Total'],
                    vmin=vmin,vmax=vmax,cmap=cmap,
                    marker='o', s=ms_size, linewidths=0.5,edgecolors='k',)
        plt.colorbar(pc, ax=ax)
        
        #ax.axis('equal')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.grid(lw=0.1)
        ax.set_xlabel('East (km)')
        ax.set_ylabel('North (km)')
        #plt.tight_layout()
        
        if save_plot:
            if Z:
                plt.savefig('%s/rmsZ_map.png'%(self.plot_dir),dpi=300,bbox_inches='tight')
            if VTF:
                plt.savefig('%s/rmsVTF_map.png'%(self.plot_dir),dpi=300,bbox_inches='tight')
    

    
    def plot_rms_map_components(self,figsize=(16,3.5), 
                                xlim = [-100, 100],ylim = [-100, 100], 
                                vmin=1,vmax=3,
                                ms_size = 40,
                                save_plot = True):

        print('Plotting nRMSE map ...')
        fig, axs = plt.subplots(1,4,sharey=True, figsize=figsize)
        cmap = cm.get_cmap('bone_r', 12)

        comps = ['Zxx','Zxy','Zyx','Zyy']
        for i, comp in enumerate(comps):
            
            pc=axs[i].scatter(self.mt_coords['east'],self.mt_coords['north'],c=self.rms_breakDown[comp],
                        vmin=vmin,vmax=vmax,cmap=cmap,
                        marker='o', s=ms_size,linewidths=0.5,edgecolors='k',)
            plt.colorbar(pc, ax=axs[i])
            axs[i].set_title('%s'%comps[i])
            axs[i].axis('equal')
            # plt.grid(lw=0.5)
            # plt.xlabel('East (km)')
            # plt.ylabel('North (km)')
        plt.tight_layout()    
        
        if save_plot:
            plt.savefig('%s/rms_map_comp.png'%(self.plot_dir),dpi=300,bbox_inches='tight')
            # plt.close('all')

    
    def plot_Z_fit(self, 
                 plot_Z = False,
                 xlim = [0,5],
                 ylim = [-6,1],
                 save_plot = True,
                 add_map_stats = False):
        
        
        if self.resp_Z is None:
            print("No Z data was inverted!")
            return 
        
        def format_fit(ax, 
                           xlim = [0,5], 
                           ylim = [-1,5], 
                           xlabel = 'Log$_{10}$ Period (sec)', 
                           ylabel = 'Log$_{10}$ $\\rho_{app}$ ($\Omega$m)'):
                
                ax.set_xlabel(xlabel)
                ax.set_xticks(np.arange(xlim[0],xlim[1] ))
                # ax.set_yticks(np.arange(ylim[0], ylim[1]))
                ax.set_ylim(ylim[0], ylim[1])
                ax.set_ylabel(ylabel)
                ax.grid(lw=0.3)
                ax.legend(loc='upper right',fontsize='xx-small',borderaxespad=-0.7,framealpha=1)
        
        
        print('Plotting inversion responses ...')
        for site in range(len(self.ids_Z)):
            print('   ...MT site ', self.ids_Z[site])
            
            data = self.resp_Z[self.resp_Z['StaID'] == self.ids_Z[site]]

            f = data['Freq[Hz]']                
            
            fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(8,6),
                            sharex='all')
            axs[0,0].set_title('%s'%self.ids_Z[site])
            
            if plot_Z: 
                
                self.plot_Z(axs[0,0],f,abs(data['Re(Zxx)_Obs']), data['Re(Zxx)_SD'],c='r',label='Re',obs=True)
                self.plot_Z(axs[0,0],f,abs(data['Im(Zxx)_Obs']), data['Im(Zxx)_SD'],c='b',label='Im',obs=True)
                self.plot_Z(axs[0,0],f,abs(data['Re(Zxx)_Cal']),c='r',obs=False)
                self.plot_Z(axs[0,0],f,abs(data['Im(Zxx)_Cal']),c='b',obs=False)
                format_fit(axs[0,0], xlim =xlim, ylim =ylim , xlabel = '', ylabel = 'Log$_{10}$ Zxx (ohm)')
                
                
                self.plot_Z(axs[0,1],f,abs(data['Re(Zxy)_Obs']), data['Re(Zxy)_SD'],c='r',label='Re',obs=True)
                self.plot_Z(axs[0,1],f,abs(data['Im(Zxy)_Obs']), data['Im(Zxy)_SD'],c='b',label='Im',obs=True)
                self.plot_Z(axs[0,1],f,abs(data['Re(Zxy)_Cal']),c='r',obs=False)
                self.plot_Z(axs[0,1],f,abs(data['Im(Zxy)_Cal']),c='b',obs=False)
                format_fit(axs[0,1], xlim =xlim, ylim =ylim , ylabel = 'Log$_{10}$ Zxy (ohm)')
                
                    
                self.plot_Z(axs[1,0],f,abs(data['Re(Zyx)_Obs']), data['Re(Zyx)_SD'],c='r',label='Re',obs=True)
                self.plot_Z(axs[1,0],f,abs(data['Im(Zyx)_Obs']), data['Im(Zyx)_SD'],c='b',label='Im',obs=True)
                self.plot_Z(axs[1,0],f,abs(data['Re(Zyx)_Cal']),c='r',obs=False)
                self.plot_Z(axs[1,0],f,abs(data['Im(Zyx)_Cal']),c='b',obs=False)
                format_fit(axs[1,0], xlim =xlim, ylim =ylim, xlabel = '' , ylabel = 'Log$_{10}$ Zyx (ohm)')
                
                
                self.plot_Z(axs[1,1],f,abs(data['Re(Zyy)_Obs']), data['Re(Zyy)_SD'],c='r',label='Re',obs=True)
                self.plot_Z(axs[1,1],f,abs(data['Im(Zyy)_Obs']), data['Im(Zyy)_SD'],c='b',label='Im',obs=True)
                self.plot_Z(axs[1,1],f,abs(data['Re(Zyy)_Cal']),c='r',obs=False)
                self.plot_Z(axs[1,1],f,abs(data['Im(Zyy)_Cal']),c='b',obs=False)
                format_fit(axs[1,1], xlim =xlim, ylim =ylim , ylabel = 'Log$_{10}$ Zyy (ohm)')
                
                plt.tight_layout()
                
                
                if save_plot:
                    plt.savefig('%s/Z_%s.png'%(self.plot_dir,self.ids_Z[site]),dpi=300,bbox_inches='tight')
                    plt.close('all')
                
            else:
 
                rhoXX, phyXX, drhoXX, dphyXX, log_drhoXX = self.z2rhophy(data['Freq[Hz]'], data['Re(Zxx)_Obs'],data['Im(Zxx)_Obs'],data['Re(Zxx)_SD'])
                rhoXY, phyXY, drhoXY, dphyXY, log_drhoXY = self.z2rhophy(data['Freq[Hz]'], data['Re(Zxy)_Obs'],data['Im(Zxy)_Obs'],data['Re(Zxy)_SD'])
                rhoYX, phyYX, drhoYX, dphyYX, log_drhoYX = self.z2rhophy(data['Freq[Hz]'], data['Re(Zyx)_Obs'],data['Im(Zyx)_Obs'],data['Re(Zyx)_SD'])
                rhoYY, phyYY, drhoYY, dphyYY, log_drhoYY = self.z2rhophy(data['Freq[Hz]'], data['Re(Zyy)_Obs'],data['Im(Zyy)_Obs'],data['Re(Zyy)_SD'])
                    
                rhoXX_calc, phyXX_calc, _, _, _ = self.z2rhophy(data['Freq[Hz]'], data['Re(Zxx)_Cal'],data['Im(Zxx)_Cal'],data['Re(Zxx)_SD'])
                rhoXY_calc, phyXY_calc, _, _, _ = self.z2rhophy(data['Freq[Hz]'], data['Re(Zxy)_Cal'],data['Im(Zxy)_Cal'],data['Re(Zxy)_SD'])
                rhoYX_calc, phyYX_calc, _, _, _ = self.z2rhophy(data['Freq[Hz]'], data['Re(Zyx)_Cal'],data['Im(Zyx)_Cal'],data['Re(Zyx)_SD'])
                rhoYY_calc, phyYY_calc, _, _, _ = self.z2rhophy(data['Freq[Hz]'], data['Re(Zyy)_Cal'],data['Im(Zyy)_Cal'],data['Re(Zyy)_SD'])



                self.plot_rho(axs[0,0],f,rhoXY, log_drhoXY,c='r',label='xy',obs=True)
                self.plot_rho(axs[0,0],f,rhoYX, log_drhoYX,c='b',label='yx',obs=True)
                self.plot_rho(axs[0,0],f,rhoXY_calc,c='r',label='xy',obs=False)
                self.plot_rho(axs[0,0],f,rhoYX_calc,c='b',label='yx',obs=False)
                format_fit(axs[0,0], xlim =xlim, ylim = [-1,5], xlabel = '' , ylabel = 'Log$_{10}$ $\\rho_{app}$ ($\Omega$m)')
                
                
                self.plot_phy(axs[1,0],f,phyXY, dphyXY,c='r',label='xy',obs=True)
                self.plot_phy(axs[1,0],f,phyYX, dphyYX,c='b',label='yx',obs=True)
                self.plot_phy(axs[1,0],f,phyXY_calc,c='r',obs=False)
                self.plot_phy(axs[1,0],f,phyYX_calc,c='b',obs=False)
                format_fit(axs[1,0], xlim =xlim, ylim =[-180,180] , ylabel = 'Phase ($^\circ$)')
                
                    
                self.plot_rho(axs[0,1],f,rhoXX, log_drhoXX,c='r',label='xx',obs=True)
                self.plot_rho(axs[0,1],f,rhoYY, log_drhoYY,c='b',label='yy',obs=True)
                self.plot_rho(axs[0,1],f,rhoXX_calc,c='r',obs=False)
                self.plot_rho(axs[0,1],f,rhoYY_calc,c='b',obs=False)
                format_fit(axs[0,1], xlim =xlim, ylim =[-1,5] , xlabel = '', ylabel = '')
                
                
                self.plot_phy(axs[1,1],f,phyXX, dphyXX,c='r',label='xx',obs=True)
                self.plot_phy(axs[1,1],f,phyYY, dphyYY,c='b',label='yy',obs=True)
                self.plot_phy(axs[1,1],f,phyXX_calc,c='r',obs=False)
                self.plot_phy(axs[1,1],f,phyYY_calc,c='b',obs=False)
                format_fit(axs[1,1], xlim =xlim, ylim = [-180,180] , ylabel = '')

                plt.tight_layout()

                    
            if add_map_stats:
                plt.subplots_adjust(right=0.7)
                ax_map = fig.add_axes([0.75,.55,.2,.3])
                ax_stats = fig.add_axes([0.75,.1,.2,.3])
                
                self.plot_loc_map(self.ids_Z[site], ax_map)
                self.plot_inv_stats_Z(self.ids_Z[site], ax_stats)
  
                
            if save_plot:
                plt.savefig('%s/Z_%s.png'%(self.plot_dir,self.ids_Z[site]),dpi=100,bbox_inches='tight')
                plt.close('all')
                

    def plot_loc_map(self, site_id, ax):
        ax.axis('equal')
        ax.set_xticks([]) 
        ax.set_yticks([])
        ax.plot(self.mt_coords['east'], self.mt_coords['north'], 'ko',ms=2)
        coord_site = self.mt_coords[self.mt_coords['id'] == site_id]
        ax.plot(coord_site['east'], coord_site['north'], 'ro',ms=3)
        
        
        
    def plot_inv_stats_Z(self, site_id, ax):
        
        ax.axis("off")
        ax.set_xlim(0,10)
        ax.set_ylim(0,10)
        
        rms_site = self.rms_breakDown_Z[self.rms_breakDown_Z['StaID'] == site_id]
        
        ax.text(0, 9, 'RMS site = %.2f'%rms_site['Total'].values)
        ax.text(0, 7, 'RMS Zxx  = %.2f'%rms_site['Zxx'].values)
        ax.text(0, 5, 'RMS Zxy  = %.2f'%rms_site['Zxy'].values)
        ax.text(0, 3, 'RMS Zyx  = %.2f'%rms_site['Zyx'].values)
        ax.text(0, 1, 'RMS Zyy  = %.2f'%rms_site['Zyy'].values)
 

    def plot_inv_stats_VTF(self, site_id, ax):
        
        ax.axis("off")
        ax.set_xlim(0,10)
        ax.set_ylim(0,10)
        
        rms_site = self.rms_breakDown_VTF[self.rms_breakDown_VTF['StaID'] == site_id]
        
        ax.text(0, 9, 'RMS VTF site = %.2f'%rms_site['Total'].values)
        ax.text(0, 7, 'RMS Tzx  = %.2f'%rms_site['Tzx'].values)
        ax.text(0, 5, 'RMS Tzy  = %.2f'%rms_site['Tzy'].values)
     
 
        
    def plot_induction_arrows(self,
                              frequencies,
                              # real = True,
                              # imag = False,
                              inv_response = False,
                              scale = 1,
                              xlim = [-100, 100],
                              ylim = [-100, 100], 
                              save_plot = True):
        
        if self.resp_VTF is None:
            print("No VTF data was inverted!")
            return 
        
        from scipy.spatial.distance import pdist
        dist_sites = min(pdist(np.array([self.mt_coords['east'],self.mt_coords['north']])).T)
        # print(dist_sites)
        scale = 0.5*dist_sites * scale
        
        self.mt_coords_VTF = pd.DataFrame(columns=self.mt_coords.columns)
        for rx in range(len(self.ids_VTF)):
            ind = np.where(self.ids_VTF[rx] == self.mt_coords['id'].values)[0][0]
            #self.mt_coords_VTF = self.mt_coords_VTF.append(self.mt_coords.loc[ind])
            self.mt_coords_VTF  = pd.concat([self.mt_coords_VTF, pd.DataFrame([self.mt_coords.loc[ind]])], ignore_index=True)

        fig, axs = plt.subplots(2, len(frequencies),sharex=False, sharey=False,figsize=(6*len(frequencies),10))

        for fr in range(len(frequencies)):
    
            ia = self.resp_VTF.iloc[np.where(self.resp_VTF['Freq[Hz]'] == frequencies[fr])[0]]
            ia = ia.reset_index()
        
            for i in range(len(ia)):
                
                axs[0,fr].scatter(self.mt_coords['east'],self.mt_coords['north'],c='k', marker='o', s=1)
                axs[1,fr].scatter(self.mt_coords['east'],self.mt_coords['north'],c='k', marker='o', s=1)
        
                x = self.mt_coords_VTF.iloc[np.where(self.mt_coords_VTF['id'] == ia['StaID'][i])[0]]['east'].values[0]
                y = self.mt_coords_VTF.iloc[np.where(self.mt_coords_VTF['id'] == ia['StaID'][i])[0]]['north'].values[0]
                                
                #axs[0].annotate("",xy=(x + ia['Re(Tzx)_Obs'][i] * scale , y + ia['Re(Tzy)_Obs'][i] * scale),
                #                xytext=(x,y),
                 #               arrowprops=dict(facecolor='black',shrink=0.001,lw=0.01))
        
                axs[0,fr].arrow(x,y,ia['Re(Tzy)_Obs'][i] * scale ,ia['Re(Tzx)_Obs'][i] * scale, width=0.1,
                            ec='k', fc='k')
                axs[0,fr].set_title('Re Induction arrows T = %.3f s'%(1/frequencies[fr]))
        
                if inv_response:
                    axs[0,fr].arrow(x,y,ia['Re(Tzy)_Cal'][i] * scale ,ia['Re(Tzx)_Cal'][i] * scale, width=0.1, alpha=0.8,
                            ec='r', fc='r')
        
                axs[1,fr].arrow(x,y,ia['Im(Tzy)_Obs'][i] * scale ,ia['Im(Tzx)_Obs'][i] * scale, width=0.1,
                            ec='k', fc='k')
                axs[1,fr].set_title('Imag Induction arrows T = %.3f s'%(1/frequencies[fr]))
        
                if inv_response:
                    axs[1,fr].arrow(x,y,ia['Im(Tzy)_Cal'][i] * scale ,ia['Im(Tzx)_Cal'][i] * scale, width=0.1, alpha=0.8,
                        ec='r', fc='r')
                
                axs[0,fr].set_xlim(xlim)
                axs[0,fr].set_ylim(ylim)
                axs[0,fr].grid(lw=0.1)
                axs[0,fr].set_xlabel('East (km)')
                axs[0,fr].set_ylabel('North (km)')
                
                axs[1,fr].set_xlim(xlim)
                axs[1,fr].set_ylim(ylim)
                axs[1,fr].grid(lw=0.1)
                axs[1,fr].set_xlabel('East (km)')
                axs[1,fr].set_ylabel('North (km)')
                
        if save_plot:
            plt.savefig('%s/IA.png'%(self.plot_dir),dpi=300,bbox_inches='tight')
            plt.close('all')
        
        
    def plot_VTF_fit(self,   xlim = [0,5],
                         ylim = [-6,1],
                         save_plot = True,
                         add_map_stats = False):
        
        
        if self.resp_VTF is None:
            print("No VTF data was inverted!")
            return 
        
        def format_fit(ax, 
                           xlim = [0,5], 
                           ylim = [-1.1,1.1], 
                           xlabel = 'Log$_{10}$ Period (sec)', 
                           ylabel = 'Tzx'):
                
                ax.set_xlabel(xlabel)
                ax.set_xticks(np.arange(xlim[0],xlim[1] ))
                # ax.set_yticks(np.arange(ylim[0], ylim[1]))
                ax.set_ylim(ylim[0], ylim[1])
                ax.set_ylabel(ylabel)
                ax.grid(lw=0.3)
                ax.legend(loc='upper right',fontsize='xx-small',borderaxespad=-0.7,framealpha=1)
        
        
        print('Plotting VTF responses ...')
        for site in range(len(self.ids_VTF)):
            print('   ...MT site ', self.ids_VTF[site])
            
            data = self.resp_VTF[self.resp_VTF['StaID'] == self.ids_VTF[site]]

            f = data['Freq[Hz]']                
            
            fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(5,5),
                            sharex='all')
            axs[0].set_title('%s'%self.ids_VTF[site])
            
            self.plot_tz(axs[0], f, data['Re(Tzx)_Obs'], c='r',label='Re',obs=True)
            self.plot_tz(axs[0], f, data['Im(Tzx)_Obs'], c='b',label='Im',obs=True)
            self.plot_tz(axs[0], f, data['Re(Tzx)_Cal'], c='r',obs=False)
            self.plot_tz(axs[0], f, data['Im(Tzx)_Cal'], c='b',obs=False)
            format_fit(axs[0], xlim =xlim, ylim =ylim , xlabel = '', ylabel = 'Tzx')
   
            self.plot_tz(axs[1], f, data['Re(Tzy)_Obs'], c='r',label='Re',obs=True)
            self.plot_tz(axs[1], f, data['Im(Tzy)_Obs'], c='b',label='Im',obs=True)
            self.plot_tz(axs[1], f, data['Re(Tzy)_Cal'], c='r',obs=False)
            self.plot_tz(axs[1], f, data['Im(Tzy)_Cal'], c='b',obs=False)
            format_fit(axs[1], xlim =xlim, ylim =ylim , xlabel = '', ylabel = 'Tzy')

            plt.tight_layout()

            if add_map_stats:
                plt.subplots_adjust(right=0.7)
                ax_map = fig.add_axes([0.75,.55,.2,.3])
                ax_stats = fig.add_axes([0.75,.1,.2,.3])
                
                self.plot_loc_map(self.ids_VTF[site], ax_map)
                self.plot_inv_stats_VTF(self.ids_VTF[site], ax_stats)
  
                
            if save_plot:
                plt.savefig('%s/VTF_%s.png'%(self.plot_dir,self.ids_VTF[site]),dpi=100,bbox_inches='tight')
                plt.close('all')        
        
    
    def plot_Z(self, ax,f,Z,dZ=None,c='r',label='xy',obs=True):
        if obs:
            ax.errorbar(np.log10(1/f), np.log10(Z), yerr = np.log10(Z+dZ) - np.log10(Z-dZ), 
                            fmt='%s.'%c,label= label,zorder=32, 
                            elinewidth=0.6,markersize=8 ,
                            capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
        else:
            ax.plot(np.log10(1/f), np.log10(Z),c=c)
            
            
    
    def plot_rho(self, ax,f,rho,log_drho=None,c='r',label='xy',obs=True):
        if obs:
            ax.errorbar(np.log10(1/f), np.log10(rho), yerr = log_drho, 
                            fmt='%s.'%c,label= label,zorder=32, 
                            elinewidth=0.6,markersize=8 ,
                            capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
        else:
            ax.plot(np.log10(1/f), np.log10(rho),c=c)
    
    
    
    def plot_phy(self, ax,f,phy,dphy=None,c='r',label='xy',obs=True):
        if obs:
            ax.errorbar(np.log10(1/f), phy, yerr = dphy, 
                            fmt='%s.'%c,label= label,zorder=32, 
                            elinewidth=0.6,markersize=8 ,
                            capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
        else:
            ax.plot(np.log10(1/f), phy,c=c)
    
    
    def plot_tz(self, ax,f,tz,dtz=None,c='r',label='tzx',obs=True):
        if obs:
            ax.errorbar(np.log10(1/f), tz, yerr = dtz, 
                            fmt='%s.'%c,label= label,zorder=32, 
                            elinewidth=0.6,markersize=8 ,
                            capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
        else:
            ax.plot(np.log10(1/f), tz,c=c)  
        
 