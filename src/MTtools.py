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
    
import numpy as np
import re
import pandas as pd
import os
from math import isnan
from datetime import datetime
import csv
import matplotlib.pyplot as plt
import copy


def read(file_name):
    
    """
    Input:
    - complete path to the EDI file
    
    Output:
    - data (EDI pandas DF format)
    - site_id (string)
    - coord (panda DF of 7 values (['Lat_deg']['Lat_min']['Lat_sec']['Long_deg']['Long_min']['Long_sec']['Elev'])    
    """   


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





def plot_rho(ax,f,rho,log_drho=None,c='r',label='xy',obs=True):
    if obs:
        ax.errorbar(np.log10(1/f), np.log10(rho), yerr = log_drho, 
                        fmt='%s.'%c,label= label,zorder=32, 
                        elinewidth=0.6,markersize=8 ,
                        capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
    else:
        ax.plot(np.log10(1/f), np.log10(rho),c=c)



def plot_phy(ax,f,phy,dphy=None,c='r',label='xy',obs=True):
    if obs:
        ax.errorbar(np.log10(1/f), phy, yerr = dphy, 
                        fmt='%s.'%c,label= label,zorder=32, 
                        elinewidth=0.6,markersize=8 ,
                        capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
    else:
        ax.plot(np.log10(1/f), phy,c=c)


def plot_tz(ax,f,tz,dtz=None,c='r',label='tzx',obs=True):
    if obs:
        ax.errorbar(np.log10(1/f), tz, yerr = dtz, 
                        fmt='%s.'%c,label= label,zorder=32, 
                        elinewidth=0.6,markersize=8 ,
                        capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
    else:
            ax.plot(np.log10(1/f), tz,c=c)  


def plot_loc_map(femticObj, site_id, ax):
    ax.axis('equal')
    ax.set_xticks([]) 
    ax.set_yticks([])
    ax.plot(femticObj.mt_coords['east'], femticObj.mt_coords['north'], 'ko',ms=2)
    coord_site = femticObj.mt_coords[femticObj.mt_coords['id'] == site_id]
    ax.plot(coord_site['east'], coord_site['north'], 'ro',ms=3)
        


def z2rhophy(FREQ,ZR,ZI,dZ):
    
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




def plot_Z(femticObj, plot_Z = False,
                xlim = [0,5],
                ylim = [-6,1],
                save_plot = True,
                add_map = False):
    
    
    if femticObj.data_Z is None:
        print("No Z data was inverted!")
        return 
    
    def format_fit(ax, 
                        xlim = [-2.5,4], 
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
    
    
    print('Plotting observed data ...')
    for site in range(len(femticObj.ids_Z)):
        print('   ...MT site ', femticObj.ids_Z[site])
        
        #data = self.resp_Z[self.resp_Z['StaID'] == self.ids_Z[site]]
        data = copy.deepcopy(femticObj.data_Z[site])

        f = data['FREQ']                
        
        fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(8,6),
                        sharex='all')
        axs[0,0].set_title('%s'%femticObj.ids_Z[site])
        
        if plot_Z: 
            
            plot_Z(axs[0,0],f,abs(data['Re(Zxx)_Obs']), data['Re(Zxx)_SD'],c='r',label='Re',obs=True)
            plot_Z(axs[0,0],f,abs(data['ZXXI']), data['Im(Zxx)_SD'],c='b',label='Im',obs=True)
            plot_Z(axs[0,0],f,abs(data['Re(Zxx)_Cal']),c='r',obs=False)
            plot_Z(axs[0,0],f,abs(data['Im(Zxx)_Cal']),c='b',obs=False)
            format_fit(axs[0,0], xlim =xlim, ylim =ylim , xlabel = '', ylabel = 'Log$_{10}$ Zxx (ohm)')
            
            
            plot_Z(axs[0,1],f,abs(data['Re(Zxy)_Obs']), data['Re(Zxy)_SD'],c='r',label='Re',obs=True)
            plot_Z(axs[0,1],f,abs(data['Im(Zxy)_Obs']), data['Im(Zxy)_SD'],c='b',label='Im',obs=True)
            plot_Z(axs[0,1],f,abs(data['Re(Zxy)_Cal']),c='r',obs=False)
            plot_Z(axs[0,1],f,abs(data['Im(Zxy)_Cal']),c='b',obs=False)
            format_fit(axs[0,1], xlim =xlim, ylim =ylim , ylabel = 'Log$_{10}$ Zxy (ohm)')
            
                
            plot_Z(axs[1,0],f,abs(data['Re(Zyx)_Obs']), data['Re(Zyx)_SD'],c='r',label='Re',obs=True)
            plot_Z(axs[1,0],f,abs(data['Im(Zyx)_Obs']), data['Im(Zyx)_SD'],c='b',label='Im',obs=True)
            plot_Z(axs[1,0],f,abs(data['Re(Zyx)_Cal']),c='r',obs=False)
            plot_Z(axs[1,0],f,abs(data['Im(Zyx)_Cal']),c='b',obs=False)
            format_fit(axs[1,0], xlim =xlim, ylim =ylim, xlabel = '' , ylabel = 'Log$_{10}$ Zyx (ohm)')
            
            
            plot_Z(axs[1,1],f,abs(data['Re(Zyy)_Obs']), data['Re(Zyy)_SD'],c='r',label='Re',obs=True)
            plot_Z(axs[1,1],f,abs(data['Im(Zyy)_Obs']), data['Im(Zyy)_SD'],c='b',label='Im',obs=True)
            plot_Z(axs[1,1],f,abs(data['Re(Zyy)_Cal']),c='r',obs=False)
            plot_Z(axs[1,1],f,abs(data['Im(Zyy)_Cal']),c='b',obs=False)
            format_fit(axs[1,1], xlim =xlim, ylim =ylim , ylabel = 'Log$_{10}$ Zyy (ohm)')
            
            plt.tight_layout()
            
            
            if save_plot:
                plt.savefig('%s/Z_%s.png'%(femticObj.plot_dir,femticObj.ids_Z[site]),dpi=300,bbox_inches='tight')
                plt.close('all')
            
        else:

            rhoXX, phyXX, drhoXX, dphyXX, log_drhoXX = z2rhophy(data['FREQ'], data['ZXXR'],data['ZXXI'],data['ZXX.VAR']**0.5)
            rhoXY, phyXY, drhoXY, dphyXY, log_drhoXY = z2rhophy(data['FREQ'], data['ZXYR'],data['ZXYI'],data['ZXY.VAR']**0.5)
            rhoYX, phyYX, drhoYX, dphyYX, log_drhoYX = z2rhophy(data['FREQ'], data['ZYXR'],data['ZYXI'],data['ZYX.VAR']**0.5)
            rhoYY, phyYY, drhoYY, dphyYY, log_drhoYY = z2rhophy(data['FREQ'], data['ZYYR'],data['ZYYI'],data['ZYY.VAR']**0.5)
            
            plot_rho(axs[0,0],f,rhoXY, log_drhoXY,c='r',label='xy',obs=True)
            plot_rho(axs[0,0],f,rhoYX, log_drhoYX,c='b',label='yx',obs=True)
            format_fit(axs[0,0], xlim =xlim, ylim = [-1,5], xlabel = '' , ylabel = 'Log$_{10}$ $\\rho_{app}$ ($\Omega$m)')
            
            
            plot_phy(axs[1,0],f,phyXY, dphyXY,c='r',label='xy',obs=True)
            plot_phy(axs[1,0],f,phyYX, dphyYX,c='b',label='yx',obs=True)
            format_fit(axs[1,0], xlim =xlim, ylim =[-180,180] , ylabel = 'Phase ($^\circ$)')
            
                
            plot_rho(axs[0,1],f,rhoXX, log_drhoXX,c='r',label='xx',obs=True)
            plot_rho(axs[0,1],f,rhoYY, log_drhoYY,c='b',label='yy',obs=True)
            format_fit(axs[0,1], xlim =xlim, ylim =[-1,5] , xlabel = '', ylabel = '')
            
            
            plot_phy(axs[1,1],f,phyXX, dphyXX,c='r',label='xx',obs=True)
            plot_phy(axs[1,1],f,phyYY, dphyYY,c='b',label='yy',obs=True)
            format_fit(axs[1,1], xlim =xlim, ylim = [-180,180] , ylabel = '')

            plt.tight_layout()


        if add_map:
            plt.subplots_adjust(right=0.7)
            ax_map = fig.add_axes([0.75,.55,.2,.3])            
            plot_loc_map(femticObj, femticObj.ids_Z[site], ax_map)
            
        if save_plot:
            if not os.path.isdir(femticObj.plot_dir):
                os.mkdir(femticObj.plot_dir)

            plt.savefig('%s/Z_%s.png'%(femticObj.plot_dir,femticObj.ids_Z[site]),dpi=100,bbox_inches='tight')
            plt.close('all')
            



def plot_VTF(femticObj,   xlim = [0,5],
                        ylim = [-6,1],
                        save_plot = True,
                        add_map = False):
    
    
    if femticObj.data_VTF is None:
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
    for site in range(len(femticObj.ids_VTF)):
        print('   ...MT site ', femticObj.ids_VTF[site])
        #data = femticObj.resp_VTF[femticObj.resp_VTF['StaID'] == femticObj.ids_VTF[site]]
        data = copy.deepcopy(femticObj.data_VTF[site])
        f = data['FREQ']                
        
        fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(5,5),
                        sharex='all')
        axs[0].set_title('%s'%femticObj.ids_VTF[site])
        
        plot_tz(axs[0], f, data['TXR.EXP'], dtz=data['TXVAR.EXP']**0.5, c='r',label='Re',obs=True)
        plot_tz(axs[0], f, data['TXI.EXP'], dtz=data['TXVAR.EXP']**0.5, c='b',label='Im',obs=True)
        #plot_tz(axs[0], f, data['Re(Tzx)_Cal'], c='r',obs=False)
        #plot_tz(axs[0], f, data['Im(Tzx)_Cal'], c='b',obs=False)
        format_fit(axs[0], xlim =xlim, ylim =ylim , xlabel = '', ylabel = 'Tzx')

        plot_tz(axs[1], f, data['TYR.EXP'], dtz=data['TYVAR.EXP']**0.5, c='r',label='Re',obs=True)
        plot_tz(axs[1], f, data['TYI.EXP'], dtz=data['TYVAR.EXP']**0.5, c='b',label='Im',obs=True)
        #plot_tz(axs[1], f, data['Re(Tzy)_Cal'], c='r',obs=False)
        #plot_tz(axs[1], f, data['Im(Tzy)_Cal'], c='b',obs=False)
        format_fit(axs[1], xlim =xlim, ylim =ylim , xlabel = '', ylabel = 'Tzy')

        plt.tight_layout()

        if add_map:
            plt.subplots_adjust(right=0.7)
            ax_map = fig.add_axes([0.75,.55,.2,.3])
            plot_loc_map(femticObj, femticObj.ids_VTF[site], ax_map)

            
        if save_plot:
            plt.savefig('%s/VTF_%s.png'%(femticObj.plot_dir,femticObj.ids_VTF[site]),dpi=100,bbox_inches='tight')
            plt.close('all')        
    