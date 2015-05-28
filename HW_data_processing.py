# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 15:47:01 2015

@author: jura
"""

from __future__ import division
import matplotlib.pyplot as plt
import HotWire_calc as hw
import os
import numpy as np
from multiprocessing import Pool

#######################x
U0=3.68 #[m/s] Uniform velocity 
wire_temper=260 #[°C] Wire temperature 
processor_no=4 #Number of processor cores on a PC
#######################x

prumery=[]
pozice_all=[]
cta=[]
rozmery=[]



#-------Gets all the files in all the subdirectories---------
directory='.\\Viric_data'
file_paths=[]
ax=[]
vir=[]
natoceni=[]

for root, directories, files in os.walk(directory):
  for filename in files:
    # Join the two strings in order to form the full filepath.
    if '.txt' in filename and 'viric_' in root:
      filepath = os.path.join(root, filename)
      file_paths.append(filepath)  # Add it to the list.
      ax.append(root.split('\\')[-1]) #Get all the axial distances
      vir.append(root.split('\\')[-2]) #Get all the swirler types
      if '_' in filename:
        natoceni.append(filename.split('_')[-1][:-4]) # Get all the rotation angles
#----------------------------------------------------------


#------Here you can restrict range of files to be processed--------------
file_paths=file_paths[:31]
ax=ax[:31]
vir=vir[:31]
natoceni=natoceni[:31]
#------------------------------------------------------------------------

#print file_paths

##---Gets all the data from all files------
#def read_cta(filenm):
#  return hw.HotWireMeasurement(filenm,wire_temper)
#
#p=Pool(processes=processor_no)
#print len(file_paths)
#cta=p.map(read_cta,file_paths)


for f in file_paths:
  #Načte data pro daný soubor f, uloží je do třídy HotWireMeasurement a třídu uloží do pole cta
  cta.append(hw.HotWireMeasurement(f,wire_temper=wire_temper)) 

for c in cta:
  print c
  peaks=c.fft_analysis()
  #print 'Significant frequencies: ',peaks,'Hz'
  
##-----------------------------------------
  
#Eliminates repeates
ax=np.unique(ax)
vir=np.unique(vir)
natoceni=np.unique(natoceni)
#print natoceni

#swirl number

nazev_ind=[]
for v in vir:
  for ax_pozice in ax:
    if ax_pozice == '0,1D':
      ax_pozice='01D'
    elif ax_pozice == '1D':
      ax_pozice='10D'
      
    for nat in natoceni:
      txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind.append(n)
    ccs=[cta[c] for c  in nazev_ind ]
    phi=[]
    for nzv in nazev_ind:
        phi.append(int(file_paths[nzv].split('_')[-1][:-4]))
    
    S=hw.swirl_number(ccs, 150)
    #phi=[int(c) for c in natoceni]
    print file_paths[nazev_ind[-1]][:5]+': Swirl nuber:',S
    if len(phi) == len(ccs):
        V=hw.flow_rate(ccs,150,phi)
        print file_paths[nazev_ind[-1]][:5]+': Flow rate:',V,"m3/s"
    else:
        print "Cannot calculate flow rate: phi list has different length from ccs list!"
    nazev_ind=[]

nazev_ind=[]
for v in vir:
    for ax_pozice in ax:
        if ax_pozice == '0,1D':
            ax_pozice='01D'
        elif ax_pozice == '1D':
            ax_pozice='10D'
        
        for nat in natoceni:
            txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
            n=[n for n,s in enumerate(file_paths) if txt in s]
            if len(n) > 0:
                n=n[0]
                nazev_ind.append(n)
        # ccs=[cta[c] for c in nazev_ind ]
        #U=[cta[c].U_mag for c in nazev_ind ]
        U=cta[nazev_ind[0]].U_mag
        if len(nazev_ind)>1:
            for c in nazev_ind[1:]:
                U=np.vstack((U,cta[c].U_mag))
            fig, axis = plt.subplots(figsize=(5,5),subplot_kw=dict(projection='polar'))
            phi=[]
            for nzv in nazev_ind:
                phi.append(int(file_paths[nzv].split('_')[-1][:-4]))
            #phi=[int(c) for c in natoceni]
            phi=np.radians(phi)
            
            rad=150-cta[nazev_ind[-1]].x
            
            phis=np.array(phi[:-2])
            # rads=np.array(rad)
            Us=np.copy(U[:-2])
            for i in range(7):
                phis=np.append(phis,phi[:-2]+np.pi*(i+1)/4)
                Us=np.vstack((Us,U[:-2]))
            
            phis=np.append(phis,phi[-1]+np.pi*7/4)
            Us=np.vstack((Us,U[0]))
            rad, phi=np.meshgrid(rad,phis)
            CS=axis.contourf(phi,rad,Us)
            plt.colorbar(CS)
            #plt.rcParams.update({'font.size': 8})
            plt.savefig(file_paths[nazev_ind[-1]][:-16]+'U_mag_2D.png',dpi=300,bbox_inches='tight')
            nazev_ind=[]
            
            plt.close()
        
nazev_ind=[]
for v in vir:
    for ax_pozice in ax:
        if ax_pozice == '0,1D':
            ax_pozice='01D'
        elif ax_pozice == '1D':
            ax_pozice='10D'
        
        for nat in natoceni:
            txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
            n=[n for n,s in enumerate(file_paths) if txt in s]
            if len(n) > 0:
                n=n[0]
                nazev_ind.append(n)
        # ccs=[cta[c] for c in nazev_ind ]
        #U=[cta[c].U_aver for c in nazev_ind ]
        U=cta[nazev_ind[0]].U_aver
        if len(nazev_ind)>1:
            for c in nazev_ind[1:]:
                U=np.vstack((U,cta[c].U_aver))
            fig, axis = plt.subplots(figsize=(5,5),subplot_kw=dict(projection='polar'))
            phi=[]
            for nzv in nazev_ind:
                phi.append(int(file_paths[nzv].split('_')[-1][:-4]))
            #phi=[int(c) for c in natoceni]
            phi=np.radians(phi)
            
            rad=150-cta[nazev_ind[-1]].x
            
            phis=np.array(phi[:-2])
            # rads=np.array(rad)
            Us=np.copy(U[:-2])
            for i in range(7):
                phis=np.append(phis,phi[:-2]+np.pi*(i+1)/4)
                Us=np.vstack((Us,U[:-2]))
            
            phis=np.append(phis,phi[-1]+np.pi*7/4)
            Us=np.vstack((Us,U[0]))
            rad, phi=np.meshgrid(rad,phis)
            CS=axis.contourf(phi,rad,Us)
            plt.colorbar(CS)
            #plt.rcParams.update({'font.size': 8})
            plt.savefig(file_paths[nazev_ind[-1]][:-16]+'U_aver_2D.png',dpi=300,bbox_inches='tight')
            nazev_ind=[]
            
            plt.close()
            
nazev_ind=[]
for v in vir:
    for ax_pozice in ax:
        if ax_pozice == '0,1D':
            ax_pozice='01D'
        elif ax_pozice == '1D':
            ax_pozice='10D'
        
        for nat in natoceni:
            txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
            n=[n for n,s in enumerate(file_paths) if txt in s]
            if len(n) > 0:
                n=n[0]
                nazev_ind.append(n)
        # ccs=[cta[c] for c in nazev_ind ]
        #U=[cta[c].U_mag for c in nazev_ind ]
        U=cta[nazev_ind[0]].V_aver
        if len(nazev_ind)>1:
            for c in nazev_ind[1:]:
                U=np.vstack((U,cta[c].V_aver))
            fig, axis = plt.subplots(figsize=(5,5),subplot_kw=dict(projection='polar'))
            phi=[]
            for nzv in nazev_ind:
                phi.append(int(file_paths[nzv].split('_')[-1][:-4]))
            #phi=[int(c) for c in natoceni]
            phi=np.radians(phi)
            
            rad=150-cta[nazev_ind[-1]].x
            
            phis=np.array(phi[:-2])
            # rads=np.array(rad)
            Us=np.copy(U[:-2])
            for i in range(7):
                phis=np.append(phis,phi[:-2]+np.pi*(i+1)/4)
                Us=np.vstack((Us,U[:-2]))
            
            phis=np.append(phis,phi[-1]+np.pi*7/4)
            Us=np.vstack((Us,U[0]))
            rad, phi=np.meshgrid(rad,phis)
            CS=axis.contourf(phi,rad,Us)
            plt.colorbar(CS)
            #plt.rcParams.update({'font.size': 8})
            plt.savefig(file_paths[nazev_ind[-1]][:-16]+'V_aver_2D.png',dpi=300,bbox_inches='tight')
            nazev_ind=[]
            
            plt.close()


nazev_ind=[]
for v in vir:
    for ax_pozice in ax:
        if ax_pozice == '0,1D':
            ax_pozice='01D'
        elif ax_pozice == '1D':
            ax_pozice='10D'
        
        for nat in natoceni:
            txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
            n=[n for n,s in enumerate(file_paths) if txt in s]
            if len(n) > 0:
                n=n[0]
                nazev_ind.append(n)
        # ccs=[cta[c] for c in nazev_ind ]
        #U=[cta[c].U_aver for c in nazev_ind ]
        if len(nazev_ind)>1:
            U=cta[nazev_ind[0]].U_aver
            V=cta[nazev_ind[0]].U_aver
            kt=0.5*((U**2+V**2)**0.5)
            for c in nazev_ind[1:]:
                U=cta[c].U_aver
                V=cta[c].V_aver
                kt=np.vstack((kt,0.5*((U**2+V**2)**0.5)))
            fig, axis = plt.subplots(figsize=(5,5),subplot_kw=dict(projection='polar'))
            phi=[]
            for nzv in nazev_ind:
                phi.append(int(file_paths[nzv].split('_')[-1][:-4]))
            #phi=[int(c) for c in natoceni]
            phi=np.radians(phi)
            
            rad=150-cta[nazev_ind[-1]].x
            
            phis=np.array(phi[:-2])
            # rads=np.array(rad)
            kts=np.copy(kt[:-2])
            for i in range(7):
                phis=np.append(phis,phi[:-2]+np.pi*(i+1)/4)
                kts=np.vstack((kts,kt[:-2]))
            
            phis=np.append(phis,phi[-1]+np.pi*7/4)
            kts=np.vstack((kts,kt[0]))
            rad, phi=np.meshgrid(rad,phis)
            CS=axis.contourf(phi,rad,kts)
            plt.colorbar(CS)
            #plt.rcParams.update({'font.size': 8})
            plt.savefig(file_paths[nazev_ind[-1]][:-16]+'kt_2D.png',dpi=300,bbox_inches='tight')
            nazev_ind=[]
            
            plt.close()



for v in vir:
  for ax_pozice in ax:
    if ax_pozice == '0,1D':
      ax_pozice='01D'
    elif ax_pozice == '1D':
      ax_pozice='10D'
      
    for nat in natoceni:
      txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind.append(n)
    ccs=[cta[c] for c  in nazev_ind ]
    S=hw.swirl_number(ccs, 150)
    print file_paths[nazev_ind[-1]][:5]+': Swirl nuber:',S
    nazev_ind=[]



#-------Zde definuje která natočení se budou tisknout do grafů------
natoceni=['00','30', '40']
#---------------------------------------------------------


#print ax, vir, natoceni
nazev_ind=0

#složený graf U_mag přes natočení i vzdálenosti
for v in vir:
  fig, axes = plt.subplots(5, 1, figsize=(5,15),sharex=True, sharey=True)
  for xi,ax_pozice in enumerate(ax):
   
   
    xi=4-xi 
    for nat in natoceni:
      txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
        r_D=(150-cta[n].x)/300.
        U_U0=cta[n].U_mag/U0
        axes[xi].plot(r_D,U_U0,label=nat, marker='+')
        
    if xi==4:
        axes[xi].set_xlabel('r/D [-]',fontsize=12)
        axes[xi].set_ylim([0,5.5])
   
    axes[xi].legend(loc=2)
    axes[xi].grid()
    axes[xi].set_ylabel('U_mag/U0 [-]',fontsize=8)
    text='x/D='+ax_pozice
    yticks = [0,1,2,3,4,5]
    axes[xi].text(0.15, 5, text,fontsize=8,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
  fig.subplots_adjust(hspace=0)
  plt.setp([a.set_yticklabels(yticks) for a in fig.axes[:-2]], visible=True)
  plt.rcParams.update({'font.size': 8})
  print 'Writing file: '+file_paths[nazev_ind][:-21]+'U_mag.png' #('/').join([v,ax_pozice,'U_mag.png'])
  
  plt.savefig(file_paths[nazev_ind][:-21]+'U_mag.png',dpi=500,bbox_inches='tight')
  plt.close()

#složený graf U_aver přes natočení i vzdálenosti
  
for v in vir:
  fig, axes = plt.subplots(5, 1, figsize=(5, 15),sharex=True, sharey=True)
  for xi,ax_pozice in enumerate(ax):
   
   
    xi=4-xi 
    for nat in natoceni:
      txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
        r_D=(150-cta[n].x)/300.
        U_U0=cta[n].U_aver/U0
        axes[xi].plot(r_D,U_U0,label=nat, marker='+')
        
    if xi==4:
        axes[xi].set_xlabel('r/D [-]',fontsize=12)
        axes[xi].set_ylim([0,5])
    axes[xi].legend(loc=2)
    axes[xi].grid()
    axes[xi].set_ylabel('U/U0 [-]',fontsize=8)
    
    text='x/D='+ax_pozice
    yticks = [0,1,2,3,4]
    axes[xi].text(0.15, 4.45, text,fontsize=8,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
  fig.subplots_adjust(hspace=0)
  plt.setp([a.set_yticklabels(yticks) for a in fig.axes[:-2]], visible=True)
  plt.rcParams.update({'font.size': 8})
  print 'Writing file: '+file_paths[nazev_ind][:-21]+'U_aver.png' #('/').join([v,ax_pozice,'U_mag.png'])
  
  plt.savefig(file_paths[nazev_ind][:-21]+'U_aver.png',dpi=500,bbox_inches='tight')
  plt.close()

#složený graf U_aver přes natočení i vzdálenosti

for v in vir:
  fig, axes = plt.subplots(5, 1, figsize=(5, 15),sharex=True, sharey=True)
  for xi,ax_pozice in enumerate(ax):
   
   
    xi=4-xi 
    for nat in natoceni:
      txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
        r_D=(150-cta[n].x)/300.
        V_U0=cta[n].V_aver/U0
        axes[xi].plot(r_D,V_U0,label=nat, marker='+')
        
    if xi==4:
        axes[xi].set_xlabel('r/D [-]',fontsize=12)
        axes[xi].set_ylim([-1,3])
    axes[xi].legend(loc=2)
    axes[xi].grid()
    axes[xi].set_ylabel('V/U0 [-]',fontsize=8)
    
    text='x/D='+ax_pozice
    yticks = [-1,-0.5,0,0.5,1,1.5,2,2.5]
    axes[xi].text(0.13, 2.5, text,fontsize=8,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
  fig.subplots_adjust(hspace=0)
  plt.setp([a.set_yticklabels(yticks) for a in fig.axes[:-2]], visible=True)
  plt.rcParams.update({'font.size': 8})
  print 'Writing file: '+file_paths[nazev_ind][:-21]+'V_aver.png' #('/').join([v,ax_pozice,'U_mag.png'])
  
  plt.savefig(file_paths[nazev_ind][:-21]+'V_aver.png',dpi=500,bbox_inches='tight')
  plt.close()      

#složený graf Angle přes natočení i vzdálenosti

for v in vir:
  fig, axes = plt.subplots(5, 1, figsize=(5, 15),sharex=True, sharey=True)
  for xi,ax_pozice in enumerate(ax):
   
   
    xi=4-xi 
    for nat in natoceni:
      txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
        r_D=(150-cta[n].x)/300.
        Angle=cta[n].phi
        axes[xi].plot(r_D,Angle,label=nat, marker='+')
        
    if xi==4:
        axes[xi].set_xlabel('r/D [-]',fontsize=12)
        axes[xi].set_ylim([-45,45])
    axes[xi].legend(loc=3)
    axes[xi].grid()
    axes[xi].set_ylabel('Angle [deg]',fontsize=8)
    
    
    text='x/D='+ax_pozice
    yticks = [-45,-40,-20,0,20, 40]
    axes[xi].text(0.02, 35, text,fontsize=8,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
  fig.subplots_adjust(hspace=0)
  plt.setp([a.set_yticklabels(yticks) for a in fig.axes[:-2]], visible=True)
  plt.rcParams.update({'font.size': 8})
  print 'Writing file: '+file_paths[nazev_ind][:-21]+'Angle.png' #('/').join([v,ax_pozice,'U_mag.png'])
  
  plt.savefig(file_paths[nazev_ind][:-21]+'Angle.png',dpi=500,bbox_inches='tight')
  plt.close() 

#složený graf Ek přes natočení i vzdálenosti

for v in vir:
  fig, axes = plt.subplots(5, 1, figsize=(5, 15),sharex=True, sharey=True)
  for xi,ax_pozice in enumerate(ax):
   
   
    xi=4-xi 
    for nat in natoceni:
      txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
        r_D=(150-cta[n].x)/300.
        k_T=(0.5*((cta[n].U_RMS**2+cta[n].U_RMS**2)**0.5))
        axes[xi].plot(r_D,k_T,label=nat, marker='+')
        
    if xi==4:
        axes[xi].set_xlabel('r/D [-]',fontsize=12)
        axes[xi].set_ylim([0,3])
    axes[xi].legend(loc=2)
    axes[xi].grid()
    axes[xi].set_ylabel('k_T [m2/s2]',fontsize=8)
    text='x/D='+ax_pozice
    yticks = [0,0.5,1,1.5,2,2.5]
    axes[xi].text(0.13,2.5, text,fontsize=8,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
  fig.subplots_adjust(hspace=0)
  plt.setp([a.set_yticklabels(yticks) for a in fig.axes[:-2]], visible=True)
  plt.rcParams.update({'font.size': 8})
  print 'Writing file: '+file_paths[nazev_ind][:-21]+'k_T.png' #('/').join([v,ax_pozice,'U_mag.png'])
  
  plt.savefig(file_paths[nazev_ind][:-21]+'k_T.png',dpi=500,bbox_inches='tight')
  plt.close() 

#grafy jen pro natočení

#for v in vir:
  #for ax_pozice in ax:
#    if ax_pozice == '0,1D':
#      ax_pozice='01D'
#    elif ax_pozice == '1D':
#      ax_pozice='10D'
      
    #for nat in natoceni:
      #txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      #n=[n for n,s in enumerate(file_paths) if txt in s]
      #if len(n) > 0:
        #n=n[0]
       # nazev_ind=n
       # r_D=(150-cta[n].x)/300.
       # U_U0=cta[n].U_aver/U0
       # plt.plot(r_D,U_U0,label=nat, marker='+')
    #plt.xlabel('r/D [-]',fontsize=16)
    #plt.ylabel('U/U0 [-]',fontsize=16)
    #plt.legend(loc=2)
    #plt.grid()
    #plt.rcParams.update({'font.size': 12})
    #plt.savefig(file_paths[nazev_ind][:-16]+'U_aver.png',dpi=500,bbox_inches='tight')
    #print 'Writing file: '+file_paths[nazev_ind][:-16]+'U_aver.png' #('/').join([v,ax_pozice,'U_aver.png'])
    
    #plt.close()
    
    #for nat in natoceni:
    #  txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
    #  n=[n for n,s in enumerate(file_paths) if txt in s]
    #  if len(n) > 0:
    #    n=n[0]
    #    nazev_ind=n
    #    r_D=(150-cta[n].x)/300.
    #    V_U0=cta[n].V_aver/U0
    #    plt.plot(r_D,V_U0,label=nat, marker='+')
    #plt.xlabel('r/D [-]',fontsize=16)
    #plt.ylabel('V/U0 [-]',fontsize=16)
    #plt.legend(loc=2)
    #plt.grid()
    #plt.rcParams.update({'font.size': 12})
    #plt.savefig(file_paths[nazev_ind][:-16]+'V_aver.png',dpi=500,bbox_inches='tight')
    #print 'Writing file: '+file_paths[nazev_ind][:-16]+'V_aver.png' #('/').join([v,ax_pozice,'V_aver.png'])
    #plt.close()
    

    #for nat in natoceni:
    #  txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
    #  n=[n for n,s in enumerate(file_paths) if txt in s]
    #  if len(n) > 0:
    #    n=n[0]
    #    nazev_ind=n
    #    r=(150-cta[n].x)
    #    Ek=(0.5*((cta[n].U_RMS**2+cta[n].U_RMS**2)**0.5))
    #    plt.plot(r,Ek,label=nat, marker='+')
    #plt.xlabel('r [mm]',fontsize=16)
    #plt.ylabel('Ek [m2/s2]',fontsize=16)
    #plt.xlim([0,150.2])
    #plt.legend(loc=2)
    #plt.grid()
    #plt.rcParams.update({'font.size': 12})
    #plt.savefig(file_paths[nazev_ind][:-16]+'Ek.png',dpi=500,)
    #print 'Writing file: '+file_paths[nazev_ind][:-16]+'Ek.png' #('/').join([v,ax_pozice,'V_aver.png'])
    #plt.close()
    
    #for nat in natoceni:
    #  txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
    #  n=[n for n,s in enumerate(file_paths) if txt in s]
    #  if len(n) > 0:
    #    n=n[0]
    #    nazev_ind=n
    #    r=(150-cta[n].x)
    #    phi=cta[n].phi
    #    plt.plot(r,phi,label=nat,marker='+')
    #plt.xlabel('r [mm]',fontsize=16)
    #plt.ylabel('angle [deg]',fontsize=16)
    #plt.xlim([0,150.2])
    #plt.legend(loc=2)
    #plt.grid()
    #plt.rcParams.update({'font.size': 12})
    #plt.savefig(file_paths[nazev_ind][:-16]+'angle.png',dpi=500,)
    #print 'Writing file: '+file_paths[nazev_ind][:-16]+'angle.png' #('/').join([v,ax_pozice,'V_aver.png'])
    #plt.close()

    
##------------------------Tisk grafů U_mag vs x/D--------------------
natoceni=['00','30','40']
nn=[]
x=40 #Pro jaký bod (poloměr) vykreslit axiální závislost
ap=np.array([])
un=np.array([])
phi=np.array([])
Ek1=np.array([])
Ek2=np.array([])
for v in vir:
  for nat in natoceni:
  
    for ax_pozice in ax:
      if ax_pozice == '0,1D':
        txt=('D\\'+v[6:]+'_'+'01'+'_'+nat+'.txt').replace(',','')
        
      elif ax_pozice == '1D':
        txt=('D\\'+v[6:]+'_'+'10'+'_'+nat+'.txt').replace(',','')
      else:
          txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
        ui=np.where(cta[n].x==x)
        if len(ui)>0:
          un=np.append(un,cta[n].U_mag[ui])
          ap=np.append(ap,float(ax_pozice[:-1].replace(',','.')))
        
    un_i=np.argsort(ap)
    un=un[un_i]
    ap=ap[un_i]
    U_U0=un/U0
    plt.plot(ap,U_U0,label=nat,marker='+')
    ap=np.array([]) #TODO nesmaže to předchozí graf?
    un=np.array([])
    
  plt.xlabel('x/D [-]',fontsize=16)
  plt.ylabel('U_mag/U0 [-]',fontsize=16)
  plt.ylim([0,4])
  plt.legend(loc=0)
  text='r ='+str(150-x)
  plt.text(0.2,3.5, text, fontsize=8,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
  plt.grid()
  plt.rcParams.update({'font.size': 12})
  plt.savefig(file_paths[nazev_ind][:26]+'U_mag_ax.png',dpi=500,bbox_inches='tight')
  print 'Writing file: '+file_paths[nazev_ind][:26]+'U_mag_ax.png' #('/').join([v,ax_pozice,'U_mag_ax.png'])
  plt.close()   



for v in vir:
  for nat in natoceni:
    
    for ax_pozice in ax:
      if ax_pozice == '0,1D':
        txt=('D\\'+v[6:]+'_'+'01'+'_'+nat+'.txt').replace(',','')
        
      elif ax_pozice == '1D':
        txt=('D\\'+v[6:]+'_'+'10'+'_'+nat+'.txt').replace(',','')
      else:
        txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
        ui=np.where(cta[n].x==x)
        if len(ui)>0:
          phi=np.append(phi,cta[n].phi[ui])
          ap=np.append(ap,float(ax_pozice[:-1].replace(',','.')))
        
    phi_i=np.argsort(ap)
    phi=phi[phi_i]
    ap=ap[phi_i]
    plt.plot(ap,phi,label=nat,marker='+')
    ap=np.array([]) #TODO nesmaže to předchozí graf?
    phi=np.array([])
  plt.xlabel('x/D [-]',fontsize=16)
  plt.ylabel('angle [deg]',fontsize=16)
  plt.ylim([-45,45])
  plt.legend(loc=0)
  text='r ='+str(150-x)
  plt.text(0.2,35, text, fontsize=8,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
  plt.grid()
  plt.rcParams.update({'font.size': 12})
  plt.savefig(file_paths[nazev_ind][:26]+'angle_x.png',dpi=500,bbox_inches='tight')
  print 'Writing file: '+file_paths[nazev_ind][:26]+'angle_x.png' #('/').join([v,ax_pozice,'U_mag_ax.png'])
  plt.close() 

for v in vir:
  for nat in natoceni:

    for ax_pozice in ax:
      if ax_pozice == '0,1D':
        txt=('D\\'+v[6:]+'_'+'01'+'_'+nat+'.txt').replace(',','')
        
      elif ax_pozice == '1D':
        txt=('D\\'+v[6:]+'_'+'10'+'_'+nat+'.txt').replace(',','')
      else:
          txt=('D\\'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
        ui=np.where(cta[n].x==x)
        if len(ui)>0:
          Ek1=np.append(Ek1,cta[n].U_RMS[ui])
          Ek2=np.append(Ek2,cta[n].V_RMS[ui])
          ap=np.append(ap,float(ax_pozice[:-1].replace(',','.')))
    ap_i=np.argsort(ap)
    Ek1=Ek1[ap_i]
    Ek2=Ek2[ap_i]
    ap=ap[ap_i]
    Ek=(0.5*((Ek1**2+Ek2**2)**0.5))
    plt.plot(ap,Ek,label=nat,marker='+')
    ap=np.array([]) #TODO nesmaže to předchozí graf?
    Ek1=np.array([])
    Ek2=np.array([])
  plt.xlabel('x/D [-]',fontsize=16)
  plt.ylabel('k_T [m2/s2]',fontsize=16)
  plt.ylim([0,3])
  plt.legend(loc=0)
  text='r ='+str(150-x)
  plt.text(0.2,2.7, text, fontsize=8,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
  plt.grid()
  plt.rcParams.update({'font.size': 12})
  plt.savefig(file_paths[nazev_ind][:26]+'k_T-x.png',dpi=500,bbox_inches='tight')
  print 'Writing file: '+file_paths[nazev_ind][:26]+'k_T-x.png' #('/').join([v,ax_pozice,'U_mag_ax.png'])
  plt.close()      
## ----------------------------------------------------------------
    
    
    
def get_sizes(filename):
  '''
  From the file name extracts dimensions
  '''
  if '.txt' in filename:
    fnm=os.path.basename(filename)
    rozmery=fnm[:-4].split('_')
    return rozmery
  elif 'viric_' in filename:
    dirnm=os.path.dirname(filename).split('\\')
    m=[n for n,s in enumerate(dirnm) if 'viric_' in s][0]
    vir=dirnm[m][6:-4]
    lopatky=dirnm[m][-2:-1]
    ax=dirnm[m+1][:-2]
    if ax == '0,1':
      ax='01'
    elif ax == '1':
      ax ='10'
    return [vir,lopatky,ax,'']
    
    
    
    
    
    
    
    
  #    print lst

#print 'D/'+v[6:]+'_'+a[-2]+'_'
#n=[n for n,s in enumerate(file_paths) if '01_45.txt' in s]
#print n
#lst=[file_paths[i] for i in n]
#print lst



#for v in vir:
#  virice=[s for s in file_paths if v in s]
#  for a in ax:
#    ax_pozice=[s for s in virice if '/'+a+'/' in s]
#    for n in natoceni:
#      nn=[s for s in ax_pozice if n+'.txt' in s]
#      if len(nn) > 0:
#        nazvy.append(nn[0])
#      st=[v,'/'+a+'/',n+'.txt']
#    print nazvy
#    nazvy=[]
#    print
    
    
    
    
    
    
#print nazvy
      
#      print st
#      if any(x in str for x in file_paths):
#        soubor.append(file_paths)
#for fp in file_paths:
#  if all(x in fp for x in st):
#    print fp
    
#if any("45.txt" in s for s in file_paths):
#  print "Je tam"
#  
#matching = [s for s in file_paths if vir[0] in s]
#print matching

#matching = [fp for fp in file_paths if all(x in fp for x in st)]
#print st
#print matching
#
#matching = [fp for fp in file_paths if "280_55/" in fp]
##print matching
#print os.path.dirname(matching[0])







#file_paths=file_paths[0:1]
#print file_paths
#
#for file_in in file_paths:
#  print file_in
#  rozmery.append(os.path.basename(file_in)[:-4].split('_'))
#  
#  cta.append(hw.HotWireMeasurement(file_in)) #Načte data pro daný soubor file_in, uloží je do třídy HotWireMeasurement a třídu uloží do pole cta
#  print cta[-1] #Vytiskne informace o daném měření
#
#for c,r in zip(cta,rozmery):
#  
#  D_v,uhel_lopatek,dist,natoceni=r
#  
#  plt.grid()
##  plt.axis([-0.05,0.5,1,4.5]) #0.35,1.8
#  plt.plot( c.x ,c.U_aver,label=natoceni) #prumer/U0
#  r=1-c.x
#  plt.xlabel('r [mm]')
#  plt.ylabel('U [m/s]')
#  plt.legend()
#  
##  plt.rcParams.update({'font.size': 22})
#  plt.savefig(file_in[:-4]+'_U.png',bbox_inches='tight')
#  plt.close()
#
#for c,r in zip(cta,rozmery):
#  
#  D_v,uhel_lopatek,dist,natoceni=r
#  
#  plt.ylabel('V [m/s]')
#  plt.plot( c.x ,c.V_aver) #prumer/U0  
##  plt.rcParams.update({'font.size': 22})
#  plt.savefig(file_in[:-4]+'_V.png',bbox_inches='tight')
#  print "Ukládám obrázek:",file_in[:-4]+".png"
#  plt.close()
