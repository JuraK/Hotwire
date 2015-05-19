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


#######################x
U0=3.5 #[m/s] Uniform velocity 
#######################x

prumery=[]
pozice_all=[]
cta=[]
rozmery=[]

def get_sizes(filename):
  if '.txt' in filename:
    fnm=os.path.basename(filename)
    rozmery=fnm[:-4].split('_')
    return rozmery
  elif 'viric_' in filename:
    dirnm=os.path.dirname(filename).split('/')
    m=[n for n,s in enumerate(dirnm) if 'viric_' in s][0]
    vir=dirnm[m][6:-4]
    lopatky=dirnm[m][-2:-1]
    ax=dirnm[m+1][:-2]
    if ax == '0,1':
      ax='01'
    elif ax == '1':
      ax ='10'
    return [vir,lopatky,ax,'']

#-------Gets all the files in all the subdirectories---------
directory='./Viric_data'
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
      ax.append(root.split('/')[-1]) #Get all the axial distances
      vir.append(root.split('/')[-2]) #Get all the swirler types
      if '_' in filename:
        natoceni.append(filename.split('_')[-1][:-4]) # Get all the rotation angles
#----------------------------------------------------------
        
##---Gets all the data from all files------
#for f in file_paths:
#  #Načte data pro daný soubor f, uloží je do třídy HotWireMeasurement a třídu uloží do pole cta
#  cta.append(hw.HotWireMeasurement(f)) 
##-----------------------------------------
  
#Eliminates repeates
ax=np.unique(ax)
vir=np.unique(vir)
natoceni=np.unique(natoceni)
print natoceni

#Zde definuje které natočení mě budou zajímat
natoceni=['05','15','25', '40']
print ax, vir, natoceni

nazev_ind=0
for v in vir:
  for ax_pozice in ax:
    if ax_pozice == '0,1D':
      ax_pozice='01D'
    elif ax_pozice == '1D':
      ax_pozice='10D'
      
    for nat in natoceni:
      txt=('D/'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
#        r_D=(150-cta[n].x)/300.
#        U_U0=cta[n].U_mag/U0
#        plt.plot(r_D,U_U0,label=nat)
#    plt.xlabel('r/D [-]')
#    plt.ylabel('U/U0 [-]')
#    plt.legend()
#    plt.savefig(file_paths[nazev_ind][:-16]+('/').join([v,ax_pozice,'U_mag.eps'])
    print 'Writing file: '+file_paths[nazev_ind][:-16]+'U_mag.eps' #('/').join([v,ax_pozice,'U_mag.eps'])
#    plt.close()
        
    for nat in natoceni:
      txt=('D/'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
#        r_D=(150-cta[n].x)/300.
#        U_U0=cta[n].U_aver/U0
#        plt.plot(r_D,U_U0,label=nat)
#    plt.xlabel('r/D [-]')
#    plt.ylabel('U/U0 [-]')
#    plt.legend()
#    plt.savefig(file_paths[nazev_ind][:-16]+('/').join([v,ax_pozice,'U_aver.eps'])
    print 'Writing file: '+file_paths[nazev_ind][:-16]+'U_aver.eps' #('/').join([v,ax_pozice,'U_aver.eps'])
#    plt.close()
    
    for nat in natoceni:
      txt=('D/'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
#        r_D=(150-cta[n].x)/300.
#        V_U0=cta[n].V_aver/U0
#        plt.plot(r_D,V_U0,label=nat)
#    plt.xlabel('r/D [-]')
#    plt.ylabel('V/U0 [-]')
#    plt.legend()
#    plt.savefig(file_paths[nazev_ind][:-16]+('/').join([v,ax_pozice,'V_aver.eps'])
    print 'Writing file: '+file_paths[nazev_ind][:-16]+'V_aver.eps' #('/').join([v,ax_pozice,'V_aver.eps'])
#    plt.close()
    
    
        
#    print file_paths[nazev_ind][:-16]+'U_mag.eps'
    
##------------------------Tisk grafů U_mag vs x/D--------------------
natoceni=['00','20']
nn=[]
x=10 #Pro jaký bod (poloměr) vykreslit axiální závislost
ap=[]
un=np.array([])
for v in vir:
  for nat in natoceni:
  
    for ax_pozice in ax:
      if ax_pozice == '0,1D':
        ax_pozice='01D'
      elif ax_pozice == '1D':
        ax_pozice='10D'

      txt=('D/'+v[6:]+'_'+ax_pozice[:-1]+'_'+nat+'.txt').replace(',','')
      n=[n for n,s in enumerate(file_paths) if txt in s]
      if len(n) > 0:
        n=n[0]
        nazev_ind=n
#        ui=np.where(cta[n].x==x)
#        if ui.shape()[0]>0:
#          un=np.append(un,cta[n].U_mag[ui])
        ap.append(float(ax_pozice[:-1].replace(',','.')))
        
    
    U_U0=un/U0
#    plt.plot(ap,U_U0,label=nat)
    ap=[] #TODO nesmaže to předchozí graf?
    un=np.array([])
#  plt.xlabel('x/D [-]')
#  plt.ylabel('U/U0 [-]')
#  plt.legend()
#  plt.savefig(file_paths[nazev_ind][:26]+('/').join([v,ax_pozice,'U_mag_ax.eps'])
  print 'Writing file: '+file_paths[nazev_ind][:26]+'U_mag_ax.eps' #('/').join([v,ax_pozice,'U_mag_ax.eps'])
#  plt.close()      
## ----------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    
    
    
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
#  plt.savefig(file_in[:-4]+'_U.eps',bbox_inches='tight')
#  plt.close()
#
#for c,r in zip(cta,rozmery):
#  
#  D_v,uhel_lopatek,dist,natoceni=r
#  
#  plt.ylabel('V [m/s]')
#  plt.plot( c.x ,c.V_aver) #prumer/U0  
##  plt.rcParams.update({'font.size': 22})
#  plt.savefig(file_in[:-4]+'_V.eps',bbox_inches='tight')
#  print "Ukládám obrázek:",file_in[:-4]+".eps"
#  plt.close()
