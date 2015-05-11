# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 15:47:01 2015

@author: jura
"""

import matplotlib.pyplot as plt
import HotWire_calc as hw

#---Temperature correction---

#######################x
T1 = False #Temperature during experimet, put False to turn it off
#######################x

#Délka sondy + držáku
Lsondy=255.
Pozice_korekce=Lsondy-250-3

prumery=[]
pozice_all=[]
cta=[]

#directory="./2015_01_28/"
directory='Viric_data/viric_240_25/0,1D/'

#files=['Tryska_Alkion_300mm.txt','Tryska_Alkion_500mm.txt','Tryska_Tuma_v00-300mm.txt','Tryska_Tuma_v00-500mm.txt']
#files=['GlobalExport_vir_01.txt']
files=['240_25_01_00.txt','240_25_01_05.txt','240_25_01_10.txt']
rozmery=[]

for file_in in files:
  rozmery.append(file_in[:-4].split('_'))
  
  cta.append(hw.HotWireMeasurement(directory+file_in)) #./2015_01_28/GlobalExport_vir_01.txt Načte data pro daný soubor file_in, uloží je do třídy HotWireMeasurement a třídu uloží do pole cta
  #cta[-1].calibration([23.627550,-47.563171,36.200703,-13.087972,2.099148,0])
  #cta[-1].prumery=np.array([np.average(i.data) for i in cta[-1].data_points]) #Průměrné hodnoty rychlostí
  #cta[-1].RMS=np.array([np.std(i.data) for i in cta[-1].data_points]) 
  
#  cta[0].directional_calibration()
#  print cta[-1] #Vytiskne informace o daném měření

for c,r in zip(cta,rozmery):
  
  D_v,uhel_lopatek,dist,natoceni=r
  
  plt.grid()
#  plt.axis([-0.05,0.5,1,4.5]) #0.35,1.8
  plt.plot( c.x ,c.U_aver,label=natoceni) #prumer/U0
  

  plt.xlabel('r [mm]')
  plt.ylabel('U [m/s]')
  plt.legend()
  
#  plt.rcParams.update({'font.size': 22})
  plt.savefig(directory+file_in[:-4]+'_U.eps',bbox_inches='tight')
plt.show()
plt.close()

for c,r in zip(cta,rozmery):
  
  D_v,uhel_lopatek,dist,natoceni=r
  
  plt.ylabel('V [m/s]')
  plt.plot( c.x ,c.V_aver) #prumer/U0  
#  plt.rcParams.update({'font.size': 22})
  plt.savefig(directory+file_in[:-4]+'_V.eps',bbox_inches='tight')
  
  #print "Ukládám obrázek:",file_in[:-4]+".eps"
  ##plt.show()
  #plt.close()

  ##plt.plot(data_time[0], data[0])
  ##plt.show()

#print 'Srovnání průměrů:'
#for p in range(cta[0].prumery.shape[0]):
  #print
  #for i in cta:  
    #print i.prumery[p],



#plt.grid()
#plt.plot( cta[0].x ,cta[0].U_aver,'x-',label=files[0][:-4]) #prumer/U0
#plt.plot( cta[0].x ,cta[0].V_aver,'x-',label=files[0][:-4]+' 2') #prumer/U0
##plt.plot( cta[1].x ,cta[1].prumery,'x-',label=files[1]) #prumer/U0
##plt.plot( cta[2].x ,cta[2].prumery,'x-',label=files[2]) #prumer/U0
#plt.xlabel('r [mm]')
#plt.ylabel('U [m/s]')
#plt.legend()
#plt.rcParams.update({'font.size': 22})
###plt.savefig(directory+'srovnani_300.eps',bbox_inches='tight')
#plt.show()

#plt.close()
#plt.grid()
##plt.plot( cta[1].x ,cta[1].prumery,'x-') #prumer/U0
##plt.plot( cta[3].x ,cta[3].prumery,'x-') #prumer/U0
##plt.xlabel('r [mm]')
##plt.ylabel('U [m/s]')
##plt.rcParams.update({'font.size': 22})
###plt.savefig('srovnani_500.eps',bbox_inches='tight')
##plt.show()
##plt.close()