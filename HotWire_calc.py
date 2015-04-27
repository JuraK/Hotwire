#!/usr/bin/python
# -*- coding: UTF-8 -*-

from __future__ import division
from math import *
import numpy as np
import matplotlib.pyplot as plt

class Data_Point:
  '''
  Data structure for storage all the data and information from one measured point.
  data_time - np.array of time records corresponding to the measured data
  data - np.array of all the raw data from measurement at the position
  ac_date - acquisition date
  ac_time - acquisition time
  x,y,z - coordinates of measurement point
  position_no - number (label) of this position
  data2 - np.array of data from second wire
  '''
  def __init__(self,data_time, data, ac_date, ac_time, x, y, z, position_no, data2=False):
    self.data_time=data_time
    self.data=data
    self.data2=data2
    self.ac_date=ac_date
    self.ac_time=ac_time
    self.x=x
    self.y=y
    self.z=z
    self.position_no=position_no


class CTA_Mereni:
  '''
  This class is responsible for all the data processing
  
  soubor='./data/soubor.txt'
  kalibrace=[0, 1, 2, 3, 4, 5] # nebo =[ [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5] ] pro dvoudrátkovou sondu
  '''
  def __init__(self, soubor, calibr_coeff=-1, teplota_kalibrace=20, teplota_experiment=7, teplota_zhaveni=260):
    self.data_points=[] #Pole tříd jednotlivých měření (data + informace o daném měření - čas, datum, pozice, číslo pozice)
    self.sensors_calibration=[] #Kalibrační konstanty pro každou drátek sondy
    self.sensors_no=-1 #Počet senzorů
    self.corr_T=1
    self.probes_no=-1 #Počet sond
    self.positions=-1 #Počet měřicích pozic
    self.temperature_probe=False
    self.directional_calib_made=False
    
    #------Pro směrovou kalibraci------
    self.U1_aver=np.array([]) #Ve směru drátku
    self.U2_aver=np.array([]) #Ve směru drátku
    self.U1_RMS=np.array([]) #Ve směru drátku
    self.U2_RMS=np.array([]) #Ve směru drátku
    self.U_aver=np.array([]) #Ve směru sondy
    self.V_aver=np.array([]) #Ve směru sondy
    self.U_RMS=np.array([]) #Ve směru sondy
    self.V_RMS=np.array([]) #Ve směru sondy
    
    self.temper_corr(teplota_experiment, teplota_zhaveni, teplota_kalibrace) #Vypočte teplotní korekci
    self.read_data(soubor) #Načte data
    
    print "calibr_coeff =", calibr_coeff
    if calibr_coeff != -1:
      print 'Adjusting data with given calibration coefficients.' #'Provádím kalibraci dle externě zadaných koeficientů.'
      if (len(calibr_coeff)==2 and self.sensors_no == 2):
        self.sensors_calibration=calibr_coeff
      elif len(calibr_coeff) == 6:
        self.sensors_calibration=[calibr_coeff]
      
    self.kalibrace(self.sensors_calibration) #Provede přepočet z napětí na rychlosti dle kalibračních koeficientů
    
    self.x=np.array([i.x for i in self.data_points]) #Pozice měření (uvažuje se 1D posuv)
    self.prumery=np.array([np.average(i.data) for i in self.data_points]) #Průměrné hodnoty rychlostí
    self.RMS=np.array([np.std(i.data) for i in self.data_points]) 
    if self.sensors_no >1:
      self.prumery2=np.array([np.average(i.data2) for i in self.data_points]) #Průměrné hodnoty rychlostí
      self.RMS2=np.array([np.std(i.data2) for i in self.data_points]) 
      
  
  def temper_corr(self, T1, Tw, T0):
    '''
    T1 - [°C] temperature during experiment
    Tw - [°C] wire temperature
    T0 - [°C] reference temperature (during calibration)
    
    Ještě jiný postup je popsán v 
    [1]   HULTMARK, M. and SMITS, A.J. Temperature corrections for constant temperature and constant current hot-wire anemometers. In: Measurement Science and Technology. 1 October 2010, Vol. 21, no. 10, pp. 105404. DOI 10.1088/0957-0233/21/10/105404. .
    
    U/nu=f(E^2/(k*dT))
    dT=Tw-Ta 
    '''
    if not((-20 < T1 < 300) and (0<Tw<800) and (-20<T0<300)):
      print "Error!: Wrong temperature for temperature correction"
      self.corr_T = -1
      return -1    
    
    if T1 is not False:
      #m = 2e-5+5e-8*T1+3e-11*T1**2 # Temperature loading factor - from physical library [StreamWare Pro, Installation and User Guide], temperature input in [°C], page no. 266 - WRONG!!
      m=0.2
      T1+=273.15
      Tw+=273.15
      T0+=273.15
      print "m =",m
      #m=0.1
      if T0 > T1:
        self.corr_T = ((Tw-T0)/(Tw-T1))**(0.5*(1-m))
      else:
        self.corr_T = ((Tw-T0)/(Tw-T1))**(0.5*(1+m))
    else:
      self.corr_T = 1 
      
  def __repr__(self):
    print "Coefficient for temperature adjustment:", self.corr_T,"->resulting change by",100-self.corr_T*100,"%" #Koeficient pro teplotní kalibraci
    print "Calibration coefficients C0..C5:",self.sensors_calibration# C0,C1,C2,C3,C4,C5
    print 'Number of measurement points x:',len(self.data_points) #data.shape[0]
    print 'Number of samples at one measurement point:',self.data_points[0].data.shape[0] #data.shape[1]
    print
    print "Number of sensors (wires):", self.sensors_no
    if self.sensors_no == 1:
      print '{0:>5s} {1:>9s} {2:>7s} {3:>7s}'.format('x', 'Average', 'RMS', 'I')
      print '{0:>5s} {1:>7s} {2:>7s} {3:>7s}'.format('[mm]', '[m/s]', '[m/s]', '[%]')
      for i in range(len(self.data_points)):
        print '{0:>5g} {1:>7.3f} {2:>7.3f} {3:>7.3f}'.format(self.x[i], self.prumery[i], self.RMS[i], self.RMS[i]/self.prumery[i]*100)
    elif self.sensors_no > 1:
      print '{0:>5s} {1:>9s} {2:>7s} {3:>7s} {4:>9s} {5:>7s} {6:>7s}'.format('x', 'Average', 'RMS', 'I', 'AVerage2', 'RMS2', 'I2')
      print '{0:>5s} {1:>7s} {2:>7s} {3:>7s} {4:>7s} {6:>7s} {6:>7s}'.format('[mm]', '[m/s]', '[m/s]', '[%]', '[m/s]', '[m/s]', '[%]')
      for i in range(len(self.data_points)):
        print '{0:>5g} {1:>7.3f} {2:>7.3f} {3:>7.3f} {4:>7.3f} {5:>7.3f} {6:>7.3f}'.format(self.x[i], self.prumery[i], self.RMS[i], self.RMS[i]/self.prumery[i]*100, self.prumery2[i], self.RMS2[i], self.RMS2[i]/self.prumery2[i]*100)
      if self.directional_calib_made:
        print
        print '{0:>5s} {1:>9s} {2:>7s} {3:>7s} {4:>9s} {5:>7s} {6:>7s}'.format('x', 'U', 'U_RMS', 'U_I', 'V', 'V_RMS', 'V_I')
        print '{0:>5s} {1:>7s} {2:>7s} {3:>7s} {4:>7s} {6:>7s} {6:>7s}'.format('[mm]', '[m/s]', '[m/s]', '[%]', '[m/s]', '[m/s]', '[%]')
        for i in range(len(self.data_points)):
          print '{0:>5g} {1:>7.3f} {2:>7.3f} {3:>7.3f} {4:>7.3f} {5:>7.3f} {6:>7.3f}'.format(self.x[i], self.U_aver[i], self.U_RMS[i], self.U_RMS[i]/self.U_aver[i]*100, self.V_aver[i], self.V_RMS[i], self.V_RMS[i]/self.V_aver[i]*100)
    
    return ''
  
  
  def read_data(self,soubor):
    hodnoty=np.array([])
    hodnoty2=np.array([])
    hodnoty_time=np.array([])
    sw=-1
    sens_order=-1
    
    print
    print '----File: ',file_in,'----'

    with open(soubor,'r') as fh:
      for line in fh:
        
        if line[0:14]=='[FILE HEADER]:':
          sw=0
        elif line[0:14]=='[USER HEADER]:':
          sw=1
        elif line[0:15]=='[PROBE HEADER]:':
          sw=2
          temperature_probe=False
        elif line[0:14]=='[DATA HEADER]:':
          sw=3
          if hodnoty.shape[0]>0:
            self.data_points.append( Data_Point(hodnoty_time, hodnoty, ac_date, ac_time, x, y, z, position_no, hodnoty2) )
            #print "Hodnoty:", hodnoty[0]
            hodnoty=np.array([])
            hodnoty2=np.array([])
            hodnoty_time=np.array([])
        elif '[DATA BLOCK' in line:
          sw=4
        
        if sw==0:
          if 'Raw Data Event Date:' in line:
            ac_date=line.split()[-1]
          elif 'Number of probes selected:' in line:
            self.probes_no=line.split()[-1]
          elif 'Number of positions:' in line:
            self.positions=line.split()[-1]
        
        if sw == 2:
          if 'teplot' in line:
            temperature_probe=True
            self.temperature_probe=True
          
          if not temperature_probe:
            
            if self.sensors_no==2 and sens_order==1:
              #Načte kalibrační koeficienty C0..C5 pro drátek 2
              if line[0:3] == "C0:" and sens_order == 1:
                C0=float(line.split()[1])
              elif line[0:3] == "C1:" and sens_order == 1:
                C1=float(line.split()[1])
              elif line[0:3] == "C2:" and sens_order == 1:
                C2=float(line.split()[1])
              elif line[0:3] == "C3:" and sens_order == 1:
                C3=float(line.split()[1])
              elif line[0:3] == "C4:" and sens_order == 1:
                C4=float(line.split()[1])
              elif line[0:3] == "C5:" and sens_order == 1:
                C5=float(line.split()[1])
                self.sensors_calibration.append([C0,C1,C2,C3,C4,C5])
                #print "Kalibrační koeficienty2: ", self.sensors_calibration
            
            if 'No. of Sensors' in line:
              self.sensors_no=float(line.split()[-1])
              sens_order=0
            #Načte kalibrační koeficienty C0..C5 pro drátek 1
            elif line[0:3] == "C0:" and sens_order == 0:
              C0=float(line.split()[1])
            elif line[0:3] == "C1:" and sens_order == 0:
              C1=float(line.split()[1])
            elif line[0:3] == "C2:" and sens_order == 0:
              C2=float(line.split()[1])
            elif line[0:3] == "C3:" and sens_order == 0:
              C3=float(line.split()[1])
            elif line[0:3] == "C4:" and sens_order == 0:
              C4=float(line.split()[1])
            elif line[0:3] == "C5:" and sens_order == 0:
              C5=float(line.split()[1])
              sens_order=1 #Můžou se začít načítat data o kalibraci druhého drátku
              self.sensors_calibration.append([C0,C1,C2,C3,C4,C5])
              #print "Kalibrační koeficienty: ", self.sensors_calibration
            
            
            
        if sw==3: #DATA HEADER
          #Načte informaci o poloze sondy v aktuálním měřicím bodě
          if line[0:15] == "(mm,mm,mm,deg):":
            x,y,z,deg = [float(i) for i in line.split()[1:]]
            #print "x,y,z:", x,y,z
          #Načte obecné informace o této sadě měření (ID, čas, datum, ...)
          elif 'Acquisition time:' in line:
            ac_time = line.split()[-1]
          
          elif 'Acquisition date:' in line:
            if not '.' in line.split()[-1]:
              print 'Není tam čas !!!!!'
              ac_date = line.split()[-1]
                        
          elif 'Position no.:' in line:
            position_no= int(line.split()[-1])
            
        if sw == 4: #Jsme v bloku naměřených dat
          if len(line.strip())>0 and line.strip()[0].isdigit(): # testuje zda první pozice na řádku je číslo
            E = float(line.strip().split()[1]) * self.corr_T # Včetně teplotní korekce
            if self.sensors_no > 1:
              E2 = float(line.strip().split()[2]) * self.corr_T # Včetně teplotní korekce
            #C0,C1,C2,C3,C4,C5=self.sensors_calibration[0]
            #U=C0+C1*E+C2*E**2+C3*E**3+C4*E**4+C5*E**5 # Aplikována kalibrace
            hodnoty=np.append(hodnoty, E)
            hodnoty2=np.append(hodnoty2, E2)
            hodnoty_time=np.append(hodnoty_time, float(line.strip().split()[0]))
            #print hodnoty[-1]
          
          #elif ('[DATA HEADER' in line) and len(hodnoty)>0:
            #self.data_points.append( Data_Point(hodnoty_time, hodnoty, ac_date, ac_time, x, y, z, position_no, hodnoty2) )
            #print "Hodnoty:", hodnoty[0]
            #hodnoty=np.array([])
            #hodnoty2=np.array([])
            #hodnoty_time=np.array([])
                      
    self.data_points.append( Data_Point(hodnoty_time, hodnoty, ac_date, ac_time, x, y, z, position_no, hodnoty2) )
    #print "Hodnoty:", hodnoty[0]
            
  def kalibrace(self,coeff):
    '''
    Provede přepočet z napětí na rychlosti dle kalibračních koeficientů    
    '''
    if len(coeff) == 2 and self.sensors_no ==2: #Dvoudrátková sonda
      
      C0,C1,C2,C3,C4,C5=coeff[0]
      C0b,C1b,C2b,C3b,C4b,C5b=coeff[1]
      for dp in self.data_points:
        dp.data=C0+C1*dp.data+C2*dp.data**2+C3*dp.data**3+C4*dp.data**4+C5*dp.data**5 # Aplikována kalibrace
        dp.data2=C0b+C1b*dp.data2+C2b*dp.data2**2+C3b*dp.data2**3+C4b*dp.data2**4+C5b*dp.data2**5 # Aplikována kalibrace
    elif len(coeff) == 1: #Jednodrátková sonda
      C0,C1,C2,C3,C4,C5=coeff[0]
      for dp in self.data_points:
        dp.data=C0+C1*dp.data+C2*dp.data**2+C3*dp.data**3+C4*dp.data**4+C5*dp.data**5 # Aplikována kalibrace
    else:
      print 'Error!: Wrong input of calibration coefficients!' #'Chyba: Špatně zadané kalibrační koeficienty!'
      
  def smerova_kalibrace(self, k1, k2):
    '''
    Směrová kalibrace dle zadaných koeficientů - musí se volat extra. Není v defaultním běhu!
    '''
    for dp in self.data_points:
      U1=((2**0.5)/2)*(np.abs((1+k2**2)*dp.data2**2-k2**2*dp.data**2)**0.5)
      U2=((2**0.5)/2)*(np.abs((1+k1**2)*dp.data**2-k1**2*dp.data2**2)**0.5)

      U=((2**0.5)/2)*U1+((2**0.5)/2)*U2    #rychlost ve směru drátku 1
      V=((2**0.5)/2)*U1-((2**0.5)/2)*U2     #rychlost ve směru drátku 2
      
      self.U1_aver=np.append(self.U1_aver,np.mean(U1))
      self.U2_aver=np.append(self.U2_aver,np.mean(U2))
      
      self.U1_RMS=np.append(self.U1_RMS,np.std(U1))
      self.U2_RMS=np.append(self.U2_RMS,np.std(U2))
      
      self.U_aver=np.append(self.U_aver,np.mean(U))
      self.V_aver=np.append(self.V_aver,np.mean(V))
      
      self.U_RMS=np.append(self.U_RMS,np.std(U))
      self.V_RMS=np.append(self.V_RMS,np.std(V))
    
    self.directional_calib_made=True
      
#######################x
#itoe = 3276.8166089965 # Převod integer na napětí - vždy zkontrolovat !!!!!!!
#######################x


#---Teplotní korekce---

#######################x
T1 = False #Teplota během experimentu, zadej False pro vypnutí korekce
#######################x

#Délka sondy + držáku
Lsondy=255.
Pozice_korekce=Lsondy-250-3

#if T1 is not False:
  #Tw = 215. #Teplota žhaveného drátku
  #T0 = 25. #Referenční teplota (při kalibraci)
  #m = 2e-5+5e-8*T1+3e-11*T1**2 # Temperature loading factor - from physical library
  #if T0 > T1:
    #corr_T = ((Tw-T0)/(Tw-T1))**(0.5*(1-m))
  #else:
    #corr_T = ((Tw-T0)/(Tw-T1))**(0.5*(1+m))
#else:
  #corr_T = 1

#a = np.genfromtxt('GlobalExport1.txt',skip_header=1)

prumery=[]
pozice_all=[]
cta=[]

#directory="./2015_01_28/"
directory="./2015_02_12_2wire/"

#files=['Tryska_Alkion_300mm.txt','Tryska_Alkion_500mm.txt','Tryska_Tuma_v00-300mm.txt','Tryska_Tuma_v00-500mm.txt']
#files=['GlobalExport_vir_01.txt']
files=['01_006_short.txt']

for file_in in files:
  
  #cta.append(CTA_Mereni(directory+file_in, [23.627550,-47.563171,36.200703,-13.087972,2.099148,0], 19.5, 7, 260 )) #./2015_01_28/GlobalExport_vir_01.txt Načte data pro daný soubor file_in, uloží je do třídy CTA_Mereni a třídu uloží do pole cta
  cta.append(CTA_Mereni(directory+file_in,[[23.627550,-47.563171,36.200703,-13.087972,2.099148,0],[23.627550,-47.563171,36.200703,-13.087972,2.099148,0]]))
  #cta[-1].kalibrace([23.627550,-47.563171,36.200703,-13.087972,2.099148,0])
  #cta[-1].prumery=np.array([np.average(i.data) for i in cta[-1].data_points]) #Průměrné hodnoty rychlostí
  #cta[-1].RMS=np.array([np.std(i.data) for i in cta[-1].data_points]) 
  
  cta[0].smerova_kalibrace(12.18**0.5, 15.6**0.5)
  print cta[-1] #Vytiskne informace o daném měření
  
  ##print 'U0 =', U0
  #plt.grid()
  ##plt.axis([-0.05,0.5,1,4.5]) #0.35,1.8
  ##plt.plot( (150-(pozice+Pozice_korekce))/300 ,prumer,'o-') #prumer/U0
  #plt.plot( cta[-1].x ,cta[-1].prumery,'o-') #prumer/U0

  #plt.xlabel('r [mm]')
  #plt.ylabel('U [m/s]')
  
  #plt.rcParams.update({'font.size': 22})
  ##plt.savefig(file_in[:-4]+'.eps',bbox_inches='tight')
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



plt.grid()
plt.plot( cta[0].x ,cta[0].U_aver,'x-',label=files[0][:-4]) #prumer/U0
plt.plot( cta[0].x ,cta[0].V_aver,'x-',label=files[0][:-4]+' 2') #prumer/U0
#plt.plot( cta[1].x ,cta[1].prumery,'x-',label=files[1]) #prumer/U0
#plt.plot( cta[2].x ,cta[2].prumery,'x-',label=files[2]) #prumer/U0
plt.xlabel('r [mm]')
plt.ylabel('U [m/s]')
plt.legend()
plt.rcParams.update({'font.size': 22})
##plt.savefig(directory+'srovnani_300.eps',bbox_inches='tight')
plt.show()

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