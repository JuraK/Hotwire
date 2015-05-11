#!/usr/bin/python
# -*- coding: UTF-8 -*-

from __future__ import division
from math import *
import numpy as np
from scipy import signal

class DataPoint:
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


class HotWireMeasurement:
  '''
  This class is responsible for all the data processing
  
  in_file='./data/in_file.txt'
  calibration=[0, 1, 2, 3, 4, 5] # nebo =[ [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5] ] pro dvoudrátkovou sondu
  '''
  def __init__(self, in_file, calibr_coeff=-1, calibration_temper=20, experiment_temper=7, wire_temper=260):
    self.DataPoints=[] #List of measurement points (data + info about the measurement point -time, date, position, no of position)
    self.sensors_calibration=[] #Calibration coeff for each wire
    self.sensors_no=-1 #Number of wires
    self.corr_T=1
    self.probes_no=-1 #Number of probes
    self.positions=-1 #Number of measured positions
    self.temperature_probe=False
    self.directional_calib_made=False
    self.k1=0
    self.k2=0
    self.TI_UV=0
    
    #------Pro směrovou kalibraci------
    self.U1_aver=np.array([]) #Ve směru drátku
    self.U2_aver=np.array([]) #Ve směru drátku
    self.U1_RMS=np.array([]) #Ve směru drátku
    self.U2_RMS=np.array([]) #Ve směru drátku
    self.U_aver=np.array([]) #Ve směru sondy
    self.V_aver=np.array([]) #Ve směru sondy
    self.U_RMS=np.array([]) #Ve směru sondy
    self.V_RMS=np.array([]) #Ve směru sondy
    
    self.data_points=[] #List of all data points for the measurement
    self.text=''
    
    self.read_data(in_file) #Reads the data from text file    
    self.temper_corr(experiment_temper, wire_temper, calibration_temper) #Calculates the temperature correction
    
    #print "calibr_coeff =", calibr_coeff
    if calibr_coeff != -1:
      print 'Adjusting data with given calibration coefficients.' #'Provádím kalibraci dle externě zadaných koeficientů.'
      if (len(calibr_coeff)==2 and self.sensors_no == 2):
        self.sensors_calibration=calibr_coeff
      elif len(calibr_coeff) == 6:
        self.sensors_calibration=[calibr_coeff]
      
    self.calibration(self.sensors_calibration) #Converts voltage into velocity according to calibration coefficients
    
    if not (self.k1 == 0 and self.k2 == 0):
      if self.k1>5 or self.k2>5:
        self.directional_calibration(1/self.k1,1/self.k2) #TODO it is not necessary if directional calibration is correct
      else:
        self.directional_calibration(self.k1,self.k2)
    
    self.x=np.array([i.x for i in self.data_points]) #Numpy array of all positions (assumed only 1D traversing) TODO: Extend it into 3D
    self.prumery=np.array([np.average(i.data) for i in self.data_points]) #Numpy array of average velocities in each point
    self.RMS=np.array([np.std(i.data) for i in self.data_points]) #Numpy array of root mean square (np.std() ) in each point
    if self.sensors_no >1:
      self.prumery2=np.array([np.average(i.data2) for i in self.data_points]) #Numpy array of second wire average velocities in each point
      self.RMS2=np.array([np.std(i.data2) for i in self.data_points]) 
      self.TI=np.sqrt(1/2.*self.RMS**2+self.RMS2**2)/(np.sqrt(self.prumery**2+self.prumery2**2))
        
  def temper_corr(self, T1, Tw, T0):
    '''
    Performes data correction for different gas temperature than during calibration.
    
    T1 - [°C] temperature during experiment
    Tw - [°C] wire temperature
    T0 - [°C] reference temperature (during calibration)
    
    Ještě jiný postup je popsán v 
    [1]   HULTMARK, M. and SMITS, A.J. Temperature corrections for constant temperature and constant current hot-wire anemometers. In: Measurement Science and Technology. 1 October 2010, Vol. 21, no. 10, pp. 105404. DOI 10.1088/0957-0233/21/10/105404. .
    
    U/nu=f(E^2/(k*dT))
    dT=Tw-Ta 
    '''
    
    #Check for reasonable temperature user input
    if not((-20 < T1 < 300) and (0<Tw<800) and (-20<T0<300)):
      self.text += "Error!: Wrong temperature for temperature correction\n"
      self.corr_T = -1
      return -1    
    
    if T1 is not False:
      #m = 2e-5+5e-8*T1+3e-11*T1**2 # Temperature loading factor - from physical library [StreamWare Pro, Installation and User Guide], temperature input in [°C], page no. 266 - WRONG!!
      m=0.2
      #m=0.1 #Alternative option
      T1+=273.15
      Tw+=273.15
      T0+=273.15
      self.text += "m = "+str(m)+'\n'
            
      if T0 > T1:
        self.corr_T = ((Tw-T0)/(Tw-T1))**(0.5*(1-m))
      else:
        self.corr_T = ((Tw-T0)/(Tw-T1))**(0.5*(1+m))
    else:
      self.corr_T = 1 
      
  def __str__(self):
    text = '\n'
    text += ("Coefficient for temperature adjustment: %.2f->resulting change by %.2f %%\n")%(self.corr_T, 100-self.corr_T*100) #Koeficient pro teplotní kalibraci
    text += "Calibration coefficients C0..C5: "+str(self.sensors_calibration)+'\n'# C0,C1,C2,C3,C4,C5
    text += 'Number of measurement points x: '+str(len(self.data_points))+'\n' #data.shape[0]
    text += 'Number of samples at one measurement point: '+str(self.data_points[0].data.shape[0])+'\n' #data.shape[1]
    text += "Number of sensors (wires): "+ str(self.sensors_no)+'\n'
    if self.directional_calib_made:
      text += ('Directional calibration constants: k1 = %f; k2 = %f\n')%(self.k1, self.k2)
    
    if self.sensors_no == 1:
      text += '{:>5s} {:>9s} {:>7s}\n'.format('x', 'Average', 'TI')
      text += '{:>5s} {:>9s} {:>7s}\n'.format('[mm]', '[m/s]', '[%]')
      for i in range(len(self.data_points)):
        text += '{:>5g} {:>9.3f} {:>7.3f}\n'.format(self.x[i], self.prumery[i], self.RMS[i]/self.prumery[i]*100)
    elif self.sensors_no > 1:
      text += '{:>5s} {:>9s} {:>9s} {:>7s}\n'.format('x', 'Average', 'AVerage2', 'TI')
      text += '{:>5} {:>9} {:>9} {:>7}\n'.format('[mm]', '[m/s]', '[m/s]', '[%]')
      for i in range(len(self.data_points)):
        text += '{:>5g} {:>9.3f} {:>9.3f} {:>7.3f}\n'.format(self.x[i], self.prumery[i], self.prumery2[i], self.TI[i]*100)
      if self.directional_calib_made:
        text += '\n'
        text += '{:>5s} {:>9s} {:>9s} {:>7s}\n'.format('x', 'U', 'V', 'TI')
        text += '{:>5s} {:>9s} {:>9s} {:>7s}\n'.format('[mm]', '[m/s]', '[m/s]', '[%]')
        for i in range(len(self.data_points)):
          text += '{:>5g} {:>9.3f} {:>9.3f} {:>7.3f}\n'.format(self.x[i], self.U_aver[i], self.V_aver[i], self.TI_UV[i]*100)
    
    text +='\n'
    return self.text+text
  
  
  def read_data(self,in_file):
    '''
    Loads data from text file exported from Dantec's software StreamWare
    '''
    
    values=np.array([])
    values2=np.array([])
    values_time=np.array([])
    sw=-1
    sens_order=-1
    
    self.text +='\n'
    self.text += '----File: '+in_file+' ----'

    with open(in_file,'r') as fh:
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
          if values.shape[0]>0:
            self.data_points.append( DataPoint(values_time, values, ac_date, ac_time, x, y, z, position_no, values2) )
            #print "Hodnoty:", values[0]
            values=np.array([])
            values2=np.array([])
            values_time=np.array([])
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
              elif "	k =" in line:
                self.k1=float(line.split()[-2]) #TODO - should be k1**0.5?
                self.k2=float(line.split()[-1]) #TODO - should be k2**0.5?
                #print ("Načteno k1 = %f, k2 = %f")%(self.k1,self.k2)
            
            if 'No. of Sensors' in line:
              self.sensors_no=float(line.split()[-1])
              sens_order=0
            #Reads calibration coefficients C0..C5 for wire no 1
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
              sens_order=1 #Calibration data of the second wire can be loaded in next step
              self.sensors_calibration.append([C0,C1,C2,C3,C4,C5])
              #print "CAlibration coefficients: ", self.sensors_calibration
            
            
            
        if sw==3: #DATA HEADER
          #Reads info about probe position in current measurement point
          if line[0:15] == "(mm,mm,mm,deg):":
            x,y,z,deg = [float(i) for i in line.split()[1:]]
            #print "x,y,z:", x,y,z
          #Reads general info about the measurement set (ID, time, date, ...)
          elif 'Acquisition time:' in line:
            ac_time = line.split()[-1]
          
          elif 'Acquisition date:' in line:
            if not '.' in line.split()[-1]:
#              print 'There is no time specified!!!!!'
              ac_date = line.split()[-1]
                        
          elif 'Position no.:' in line:
            position_no= int(line.split()[-1])
            
        if sw == 4: #We are in the block of measured data
          if len(line.strip())>0 and line.strip()[0].isdigit(): # Is the first char on line the number?
            E = float(line.strip().split()[1]) * self.corr_T # Reads voltage and makes correction for temperature
            if self.sensors_no > 1:
              E2 = float(line.strip().split()[2]) * self.corr_T # Reads voltage and makes correction for temperature on second wire
            #C0,C1,C2,C3,C4,C5=self.sensors_calibration[0]
            #U=C0+C1*E+C2*E**2+C3*E**3+C4*E**4+C5*E**5 # Aplikována calibration
            values=np.append(values, E)
            values2=np.append(values2, E2)
            values_time=np.append(values_time, float(line.strip().split()[0]))
            #print values[-1]
          
          #elif ('[DATA HEADER' in line) and len(values)>0:
            #self.data_points.append( DataPoint(values_time, values, ac_date, ac_time, x, y, z, position_no, values2) )
            #print "Hodnoty:", values[0]
            #values=np.array([])
            #values2=np.array([])
            #values_time=np.array([])
    
    #Pass all the read data into the list of all data points
    self.data_points.append( DataPoint(values_time, values, ac_date, ac_time, x, y, z, position_no, values2) )
    
            
  def calibration(self,coeff):
    '''
    Converts voltage into valocity according to calibration coefficients.
    '''
    
    #------------Two-wire probe------------
    if len(coeff) == 2 and self.sensors_no ==2: 
      
      C0,C1,C2,C3,C4,C5=coeff[0]
      C0b,C1b,C2b,C3b,C4b,C5b=coeff[1]
      for dp in self.data_points:
        dp.data=C0+C1*dp.data+C2*dp.data**2+C3*dp.data**3+C4*dp.data**4+C5*dp.data**5 # Aplikována calibration
        dp.data2=C0b+C1b*dp.data2+C2b*dp.data2**2+C3b*dp.data2**3+C4b*dp.data2**4+C5b*dp.data2**5 # Aplikována calibration
    
    #------------One-wire probe------------
    elif len(coeff) == 1: 
      C0,C1,C2,C3,C4,C5=coeff[0]
      for dp in self.data_points:
        dp.data=C0+C1*dp.data+C2*dp.data**2+C3*dp.data**3+C4*dp.data**4+C5*dp.data**5 # Aplikována calibration
    #------------------------------------
        
    else:
      print 'Error!: Wrong input of calibration coefficients!' #'Error: Špatně zadané kalibrační koeficienty!'
      #raise NameError('Wrong input of calibration coefficients')
      
      
  def directional_calibration(self, k1, k2):
    '''
    Directional calibration according given coefficients - has to be called by user. It is not incorporated in the class run.
    
    k1, k2 - coefficients from velocity calibration    
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
      
    self.TI_UV=np.sqrt(1/2.*self.U_RMS**2+self.V_RMS**2)/(np.sqrt(self.U_aver**2+self.V_aver**2))
    
    self.directional_calib_made=True
      

  def fft_analysis(self,position=None):
    '''
    Performes Fast fourier analysis on data sets and returns either five most significant frequencies or numpy arrays of 
    frequencies and power spectrum when specific position number is passed as parameter.
    '''

    if position is None:
      freqs_5=np.array([])
      
      for dp, pr in zip(self.data_points, self.prumery):
        time_step=(dp.data_time[-1]-dp.data_time[0])/dp.data.size
        ps=np.abs(np.fft.rfft(dp.data-pr))**2 #Power spectrum
      
        freqs = np.fft.rfftfreq(dp.data.size, time_step)
        peakind = signal.find_peaks_cwt(ps, np.arange(1,50))
        freqs_5 = np.append(freqs_5,freqs[peakind[:5]])
        
      return freqs[peakind[:5]] #Returns first five most significant frequencies
    
    else: #Returns arrays of frequencies and power spectrum and five most significant frequencies.
      for dp,pr in zip(self.data_points, self.prumery):
        if position == dp.position_no:
          time_step=(dp.data_time[-1]-dp.data_time[0])/dp.data.size
          ps=np.abs(np.fft.rfft(dp.data-pr))**2 #Power spectrum
          #ps=np.angle(np.fft.rfft(dp.data-pr)) #Phase spectrum
          #ps=20*np.log10(np.abs(np.fft.rfft(dp.data-pr))) #Amplitude spectrum in [dB]
        
          freqs = np.fft.rfftfreq(dp.data.size, time_step)
          peakind = signal.find_peaks_cwt(ps, np.arange(1,50))
          
          
          return freqs, ps, freqs[peakind[:5]]
      
if __name__ == "__main__":
  import matplotlib.pyplot as plt
  import HotWire_calc as hw
  
  T1 = False #Temperature during experimet, put False to turn it off
  
  cta=[]
#  directory="./2015_02_12_2wire/"
#  files=['01_006_short.txt']
  directory='Viric_data/viric_240_25/0,1D/'
  file_in='240_25_01_00.txt'
  
  
#    cta.append(hw.HotWireMeasurement(directory+file_in,[[23.627550,-47.563171,36.200703,-13.087972,2.099148,0],[23.627550,-47.563171,36.200703,-13.087972,2.099148,0]]))
#    cta.directional_calibration(12.18**0.5, 15.6**0.5)    
  cta=hw.HotWireMeasurement(directory+file_in)
  peaks=cta.fft_analysis()
#    cta.directional_calibration()
  print cta #Vytiskne informace o daném měření
  print 'Significant frequencies: ',peaks,'Hz'
    
  plt.grid()
  plt.plot( cta.x ,cta.U_aver,'x-',label=file_in[0][:-4]+' U') #prumer/U0
  plt.plot( cta.x ,cta.V_aver,'x-',label=file_in[0][:-4]+' V') #prumer/U0
  plt.xlabel('r [mm]')
  plt.ylabel('U [m/s]')
  plt.legend()
  plt.rcParams.update({'font.size': 20})
  plt.show()
  
  fr,ps=cta.fft_analysis(3)[:2]
  plt.xlabel('Frequency [Hz]')
  plt.ylabel('Power spectrum')
  plt.plot(fr,ps)
  plt.show()