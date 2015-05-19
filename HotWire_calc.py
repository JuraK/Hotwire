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
  def __init__(self, in_file, calibr_coeff=-1, calibration_temper=17, experiment_temper=False, wire_temper=260):
    #self.DataPoints=[] #List of measurement points (data + info about the measurement point -time, date, position, no of position)
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
    self.calibration_temper=calibration_temper
    self.experiment_temper=experiment_temper
    self.wire_temper=wire_temper
    
    #------Pro směrovou kalibraci------
    self.U1_aver=np.array([]) #Ve směru drátku
    self.U2_aver=np.array([]) #Ve směru drátku
    self.U1_RMS=np.array([]) #Ve směru drátku
    self.U2_RMS=np.array([]) #Ve směru drátku
    self.U_aver=np.array([]) #Ve směru sondy
    self.V_aver=np.array([]) #Ve směru sondy
    self.U_RMS=np.array([]) #Ve směru sondy
    self.V_RMS=np.array([]) #Ve směru sondy
#    self.U_mag=np.array([]) #Velocity magnitude
#    self.phi=np.array([]) #Angle between U and V velocity components
    
    self.data_points=[] #List of all data points for the measurement
    self.text=''
#    self.nan_val_index=[] #Ratio of NaN count in directionally calibrated data to the total data count for every Data Point
    
    self.read_data(in_file) #Reads the data from text file    
    self.temper_corr(self.experiment_temper, wire_temper, self.calibration_temper) #Calculates the temperature correction
    
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
      self.text += "\n\nm = "+str(m)
            
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
#      text += ('NaN index in position %d is %g %%\n')%(self.x[np.argmax(self.nan_val_index)], self.nan_val_index[np.argmax(self.nan_val_index)]*100)
    
    if self.sensors_no == 1:
      text += '\n{:>5s} {:>9s} {:>7s}\n'.format('x', 'Average', 'TI')
      text += '{:>5s} {:>9s} {:>7s}\n'.format('[mm]', '[m/s]', '[%]')
      for i in range(len(self.data_points)):
        text += '{:>5g} {:>9.3f} {:>7.3f}\n'.format(self.x[i], self.prumery[i], self.RMS[i]/self.prumery[i]*100)
    elif self.sensors_no > 1:
      text += '\n{:>5s} {:>9s} {:>9s} {:>7s}\n'.format('x', 'Average', 'AVerage2', 'TI')
      text += '{:>5} {:>9} {:>9} {:>7}\n'.format('[mm]', '[m/s]', '[m/s]', '[%]')
      for i in range(len(self.data_points)):
        text += '{:>5g} {:>9.3f} {:>9.3f} {:>7.3f}\n'.format(self.x[i], self.prumery[i], self.prumery2[i], self.TI[i]*100)
      if self.directional_calib_made:
        text += '\n'
        text += '{:>5s} {:>9s} {:>9s} {:>7s} {:>7s} {:>7s}\n'.format('x', 'U', 'V', 'TI', 'Umag', 'PHI')
        text += '{:>5s} {:>9s} {:>9s} {:>7s} {:>7s} {:>7s}\n'.format('[mm]', '[m/s]', '[m/s]', '[%]', '[m/s]', '[°]')
        for i in range(len(self.data_points)):
          text += '{:>5g} {:>9.3f} {:>9.3f} {:>7.3f} {:>7.3f} {:>7.3f}\n'.format(self.x[i], self.U_aver[i], self.V_aver[i], self.TI_UV[i]*100, self.U_mag[i], self.phi[i])
    
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
          elif 'Measurement Temperature [°C]:' in line:
            self.experiment_temper=float(line.split(':')[-1].strip())
        
        if sw == 2:
          if 'teplot' in line:
            temperature_probe=True
            self.temperature_probe=True
          
          if not temperature_probe:
            if "Cal. ref. temp.:" in line:
              ll=line.split(':')[1].strip()
              if len(ll)>0:
                self.calibration_temper=float(ll)
            
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
                self.k1=float(line.split()[-2]) 
                self.k2=float(line.split()[-1])
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
      
      U1=(((0.5*(dp.data2**2)*(1+k2))-(0.5*(dp.data**2)*k2*(1+k1)))/(1-(k1*k2)))**0.5
      U2=(((0.5*(dp.data**2)*(1+k1))-(0.5*(dp.data2**2)*k1*(1+k2)))/(1-(k1*k2)))**0.5

      #------Check and correct NaN values in calculated velocities----
      nan_i=np.append(np.where(np.isnan(U1)),np.where(np.isnan(U2)))
      nan_i=np.unique(nan_i)
      U1[nan_i]=(((0.5*(dp.data2[nan_i]**2)*(1+k2)))/(1-(k1*k2)))**0.5
      U2[nan_i]=(((0.5*(dp.data[nan_i]**2)*(1+k1)))/(1-(k1*k2)))**0.5
      #TODO check the correction of thsi behaviour - ignoring one part of equation when negativ value under sqrt.
      #--------------------------------------------------------------

      U=((2**0.5)/2)*U1+((2**0.5)/2)*U2    #rychlost ve směru osy x
      V=((2**0.5)/2)*U2-((2**0.5)/2)*U1     #rychlost ve směru osy y
      
#      nan_val_count=np.max([np.count_nonzero(np.isnan(U1)),np.count_nonzero(np.isnan(U2))])
#      self.nan_val_index.append(nan_val_count/np.shape(U1)[0]) #Index of NaN values in velocity arrays
#      
      self.U1_aver=np.append(self.U1_aver,np.nanmean(U1))
      self.U2_aver=np.append(self.U2_aver,np.nanmean(U2))
      
      self.U1_RMS=np.append(self.U1_RMS,np.nanstd(U1))
      self.U2_RMS=np.append(self.U2_RMS,np.nanstd(U2))
      
      self.U_aver=np.append(self.U_aver,np.nanmean(U))
      self.V_aver=np.append(self.V_aver,np.nanmean(V))
      
      self.U_RMS=np.append(self.U_RMS,np.nanstd(U))
      self.V_RMS=np.append(self.V_RMS,np.nanstd(V))
      
      self.U_mag = (self.U_aver**2+self.V_aver**2)**0.5
      self.phi = np.degrees(np.arctan2(self.V_aver,self.U_aver)) #TODO check if it is correct!
      
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


def swirl_number(hws,R):
  '''
  Calculates the swirl number based on supplied list of Hot Wire measurements.
  Assumes radial position x [mm] measured from the wall to the center i.e. x[0]=2mm; x[-1]=150mm, R=150mm
  radius is then calculated as r=R-x
  R [mm] - pipe radius
  hws - list of HotWireMeasurement objects
  '''
  Gy=0
  Gx=0
  try:
    p=hws[0].positions
  except AttributeError: #Single data set supplied into hws
    print "Supplied single data set for swirl nuber calculation!"
    r=R-hws.x #TODO here it is pressumed, that x is sorted!
    dr=np.array([R-r[0]])
    for i in range(0,r.shape[0]-1):
      dr=np.append(dr,r[i]-r[i+1])
    for i in range(dr.shape[0]-1):
      Gy+=dr[i]/2*(hws.U_aver[i+1]*hws.V_aver[i+1]*r[i+1]**2+hws.U_aver[i]*hws.V_aver[i]*r[i]**2)
      Gx+=dr[i]/2*(hws.U_aver[i+1]**2*r[i+1]+hws.U_aver[i]**2*r[i])
    S=Gy/(R*Gx) #Swirl number based on Weber, Roman, and Jacques Dugué. 1992. 
                 #“Combustion Accelerated Swirling Flows in High Confinements.” 
                 #Progress in Energy and Combustion Science 18 (4): 349–67. 
                 #doi:10.1016/0360-1285(92)90005-L.
    return S
  
  
  #----Compare measured points if in same locations and proper order---
  hx=hws[0].x
  for h in hws[1:]:
    if (not p==h.positions) or (not (h.x==hx).all):
      print "Error: Positions are not the same in each data set, cannot calculate Swirl number!"
      return -1
    hx=h.x
  #--------------------------------------------------------------------
    
  count=0
  r=R-hws[0].x #TODO here it is pressumed, that x is sorted!

  #---Calculates delta r  ---
  dr=np.array([R-r[0]])
  for i in range(0,r.shape[0]-1):
      dr=np.append(dr,r[i]-r[i+1])
  #--------------------------
  
  U_av=0
  V_av=0
  for h in hws:
    U_av=((U_av*count)+h.U_aver)/(count+1) #Running average from different inclinations
    V_av=((V_av*count)+h.V_aver)/(count+1) #Running average from different inclinations
    count+=1
  for i in range(dr.shape[0]-1):
    Gy+=dr[i]/2*(U_av[i+1]*V_av[i+1]*r[i+1]**2+U_av[i]*V_av[i]*r[i]**2) #Numericall integration 
    Gx+=dr[i]/2*(U_av[i+1]**2*r[i+1]+U_av[i]**2*r[i])
  S=Gy/(R*Gx) #Swirl number based on Weber, Roman, and Jacques Dugué. 1992. 
                 #“Combustion Accelerated Swirling Flows in High Confinements.” 
                 #Progress in Energy and Combustion Science 18 (4): 349–67. 
                 #doi:10.1016/0360-1285(92)90005-L.
  return S
  
  
      
      
if __name__ == "__main__":
  import matplotlib.pyplot as plt
  import HotWire_calc as hw
  
  cta=[]
#  directory="./2015_02_12_2wire/"
#  files=['01_006_short.txt']
#  directory='Viric_data/viric_240_25/0,1D/'
#  file_in='240_25_01_00.txt'
#  directory='Viric_data/smazat/'
#  file_in='ověřovací měření - zpracované.txt'
  directory='Viric_data/viric_240_25/0,1D/'
  file_in='240_25_01_00.txt'
  
#    cta.append(hw.HotWireMeasurement(directory+file_in,[[23.627550,-47.563171,36.200703,-13.087972,2.099148,0],[23.627550,-47.563171,36.200703,-13.087972,2.099148,0]]))
#    cta.directional_calibration(12.18**0.5, 15.6**0.5)    
  cta=hw.HotWireMeasurement(directory+file_in, wire_temper=260)
  peaks=cta.fft_analysis()
  print cta #Vytiskne informace o daném měření
  print 'Significant frequencies: ',peaks,'Hz'
  
  print "Swirl number:",swirl_number(cta,150)
    
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