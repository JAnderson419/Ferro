import xlrd
import re
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib.ticker as plttick
#variable declaration
first_line=1 #initially set to 1
summary_table=0 # initially set to 0
is_data_table=0
is_data_table_header=1
c=0

#path of the folder for the .dat files.This has to be changed to the location where the files reside in the destination computer
path='C:/Users/anagha/Desktop/RIT/Project/Ferro-master/Ferro-master/ferro/tests/testData/RTWhiteB'
#enter the file name required with xtension
file_name=input("enter the file name ")
#iterating through the files to know if we have the file
#this will show "file found" if the file is there else wont show anything
for items in os.listdir(path):
   if(items==file_name):
      print('file found')
	  #opening the found file
      items=os.path.join(path,file_name)
      f=open(items,'r')
      # stripping the white spaces in opened file
      #creating the directory name for the data
      print(file_name)
      dir_name=file_name[0:len(file_name)-4]
      print(dir_name)
	  #making the directory 
      try:
         print(path)
		 #changin the directory from the python directory to the intended 
         os.chdir(path)
         os.mkdir(dir_name)
         print("directory",dir_name,"created")
         path=os.path.join(path,dir_name)
         os.chdir(path)
      except FileExistsError:
         print("directory",dir_name,"already exists")
      #opening the contents of the file
      for line in f:
         if(first_line):
            if re.match(r'DynamicHysteresisResult',line):
               type=0               			
            elif re.match(r'Fatigue',line):
               type=1
            elif re.match(r'PulseResult',line):
               type=2
            elif re.match(r'LeakageResult',line):
               type=3
            print(type)
            first_line=0
         else:
            #if(summary_table==1):
            #   continue
            if re.match(r"DynamicHysteresis",line) or re.match(r"Data Table",line) or re.match(r"Pulse",line) or re.match(r"Leakage",line):
               summary_table=1
               print(summary_table)
            elif re.match(r"Table(\[*\d*.*\d*\]*)",line): 
               group=re.compile(r"Table(\[*\d*.*\d*\]*)")
               result=group.match(line)
               c=c+1
               if result:
                  table_number=str(c)
                  table=result.group(0)
                  print(table)
 
            elif re.match(r"SampleName: (.*)",line):
               group1=re.compile(r"SampleName: (.*)")
               result1=group1.match(line)			   
               if result1:
                  samplename=result1.group(0)
                  print(samplename)
            elif re.match(r"Averages: (.*)",line):
               group7=re.compile(r"Averages: (.*)")
               result7=group7.match(line)
               if result7:
                  average=result7.group(0)
                  print(average)
            elif re.match(r"Total Cycles: (\d*)",line):
               group8=re.compile(r"Total Cycles: (\d*)")
               result8=group8.match(line)
               if result8:
                  total_cycles=result8.group(0)
                  print(total_cycles)

            elif (re.match(r"^Voltage \[V\]",line)) and (type==3):
               is_data_table=1
               print(is_data_table)
            elif (re.match("^Time \[s\]",line)) and (type==0 or type==1 or type==2):
               is_data_table=1
               print(is_data_table)
            elif (is_data_table):
               
           
			   
			   
			   
			   
			   
               if(type==0):
                  datapath1=str(samplename +hysteresisfrequency +"Hz"+hysteresisamplitude +"V" +""+average+""+"Average Table"+""+table_number +".tsv")
                  datapath=''.join(e for e in datapath1 if e.isalnum())
                  #is_data_table_header=1
				  
                  print(datapath)
                  if(is_data_table_header):
                     print(path)
                     datafile=open(datapath,"w")
                     is_data_table_header=0
                     re.sub("\[","''",line)
                     re.sub("\]","''",line)
                     re.sub("\/","_per_",line)
                     re.sub("\+","plus",line)
                     re.sub("-","minus",line)
                     print("inside the if statement")
                     datafile.write(line)
                  else:
                     datafile=open(datapath,"a")
                     re.sub("\[","''",line)
                     re.sub("\]","''",line)
                     re.sub("\/","_per_",line)
                     re.sub("\+","plus",line)
                     re.sub("-","minus",line)
                     datafile.write(line)
                     print("it has jumped to the else statement")
                  if re.match(r"^$",line):
                     is_data_table=0
                     is_data_table_header=1
                     datafile.close()
                     print("it has reached this if")
               elif(type==1):
                  datapath2=str(samplename +hysteresisfrequency +"Hz"+hysteresisamplitude +"V" +average +"Average Table" + table + total_cycles+".tsv")
                  datapath=''.join(e for e in datapath2 if e.isalnum())
                  print(datapath)
                  if(is_data_table_header):
                     print(path)
                     datafile=open(datapath,"w")
                     is_data_table_header=0
                     re.sub("\[","''",line)
                     re.sub("\]","''",line)
                     re.sub("\/","_per_",line)
                     re.sub("\+","plus",line)
                     re.sub("-","minus",line)
                     print("inside the if statement")
                     datafile.write(line)
                  else:
                     datafile=open(datapath,"a")
                     re.sub("\[","''",line)
                     re.sub("\]","''",line)
                     re.sub("\/","_per_",line)
                     re.sub("\+","plus",line)
                     re.sub("-","minus",line)
                     datafile.write(line)
                     print("it has jumped to the else statement")
                  if re.match(r"^$",line):
                     is_data_table=0
                     is_data_table_header=1
                     datafile.close()
                     print("it has reached this if")
                  
               elif(type==2):

                  datapath3=str(samplename +hysteresisfrequency +"Hz"+hysteresisamplitude +"V" +average +"Average Table" + table + total_cycles+".tsv")
                  datapath=''.join(e for e in datapath3 if e.isalnum())
                  print(datapath)
                  if(is_data_table_header):
                     print(path)
                     datafile=open(datapath,"w")
                     is_data_table_header=0
                     re.sub("\[","''",line)
                     re.sub("\]","''",line)
                     re.sub("\/","_per_",line)
                     re.sub("\+","plus",line)
                     re.sub("-","minus",line)
                     print("inside the if statement")
                     datafile.write(line)
                  else:
                     datafile=open(datapath,"a")
                     re.sub("\[","''",line)
                     re.sub("\]","''",line)
                     re.sub("\/","_per_",line)
                     re.sub("\+","plus",line)
                     re.sub("-","minus",line)
                     datafile.write(line)
                     print("it has jumped to the else statement")
                  if re.match(r"^$",line):
                     is_data_table=0
                     is_data_table_header=1
                     datafile.close()
                     print("it has reached this if")
                  
                  
               elif(type==3):

                  datapath4=str(samplename +stepduration+'s steps'+ table +".tsv")
                  datapath=''.join(e for e in datapath4 if e.isalnum())
                  print(datapath)
                  if(is_data_table_header):
                     print(path)
                     datafile=open(datapath,"w")
                     is_data_table_header=0
                     re.sub("\[","''",line)
                     re.sub("\]","''",line)
                     re.sub("\/","_per_",line)
                     re.sub("\+","plus",line)
                     re.sub("-","minus",line)
                     print("inside the if statement")
                     datafile.write(line)
                  else:
                     datafile=open(datapath,"a")
                     re.sub("\[","''",line)
                     re.sub("\]","''",line)
                     re.sub("\/","_per_",line)
                     re.sub("\+","plus",line)
                     re.sub("-","minus",line)
                     datafile.write(line)
                     print("it has jumped to the else statement")
                  if re.match(r"^$",line):
                     is_data_table=0
                     is_data_table_header=1
                     datafile.close()
                     print("it has reached this if")
				  
               

				  
			   
			   
            
               

                  			   
            elif (type==3):
               group2=re.compile(r"Step Duration \[s\]: (\d*)")
               result2=group2.match(line)
               
               if result2:
                  stepduration=result2.group(0)
                  print(stepduration)
            elif (type==0 or type==1):
               group3=re.compile(r"Hysteresis Frequency \[Hz\]: (\d*)")
               group4=re.compile(r"Hysteresis Amplitude \[V\]: (\d*)")
               result4=group4.match(line)
               result3=group3.match(line)
               if result3:
                  hysteresisfrequency=result3.group(0)
                  print(hysteresisfrequency)
               elif result4:
                  hysteresisamplitude=result4.group(0)
                  print(hysteresisamplitude)
            elif(type==2):
               group5=re.compile(r"Pund Frequency \[Hz\]: (\d*)")
               group6=re.compile(r"Pund Amplitude \[V\]: (\d*)")
               result5=result5.match(line)
               result6=result6.match(line)
               if result5:
                  frequency=result5.group(0)
                  print(frequency)
               elif result6:
                  amplitude=result6.group(0)
				  amplitude=result6.group(0)
				  print(amplitude)
         
         
         			
            		    
            
            
         
  
  
   
   