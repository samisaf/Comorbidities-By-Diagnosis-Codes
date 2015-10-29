# -*- coding: utf-8 -*-
"""@author: Sami Safadi"""

from Patient import Patient
import pandas as pd
import numpy as np
import glob

demographicsFrame = pd.DataFrame(columns = ['MRN', 'age', 'gender', 'race', 'DOH'])
ptDict ={} # List containing all code values
Patients = {} # Dictonary, keys are MRN, values are Patient objects
        
def getNumPatients():
    return len(Patients)
    
def getNumCodes():
    return sum([Patients[p].codes.size for p in Patients])

def stageCKD(egfr):
    if egfr < 15: return 5
    elif egfr >= 15 and egfr < 30: return 4
    elif egfr >= 30 and egfr < 60: return 3
    elif egfr >= 60 and egfr < 90: return 2
    else: return np.NaN
     
def createPtDict(df: pd.DataFrame):
    global ptDict 
    uniqueMRN = np.unique(df.ix[:, 0])
    for i in uniqueMRN: 
        if not(i in ptDict.keys()): ptDict[i] = list()
    for i in range(len(df)):
        mrn = df.ix[i, 0]
        code = df.ix[i, 1]
        date = df.ix[i, 2]
        if not(np.isnan(mrn)): 
            ptDict[int(mrn)].append((code, str(date)))
    return ptDict

def createPts(patients: dict):
    global Patients
    for key in patients:
        codes = patients[key]
        mrn = key
        age, gender, race = 0, 0, 0
        if any(demographicsFrame.MRN == key):
            index = demographicsFrame.index[demographicsFrame.MRN == key][0]
            age = demographicsFrame.Age[index]
            gender = demographicsFrame.Gender[index]
            race = demographicsFrame.Race[index]
            encounterDate = demographicsFrame.DOH[index]
        if len(codes) > 0: Patients[mrn] = Patient(mrn, codes, age, gender, race, encounterDate)

def getMasterTable(savetofile = False):
    global Patients
    keys = list(Patients.keys())
    df = Patients[keys[0]].comorbiditiesFrame
    
    for key in keys[1:]:
        df = df.append(Patients[key].comorbiditiesFrame)
    
    outputfile = "Output/Comorbidities.csv"
    if savetofile: df.to_csv(outputfile, index = False)
    else: print(df)
    return df

def readDemographics(file: str):
    global demographicsFrame
    demographicsFrame = pd.read_csv(file)
        
def readCodes(files:[str]):
    for file in files: 
        createPtDict(pd.read_csv(file))
    createPts(ptDict)
    
def write():
    global Patients
    getMasterTable(savetofile = True)
    
if __name__ == "__main__":
    demographicsFile = "Input/Demographics.csv"
    try: readDemographics(demographicsFile)
    except: print("Demographics.csv was not found") 
    message =  "Loaded demographics file and identified {} patients...".format(demographicsFrame.MRN.size)   
    print(message)
    
    ICD9Files = glob.glob("Input/Codes*.csv")
    print("The following files are proccessed: ", ICD9Files)
    readCodes(ICD9Files)
    
    message = "Processed {} patients with {} comorbidities"\
        .format(getNumPatients(), getNumCodes(),)
    print(message)
    
    print("Writing results to disk")    
    write()
    input("Done, press enter to exit...")
