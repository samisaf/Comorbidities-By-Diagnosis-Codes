import pandas as pd
import numpy as np
import re

Comorbidities = {"MayocardialInfarction": ["410\..*", "412\..*"], 
                 "HeartFailure": ["398\.91", "402\.01", "402\.11", "402\.91", "404\.01", "404\.03", "404\.11", "404\.13", "404\.91", "404\.93", "425\..*", "428\..*"],
                 "PeripheralVascularDisease": ["093\.0", "437\.3", "440\..*", "441\..*", "443\..*", "47\.1", "557\.1", "557\.9", "V43\.4"],
                 "CerebrovascularDisease":  ["362\.34", "430\..*", "431\..*", "432\..*", "433\..*", "434\..*", "435\..*", "436\..*", "437\..*", "438\..*"], 
                 "Dementia": ["290\..*", "294.1", "331.2"], 
                 "ChronicPulmonaryDisease": ["416.8", "416.9", "490\..*", "491\..*", "492\..*", "493\..*", "494\..*", "495\..*", "496\..*", "497\..*", "498\..*", "499\..*", "500\..*", "500\..*", "501\..*", "502\..*", "503\..*", "504\..*", "505\..*", "506.4", "508.1", "508.8"], 
                 "RheumaticDisease": ["446.5", "710\..*", "714.0", "714.01", "714.2", "714.8", "725\..*"], 
                 "PepticUlcerDisease": ["531\..*", "532\..*", "533\..*", "534\..*"], 
                 "MildLiverDisease": ["070.22", "070.23", "070.32", "070.33", "070.44", "070.54", "070.6", "070.9", "570\..*", "571\..*", "573.3", "573.4", "573.8", "573.9", "V42.7"], 
                 "DiabetesWithoutComplications": ["250.0", "250.1", "250.2", "250.3", "250.8", "250.9"], 
                 "DiabetesWithComplications": ["250.4", "250.5", "250.6", "250.7"], 
                 "HemiplegiaParaplegia": ["334.1", "342\..*", "343\..*", "344\..*"], 
                 "RenalDisease": ["403.01", "403.11", "403.91", "404.02", "404.03", "404.12", "404.13", "404.92", "404.93", "582\..*", "583\..*", "585\..*", "586\..*", "588.0", "V42.0", "V45.1", "V56\..*"], 
                 "Malignancy": ["14\d\..*", "15\d\..*", "16\d\..*", "170\..*", "171\..*", "172\..*", "174\..*", "175\..*", "176\..*", "177\..*", "178\..*", "179\..*", "18\d\..*", "190\..*", "191\..*", "192\..*", "193\..*", "194\..*", "195\..*", "20\d\..*", "238.6"], 
                 "ModerateSevereLiverDisease": ["456.0", "456.1", "456.2", "572\..*"], 
                 "MetastaticCancer": ["196\..*", "197\..*", "198\..*", "199\..*"], 
                 "HIVAIDS": ["042\..*", "043\..*", "044\..*"], 
                 "HypertensionUncomplicated": ["401\..*"],  
                 "HypertensionComplicated": ["402\..*", "403\..*", "404\..*", "405\..*"],
                 "Coagulopathy":  ["286\..*", "287.1", "287.3", "287.4", "287.5"],
                 "Obesity": ["278.0"],
                 "WeightLoss": ["260\..*", "261\..*", "262\..*", "263\..*", "783.2", "799.4"],
                 "FluidElectrolyteDisorders": ["253.6", "276\..*"]
                }
                
                
Weights = {"MayocardialInfarction": 1, 
                 "HeartFailure": 1, 
                 "PeripheralVascularDisease": 1, 
                 "CerebrovascularDisease":  1, 
                 "Dementia": 1, 
                 "ChronicPulmonaryDisease": 1, 
                 "RheumaticDisease": 1, 
                 "PepticUlcerDisease": 1, 
                 "MildLiverDisease": 1, 
                 "DiabetesWithoutComplications": 1, 
                 "DiabetesWithComplications": 2, 
                 "HemiplegiaParaplegia": 2, 
                 "RenalDisease": 2, 
                 "Malignancy": 2, 
                 "ModerateSevereLiverDisease": 3, 
                 "MetastaticCancer": 6, 
                 "HIVAIDS": 6, 
                 "HypertensionUncomplicated": 0, 
                 "HypertensionComplicated": 0, 
                 "Coagulopathy": 0, 
                 "Obesity": 0, 
                 "WeightLoss": 0,  
                 "FluidElectrolyteDisorders": 0
                }
Names = ["MayocardialInfarction", "HeartFailure", "PeripheralVascularDisease", "CerebrovascularDisease", "Dementia", "ChronicPulmonaryDisease", "RheumaticDisease", "PepticUlcerDisease", 
         "MildLiverDisease", "DiabetesWithoutComplications", "DiabetesWithComplications", "HemiplegiaParaplegia", "RenalDisease", "Malignancy", "ModerateSevereLiverDisease", "MetastaticCancer", "HIVAIDS", 
         "HypertensionUncomplicated", "HypertensionComplicated", "Coagulopathy", "Obesity", "WeightLoss",  "FluidElectrolyteDisorders"]   


Nos =  ["No"] * len(Names)
         
class Patient(object):
    def __init__(self, mrn, codes = [], age = 0, gender = 0, race = 0, encounterDate = 0):
        self.mrn = mrn
        self.age = age
        self.gender = gender
        self.race = race
        self.encounterDate = pd.to_datetime(encounterDate)
        self.codes = pd.DataFrame(codes, columns=['ICD9', 'date'])
        self.codes.date = pd.to_datetime(self.codes.date)
        self.comorbidities = []
        self.score = 0
        if self.codes.size > 0: 
            self.comorbidities= self.__calcComorbidities__()
            self.score = self.__calcScoreWOAge__()
            self.scoreAge = self.__calcScoreWAge__()
            self.comorbiditiesFrame = self.__createComFrame__()
                    
    def __calcComorbidities__(self):
        categories = Comorbidities.keys() # Comorbid categories such as MI, CHF, ...
        for cat in categories:
            icd9patterns = Comorbidities[cat] # ICD9 patterns such as 410\..* and 412\..* for MI
            for icd9pat in icd9patterns:
                pat = re.compile(icd9pat)
                matches = self.codes.ICD9.apply(pat.match)
                anymatch = np.any(matches)
                if anymatch != None: self.comorbidities.append(cat)
        self.comorbidities = set(self.comorbidities)
        return list(self.comorbidities)
        
    def __createComFrame__(self):
        temp = pd.DataFrame(data=[Nos], columns= Names)
        for name in self.comorbidities:
            temp[name][0] = "Yes"
        temp['CharlsonScore'] = self.score
        temp['CharlsonScoreAge'] = self.scoreAge
        temp['Age'] = self.age
        temp['MRN'] = self.mrn
        return temp
    
    def __calcScoreWAge__(self):
        weights = [Weights[c] for c in self.comorbidities]
        agescore = 0 
        if self.age > 40: agescore = (self.age - 40) // 10
        return int(sum(weights) + agescore)

    def __calcScoreWOAge__(self):
        weights = [Weights[c] for c in self.comorbidities]
        return int(sum(weights))
        
    def __str__(self): 
        temp = "<MRN {}, Age {}, Gender {}, Race {}>"\
            .format(self.mrn, self.age, self.gender, self.race)
        return temp

    def __repr__(self): return self.__str__()
        
if __name__ == "__main__":
    print ("Running some tests")
    mrn = 1
    age = 80
    race = 'C'
    gender = 'M'
    encounterDate = "2014-03-01"
    codes = [["117.9", "2014-01-01"],
            ["197.7", "2014-01-01"],
            ["198.5", "2014-01-01"],
            ["174.1", "2014-01-01"],
            ["174.9", "2014-01-01"],
            ["197.2", "2014-01-01"],
            ["199.1", "2014-01-01"],
            ["242.90", "2014-01-01"],
            ["244.9", "2014-01-01"],
            ["380.4", "2014-01-01"],
            ["388.70", "2014-01-01"],
            ["397.0", "2014-01-01"],
            ["401.9", "2014-01-01"],
            ["425.4", "2014-01-01"],
            ["424.0", "2014-01-01"],
            ["427.89", "2014-01-01"],
            ["429.3", "2014-01-01"],
            ["427.9", "2014-01-01"],
            ["440.0", "2014-01-01"],
            ["466.0", "2014-01-01"],
            ["470", "2014-01-01"],
            ["471.0", "2014-01-01"],
            ["461.9", "2014-01-01"],
            ["471.9", "2014-01-01"],
            ["473.0", "2014-01-01"],
            ["473.1", "2014-01-01"],
            ["477.2", "2014-01-01"],
            ["478.19", "2014-01-01"],
            ["473.2", "2014-01-01"],
            ["473.3", "2014-01-01"],
            ["473.8", "2014-01-01"],
            ["477.9", "2014-01-01"],
            ["478.0", "2014-01-01"],
            ["511.81", "2014-01-01"],
            ["511.9", "2014-01-01"],
            ["512.89", "2014-01-01"],
            ["515", "2014-01-01"],
            ["518.0", "2014-01-01"],
            ["518.89", "2014-01-01"],
            ["564.00", "2014-01-01"],
            ["575.9", "2014-01-01"],
            ["573.8", "2014-01-01"],
            ["579.0", "2014-01-01"],
            ["593.2", "2014-01-01"],
            ["611.72", "2014-01-01"],
            ["719.45", "2014-01-01"],
            ["733.00", "2014-01-01"],
            ["733.90", "2014-01-01"],
            ["781.1", "2014-01-01"],
            ["780.79", "2014-01-01"],
            ["784.42", "2014-01-01"],
            ["784.49", "2014-01-01"],
            ["784.91", "2014-01-01"],
            ["783.21", "2014-01-01"],
            ["786.2", "2014-01-01"],
            ["786.50", "2014-01-01"],
            ["786.09", "2014-01-01"],
            ["787.02", "2014-01-01"],
            ["789.00", "2014-01-01"],
            ["789.06", "2014-01-01"],
            ["790.21", "2014-01-01"],
            ["790.29", "2014-01-01"],
            ["790.6", "2014-01-01"],
            ["793.2", "2014-01-01"],
            ["793.80", "2014-01-01"],
            ["794.31", "2014-01-01"],
            ["796.2", "2014-01-01"],
            ["V04.81", "2014-01-01"],
            ["V05.9", "2014-01-01"],
            ["V06.1", "2014-01-01"],
            ["V10.3", "2014-01-01"],
            ["V12.2", "2014-01-01"],
            ["V12.29", "2014-01-01"],
            ["V58.69", "2014-01-01"],
            ["V70.7", "2014-01-01"],
            ["V72.81", "2014-01-01"],
            ["V72.83", "2014-01-01"],
            ["V72.84", "2014-01-01"],
            ["V76.12", "2014-01-01"],
            ["V77.1", "2014-01-01"],
            ["V81.2", "2014-01-01"],
            ["V82.81", "2014-01-01"]]   

    p1 = Patient(mrn, codes, age, gender, race, encounterDate)       
    print(p1)
    print(p1.comorbidities)
    print(p1.score)   
    print(p1.comorbiditiesFrame)  