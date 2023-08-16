import numpy as np
import matplotlib.pyplot as plt
import sys
#Module for computing compliance with MC2010
#TODO: implement lightweight concrete


"""

"""

class creepMC2010():
    def __init__(self, Ac, u, strengthClass, cementType, RH, modulusCalculatedPerCode:bool=True, considerTemperatureEffect:bool=True, aggregateType:str=None, α_ThreeParModulus:float=None, β_ThreeParModulus:float=None, 𝜏_ThreeParModulus:float=None, verbose:bool=False):
        """
        Class constructor.
        This class defines creep-related properties for creep coefficient calculation according to Model Code 2010

        Parameters
        ----------
        Ac (int): Cross sectional area in mm2 (5.1-71d)
        u (int): Cross sectional perimeter in contact with drying, in mm (5.1-71d)
        strengthClass (str): Strength class accordingly to Model Code 2010 (check Table 5.1-3)
        cementType (str): Cement type accordingly to Model Code 2010 (check Table 5.1-9)
        Rh (float): Average relative humidity of the ambient,  in %
        modulusCalculatedPerCode (bool): Inform whether E-modulus should be calculated using MC2010 formulation. If "False", parameters of a three-parameter model must be informed
        considerTemperatureEffect (bool): Inform whether temperature effects and type of cement must be taken into account for correcting t0 (age of loading, as per Eq. B.10)
        aggregateType (str): Aggregate type accordingly to types addressed by MC2010. Influences MC2010 formulation for E-modulus estimation. Check Table 5.1-6. Options: "basalt", "quartzite", "limestone", "sandstone"
        α_ThreeParModulus (float): The alfa parameter of a three parameter model used to compute E-modulus evolution. Model must give modulus in MPa and time must be in hours. Used if modulusCalculatedPerCode is False.
        β_ThreeParModulus (float): The beta parameter of a three parameter model used to compute E-modulus evolution. Model must give modulus in MPa and time must be in hours. Used if modulusCalculatedPerCode is False.
        𝜏_ThreeParModulus (float): The tau parameter of a three parameter model used to compute E-modulus evolution. Model must give modulus in MPa and time must be in hours. Used if modulusCalculatedPerCode is False.
        verbose (bool): Controls whether some warning messages addressing implementation aspects are displayed.
        """

        #Geometrical information about the member
        self.Ac=Ac #5.1-71d, cross sectional area in mm2 
        self.u=u #5.1-71d, cross sectional permiter in contract with the atmosphere in mm
        self.h=2*self.Ac/self.u #Notional size
        
        #Validate cement type
        if cementType not in ["32.5N","32.5R","42.5N","42.5R","52.5N","52.5R"]:
            raise SystemExit("Cement type does not exist. Check Table 5.1-9")
        else:
            self.cementType = cementType

        #Concrete strengths
        #Characteristics compressive strength
        if strengthClass not in ["C12","C16","C20","C25","C30","C35","C40","C45","C50","C55","C60","C70","C80","C90","C100","C110","C120"]:
            print("Strength class does not exist. Check Table 5.1-3")
            exit()
        else:
            self.fck=int(strengthClass[-2:])  
        #Mean value of compressive strength      
        self.fcm = self.fck+8 #Eq. 5.1-1
        
        #Table 5.1-9: s-parameter
        if self.fcm <= 60:
            if self.cementType == "32.5N":
                self.s=0.38
                self.α=-1 #from Eq. 5.1-73
            elif self.cementType == "32.5R" or self.cementType == "42.5N":
                self.s=0.25
                self.α=0 #from Eq. 5.1-73
            elif self.cementType == "42.5R" or self.cementType == "52.5N" or self.cementType == "52.5R":
                self.s=0.20
                self.α=1 #from Eq. 5.1-73
        elif self.fcm > 60:
            self.s=0.20
            if self.cementType == "32.5N":
                self.α=-1 #from Eq. 5.1-73
            elif self.cementType == "32.5R" or self.cementType == "42.5N":
                self.α=0 #from Eq. 5.1-73
            elif self.cementType == "42.5R" or self.cementType == "52.5N" or self.cementType == "52.5R":
                self.α=1 #from Eq. 5.1-73

        #Calculation of modulus of elasticity
        self.modulusCalculatedPerCode=modulusCalculatedPerCode
        if self.modulusCalculatedPerCode is True:
            #Select parameter αE
            if aggregateType == "basalt":
                self.αE=1.2
            elif aggregateType == "quartzite":
                self.αE=1
            elif aggregateType == "limestone":
                self.αE=0.9
            elif aggregateType == "sandstone":
                self.αE=0.7
            else:
                print("Wrong type of aggregate. Check Table 5.1-6")

            try:
                self.Eci = 21.5*1000*(self.αE)*((self.fck+8)/10)**(1/3) #Eq. 5.1-20 (tangent E-modulus at initial of stress-strain diagram and approximately equal to the slope of the secant of the unloading branch for rapid unloading)
            except AttributeError as errorMessage:
                if errorMessage.args[0] == "'creepMC2010' object has no attribute 'αE'":
                    print("If using modulus calculated per code, you need to specificy aggregate type")
                    exit()
        else:
            self.α_ThreeParModulus=α_ThreeParModulus
            self.β_ThreeParModulus=β_ThreeParModulus
            self.𝜏_ThreeParModulus=𝜏_ThreeParModulus
            self.Eci = self.Eci_computeEvolution(28) #in MPa
            
        #Additional parameters of the code
        self.RH = RH #In %
        self.considerTemperatureEffect = considerTemperatureEffect #Boolean

        #Flags to display only once certain warning messages to the user, to clarify usage aspects
        self.creepTemperatureEffectsWarning=False
        self.verbose=verbose

    #Methods to compute creep properties
    def ϕ_compute (self, t, t0, timeHistory:list=None, temperatureHistory:list=None):
        """
        Parameters:
            timeHistory: list of floats
                List with consecutive time duration, in days, associated to each isothermal temperature in temperatureHistory, until the desired age of calculation.
                Ex: if we want to compute the ϕbc at age of 30 days of a material that initially stayed 20 days under 25C and then 10 days under 30C, then timeHistory = [20,10]
            temperatureHistory: list of floats
                List with the value of isothermal temperatures, in Celsius, during each interval of time 
                Ex: In the previous example at timeHistory description, the respective timeHistory would be [25,30]
        """
        #Temperature correction due to transient thermal creep (section 5.1.10.7.1)
        if self.considerTemperatureEffect is True:
            #Locate the current time and the loading in the temperature history
            temperatureHistoryCurrentTime,timeHistoryCurrentTime=self.locateTimeInTemperatureHistory(t, temperatureHistory, timeHistory)
            #For each temperature in temperatureHistoryCurrentTime, there is need for adding the correspondent transient thermal creep coefficient in case a temperature increase occurs in relation to the previous temperature
            #Of course, if a single temperature exist in temperatureHistoryCurrentTime, then no thermal transient occurs
            ΔϕTtrans=0
            currentTemperature = temperatureHistoryCurrentTime[0]
            for temperature in temperatureHistoryCurrentTime:
                if temperature > currentTemperature:
                    currentTemperature = temperature
                    ΔϕTtrans = ΔϕTtrans + 0.0004*(temperature-20)**2
        else:
            ΔϕTtrans=0
        return self.ϕbc_compute(timeHistory, temperatureHistory, t, t0)+self.ϕdc_compute(timeHistory, temperatureHistory, t, t0)+ΔϕTtrans

    def ϕbc_compute (self, timeHistory, temperatureHistory, t, t0):
        '''
        Compute the parameter ϕbc (basic creep coefficient) at an age t due to a loading at t0. According to Eq. 5.1-64. The age t is contained in the variable timeHistory.

        Parameters:
            timeHistory: list of floats
                List with consecutive time duration, in days, associated to each isothermal temperature in temperatureHistory, until the desired age of calculation.
                Ex: if we want to compute the ϕbc at age of 30 days of a material that initially stayed 20 days under 25C and then 10 days under 30C, then timeHistory = [20,10]
            temperatureHistory: list of floats
                List with the value of isothermal temperatures, in Celsius, during each interval of time 
                Ex: In the previous example at timeHistory description, the respective timeHistory would be [25,30]
            t: float
                Age in which the creep function is being evaluated, in days
            t0: float
                Loading age, in days.

        Returns:
            ϕbc: float
                Basic creep coefficient according to Eq. 5.1-64
        '''
        #Eq. 5.1-65
        βbc_fcm = 1.8/(self.fcm)**0.7
        #Eq. 5.1-66
        t0_adj = self.t0_adj_compute (t0, timeHistory, temperatureHistory)
        βbc_t_t0 = np.log((t-t0)*(((30/t0_adj)+0.035)**2)+1)
        #Eq. 5.1-64
        ϕbc=βbc_fcm*βbc_t_t0

        #Methods to compute temperature effects on creep (5.1.10.7)
        if self.considerTemperatureEffect is True:
            if self.creepTemperatureEffectsWarning is False and self.verbose is True:
                print("ϕbc_compute method: When considering temperature effects on creep (section 5.1.10.7), only the case of constant temperature is addressed by the code. Because of this, in case of temperature histories with variable temperatures, this code takes the time-weighted average value during the entire history until this moment.")
                self.creepTemperatureEffectsWarning = True
            #Locate the current time in the temperature history
            temperatureHistoryNew,timeHistoryNew=self.locateTimeInTemperatureHistory(t, temperatureHistory, timeHistory)
            Tavg=np.average(temperatureHistoryNew,weights=timeHistoryNew)
            ϕT=np.exp(0.015*(Tavg-20))
        else:
            ϕT=1

        return ϕbc*ϕT

    def ϕdc_compute (self, timeHistory, temperatureHistory, t, t0):
        '''
        Compute the parameter ϕdc (drying creep coefficient) at an age t due to a loading at t0. According to Eq. 5.1-67. The age t is contained in the variable timeHistory.

        Parameters:
            timeHistory: list of floats
                List with consecutive time duration, in days, associated to each isothermal temperature in temperatureHistory, until the desired age of calculation.
                Ex: if we want to compute the ϕbc at age of 30 days of a material that initially stayed 20 days under 25C and then 10 days under 30C, then timeHistory = [20,10]
            temperatureHistory: list of floats
                List with the value of isothermal temperatures, in Celsius, during each interval of time 
                Ex: In the previous example at timeHistory description, the respective timeHistory would be [25,30]
            t: float
                Age in which the creep function is being evaluated, in days
            t0: float
                Loading age, in days.
        
        Returns:
            ϕdc: float
                Drying creep coefficient according to Eq. 5.1-64
        '''
        #Eq. 5.1-68
        βdc_fcm = 412/(self.fcm)**1.4
        #Eq. 5.1-69
        βdc_RH = (1-self.RH/100)/(0.1*self.h/100)**(1/3)
        #Eq. 5.1-70
        t0_adj = self.t0_adj_compute (t0, timeHistory, temperatureHistory)
        βdc_t0 = 1/(0.1+t0_adj**0.2)

        #Eq. 5.1-71b
        γ_t0 = 1/(2.3 + 3.5/t0_adj**0.5)
        #Eq. 5.1-71d
        α_fcm = (35/self.fcm)**0.5
        #Eq. 5.1-71c
        βh = min(1.5*self.h+250*α_fcm, 1500*α_fcm)
        #Methods to compute temperature effects on creep (5.1.10.7)
        if self.considerTemperatureEffect is True:
            if self.creepTemperatureEffectsWarning is False and self.verbose is True:
                print("ϕdc_compute method: When considering temperature effects on creep (section 5.1.10.7), only the case of constant temperature is addressed by the code. Because of this, in case of temperature histories with variable temperatures, this code takes the time-weighted average value during the entire history until this moment.")
                self.creepTemperatureEffectsWarning = True
            #Locate the current time in the temperature history
            temperatureHistoryNew,timeHistoryNew=self.locateTimeInTemperatureHistory(t, temperatureHistory, timeHistory)
            Tavg=np.average(temperatureHistoryNew,weights=timeHistoryNew)
            βT=np.exp((1500/(273+Tavg))-5.12)
            ϕT=np.exp(0.015*(Tavg-20))
        else:
            βT=1
            ϕT=1
        #Eq. 5.1-71a
        βdc_t_t0 = ((t-t0)/(βh*βT+(t-t0)))**γ_t0

        #Eq. 5.1-67
        ϕdc=βdc_fcm*βdc_RH*βdc_t0*βdc_t_t0

        return ϕdc*(ϕT**1.2)
    
    def t0_adj_compute (self, t0, timeHistory, temperatureHistory):
        '''
        Method to compute t0_adj according to Eq. 5.1-73.
        The t0_adj computes the loading age t0, according to temperature effects and the effect of the type of cement on the creep coefficient of concrete.
        t0_adj is computed in days.

        Parameters:
            t0: float
                Loading age in days
            timeHistory: list of floats
                List with consecutive time duration, in days, associated to each isothermal temperature in temperatureHistory, until the desired age of calculation.
                Ex: if we want to compute the ϕbc at age of 30 days of a material that initially stayed 20 days under 25C and then 10 days under 30C, then timeHistory = [20,10]
            temperatureHistory: list of floats
                List with the value of isothermal temperatures, in Celsius, during each interval of time 
                Ex: In the previous example at timeHistory description, the respective timeHistory would be [25,30]
        
        Returns:
            t0_adj: float
                Adjusted time in days for creep and temperature effects.
        '''
        if self.considerTemperatureEffect is True:
            t0_T = self.adjustedAge(t0, timeHistory, temperatureHistory)
        else:
            t0_T = t0
        t0_adj = t0_T*((9/(2+t0_T**1.2))+1)**self.α
        t0_adj = max(t0_adj, 0.5)
        return t0_adj

    #Methods to compute ageing effects on properties
    def adjustedAge (self, t, timeHistory, temperatureHistory):
        '''
        Section 5.1.10.2: Maturity
        A method to apply the maturity to compute the adjusted age accordingly to the temperature history. The parameters timeHistory and temperatureHistory must be chronologically ordered and paired

        Parameters:
            t: float
                Age which is to be adjusted
            timeHistory: list of floats
                List with consecutive time duration, in days, associated to each isothermal temperature in temperatureHistory, until the desired age to be adjusted.
                Ex: if we want to compute the adjusted age of a material with 30 daus, that initially stayed 20 days under 25C and then 10 days under 30C, then timeHistory = [20,10]
            temperatureHistory: list of floats
                List with the value of isothermal temperatures, in Celsius, during each interval of time 
                Ex: In the previous example at timeHistory description, the respective timeHistory would be [25,30]
        
        Returns:
            adjustedAge: float
                Adjusted time in days
        '''
        #Locate the current time in the temperature history
        temperatureHistoryNew,timeHistoryNew=self.locateTimeInTemperatureHistory(t, temperatureHistory, timeHistory)
        adjustedAge = 0
        for Δti, T_Δti in zip (timeHistoryNew, temperatureHistoryNew):
            #Eq. 5.1-85
            adjustedAge = adjustedAge + Δti*np.exp(13.65-4000/(273+T_Δti))
        
        return adjustedAge
 
    def fcm_computeEvolution (self, t):
        '''
        Section 5.1.9.1: Development of strength with time 
        A method that computes fcm at a given age, in days.
        Age needs to be adjusted according to equivalent age concept (Eq. 5.1-85)
        Returns the value of fcm at the requested age.

        Parameters:
            t: int
                Desired age, in days. For correct modelling, it should be given in adjusted age.
        
        Returns:
            fcmAtAge: float
                Value of fcm at age in MPa
        '''
        #Compute βcc - Eq. 5.1-51
        βcc=self.βcc_compute(t)
        fcmAtAge =  βcc*self.fcm

        return fcmAtAge
        
    def Eci_computeEvolution (self, t):
        '''
        Section 5.1.9.3: Development of modulus of elasticity with time 
        A method that computes Eci at a given age, in days.
        Age needs to be adjusted according to equivalent age concept (Eq. 5.1-85)
        Returns the value of Eci at the requested age.

        Parameters:
            age: int
                Desired age, in days. For correct modelling, it should be given in adjusted age.
        
        Returns:
            EciAtAge: float
                Value of Eci at age in MPa
        '''
        if self.modulusCalculatedPerCode is True:
            #Compute βcc - Eq. 5.1-56
            βcc = self.βcc_compute(t)
            #Compute βe - Eq. 5.1-57
            βe = βcc**0.5
            #Compute elastic modulus evolution
            EciAtAge = βe*self.Eci
        else:
            #Because the three-parameter model is adjusted in terms of hours, here we must 
            #conver "t", which is given in days, to hours
            EciAtAge = self.α_ThreeParModulus*np.exp((-self.𝜏_ThreeParModulus/(t*24))**self.β_ThreeParModulus)

        return EciAtAge

    #Methods to compute code parameters
    def βcc_compute(self,t):
        '''
        Section 5.1.9.1: Development of strength with time
        Computes the parameter βcc, according to Eq. 5.1-51

        Parameters:
            t: int 
                Age in which βcc is calcualted, in days
        
        Returns:
            βcc: float
                Value of parameter accordingly to Eq. 5.1-51
        '''
        #Compute βcc - Eq. 5.1-51
        βcc=np.exp(self.s*(1-(28/t)**0.5))
        return βcc

    def locateTimeInTemperatureHistory(self, t, temperatureHistory, timeHistory):
        #Locate t inside timeHistory
        currentAge=0
        if timeHistory is None:
            print("You need to inform a timeHistory when using ϕ_compute in order to allow computing the adjusted age!")
            exit()
        if temperatureHistory is None:
            print("You need to inform a temperatureHistory when using ϕ_compute in order to allow computing the adjusted age!")
            exit()
        #Locate current time inside the time history
        for iteration,timeStep in enumerate(timeHistory):
            currentAge=currentAge+timeStep
            if t<currentAge:
                #We have found where t0 lies within timeHistory!
                break
        
        #Check if t is larger than the the sum of all time steps in timeHistory, which means it is outside timeHistory scope
        if t>sum(timeHistory):
            print("timeHistory associated to temperatureHistory does not encompass the time t=",t," days. Properly assess the temperature history time scope and the time scope involved in creep analysis, such as in loading age and instants in which creep is being computed, so to ensure all scopes match.")
            exit()

        #Build a vector containing the time and temperature history of loading, and store in timeHistoryLoading and temperatureHistoryLoading
        timeHistoryNew=timeHistory[0:iteration]
        if (t-timeHistory[iteration]<0):
            timeHistoryNew.append(t)
        else:
            timeHistoryNew.append(t-timeHistory[iteration])
        temperatureHistoryNew=temperatureHistory[0:iteration+1]

        return temperatureHistoryNew, timeHistoryNew
    
def readCSVFile():
    import csv
    with open("./costKelvinChains.csv", 'r') as file:
        csvreader = csv.reader(file)
        for iteration1, row in enumerate(csvreader):
            if iteration1 == 0:
                costCompliance=[[] for element in row]
            for iteration2, item in enumerate(row):
                costCompliance[iteration2].append(float(item))
    return costCompliance

if __name__ == '__main__':
    whichAnalysisToPerform = 1 #0: validate with DIANA; 1: validate with code values (Table 5.1-10); 2: validate with COST values
    
    if whichAnalysisToPerform == 0: #Validate with DIANA
        #Read all compliances obtained in Diana from Excel file
        import pandas as pd
        creepDiana = [pd.DataFrame() for element in range(0,24)] #To store the data
        plotGroups = [[2,"Ambient temperature"], [5, "Concrete class"], [9, "Aggregate type"], [12, "Cement type"], [15, "Relative humidity"], [18, "Notional size"], [23, "Age of loading"]] #To organize the plots to be generated
        dividePlotsInGroups = [item[0] for item in plotGroups] #To allow dividing the plots
        legendName = ["20C","40C","60C","C25","C35","C45","Sandstone","Quartzite","Limestone","Basalt","32.5N","32.5R","42.5S","50%","80%","100%","50 mm","150 mm","600 mm","1 day","7 day","28 day","90 day","365 day"] #To allow building legends
        #Read the Excel file
        for element in range(0,24):
            creepDiana[element]=pd.read_excel('2023080301-VALIDATION_DIANA-MC2010.xlsx',sheet_name=str(element+1).zfill(2), header=1, index_col=0)

        #Compute all creep curves with this code
        creepMC2010concrete = [np.zeros((101)) for element in range(0,24)]
        creepMC2010Configurations =   [[375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                       
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C25', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C35', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C45', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'quartzite'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'limestone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'basalt'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5R', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5R', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 50, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 100, 'sandstone'],
                                    [375*100, 3*500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500/4, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone']]
        ϕ_computeConfiguration =   [[28,[40000],[20]],
                                    [28,[40000],[40]],
                                    [28,[40000],[60]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [28,[40000],[20]],
                                    [5,[40000],[20]],
                                    [7,[40000],[20]],
                                    [28,[40000],[20]],
                                    [90,[40000],[20]],
                                    [365,[40000],[20]]]

        timeSpan = creepDiana[0].iloc[:,0] #Use the same instants of Excel
        for element in range(0,24):
            #Construct the concrete
            Ac=creepMC2010Configurations[element][0]
            u=creepMC2010Configurations[element][1]
            concreteClass=creepMC2010Configurations[element][2]
            cementType=creepMC2010Configurations[element][3]
            relativeHumidity=creepMC2010Configurations[element][4]
            aggregateType=creepMC2010Configurations[element][5]
            EC2concrete =  creepMC2010(Ac, u, concreteClass, cementType, relativeHumidity, aggregateType=aggregateType)
            #Compute the creep coeficient
            t0=ϕ_computeConfiguration[element][0]
            timeHistory=ϕ_computeConfiguration[element][1]
            temperatureHistory=ϕ_computeConfiguration[element][2]
            creepMC2010concrete[element]=[EC2concrete.ϕ_compute(t+t0,t0,timeHistory,temperatureHistory) for t in timeSpan]
        
        #Build the plots
        iterator=0
        normalized=False
        for element in range(0,24):
            if normalized is True:
                plt.plot(creepDiana[element].iloc[:,0],creepDiana[element].iloc[:,1]/max(creepMC2010concrete[element]), marker="o", label=legendName[element])
                plt.plot(timeSpan,creepMC2010concrete[element]/max(creepMC2010concrete[element]),marker="x",label=legendName[element]+"-custom")
            else:
                plt.plot(creepDiana[element].iloc[:,0],creepDiana[element].iloc[:,1], marker="o", label=legendName[element])
                plt.plot(timeSpan,creepMC2010concrete[element],marker="x",label=legendName[element]+"-custom")
            if element in dividePlotsInGroups:
                plt.title(plotGroups[iterator][1])
                plt.legend()
                plt.show()
                iterator=iterator+1

    elif whichAnalysisToPerform == 1: #Validade with MC2010 tabulated values
        perimetersToTest = [2500, 833.3333, 208.3333]
        area=250*250
        loading_age = np.array([1,7,28,90,365])
        timeSpan = 50*365
        creepCoefficientTable = np.zeros((3,5))
        
        #Do for RH=50%, ordinary structural concrete
        RH=50
        print("--------------------------------------------------")
        print("Validation according to Table 5.1-10 of MC2010")
        print("--------------------------------------------------")
        print("               Ordinary concrete                  ")
        print("--------------------------------------------------")
        print("                  RH =",RH,"%                     ")
        print("--------------------------------------------------")
        print("           |      Notional size (2Ac/u) [mm]      ")
        print(" Age(days) |--------------------------------------")
        for iteratorHorizontal, u in enumerate(perimetersToTest):
            MC2010Concrete = creepMC2010(area, u, 'C35', '42.5N', RH, modulusCalculatedPerCode=True, aggregateType="limestone", considerTemperatureEffect=False)
            for iteratorVertical, age in enumerate(loading_age):
                t0 = age
                creepCoefficientTable[iteratorHorizontal][iteratorVertical] = MC2010Concrete.ϕ_compute(timeSpan,t0)
        print("           |    ","{:.0f}".format(2*area/perimetersToTest[0]),"   |   ","{:.0f}".format(2*area/perimetersToTest[1]),"   |   ","{:.0f}".format(2*area/perimetersToTest[2]),"    ")
        print("--------------------------------------------------") 
        for index, age in enumerate(loading_age):
            print("   ","{:03d}".format(age),"   |   ","{:.1f}".format(creepCoefficientTable[0][index]),"   |   ","{:.1f}".format(creepCoefficientTable[1][index]),"   |   ","{:.1f}".format(creepCoefficientTable[2][index]))
        
        print("                                                  ")
        #Do for RH=80%, , ordinary structural concrete
        RH=80
        print("--------------------------------------------------")
        print("               Ordinary concrete                  ")
        print("--------------------------------------------------")
        print("                  RH =",RH,"%                     ")
        print("--------------------------------------------------")
        print("           |      Notional size (2Ac/u) [mm]      ")
        print(" Age(days) |--------------------------------------")
        for iteratorHorizontal, u in enumerate(perimetersToTest):
            MC2010Concrete = creepMC2010(area, u, 'C35', '42.5N', RH, modulusCalculatedPerCode=True, aggregateType="limestone", considerTemperatureEffect=False)
            for iteratorVertical, age in enumerate(loading_age):
                t0 = age
                creepCoefficientTable[iteratorHorizontal][iteratorVertical] = MC2010Concrete.ϕ_compute(timeSpan,t0)
        print("           |    ","{:.0f}".format(2*area/perimetersToTest[0]),"   |   ","{:.0f}".format(2*area/perimetersToTest[1]),"   |   ","{:.0f}".format(2*area/perimetersToTest[2]),"    ")
        print("--------------------------------------------------") 
        for index, age in enumerate(loading_age):
            print("   ","{:03d}".format(age),"   |   ","{:.1f}".format(creepCoefficientTable[0][index]),"   |   ","{:.1f}".format(creepCoefficientTable[1][index]),"   |   ","{:.1f}".format(creepCoefficientTable[2][index]))

        print("                                                  ")
        #Do for RH=50%, , high strength concrete
        RH=50
        print("--------------------------------------------------")
        print("             High-strength concrete               ")
        print("--------------------------------------------------")
        print("                  RH =",RH,"%                     ")
        print("--------------------------------------------------")
        print("           |      Notional size (2Ac/u) [mm]      ")
        print(" Age(days) |--------------------------------------")
        for iteratorHorizontal, u in enumerate(perimetersToTest):
            MC2010Concrete = creepMC2010(area, u, 'C80', '42.5N', RH, modulusCalculatedPerCode=True, aggregateType="limestone", considerTemperatureEffect=False)
            for iteratorVertical, age in enumerate(loading_age):
                t0 = age
                creepCoefficientTable[iteratorHorizontal][iteratorVertical] = MC2010Concrete.ϕ_compute(timeSpan,t0)
        print("           |    ","{:.0f}".format(2*area/perimetersToTest[0]),"   |   ","{:.0f}".format(2*area/perimetersToTest[1]),"   |   ","{:.0f}".format(2*area/perimetersToTest[2]),"    ")
        print("--------------------------------------------------") 
        for index, age in enumerate(loading_age):
            print("   ","{:03d}".format(age),"   |   ","{:.1f}".format(creepCoefficientTable[0][index]),"   |   ","{:.1f}".format(creepCoefficientTable[1][index]),"   |   ","{:.1f}".format(creepCoefficientTable[2][index]))

        print("                                                  ")
        #Do for RH=80%, , high strength concrete
        RH=80
        print("--------------------------------------------------")
        print("             High-strength concrete               ")
        print("--------------------------------------------------")
        print("                  RH =",RH,"%                     ")
        print("--------------------------------------------------")
        print("           |      Notional size (2Ac/u) [mm]      ")
        print(" Age(days) |--------------------------------------")
        for iteratorHorizontal, u in enumerate(perimetersToTest):
            MC2010Concrete = creepMC2010(area, u, 'C80', '42.5N', RH, modulusCalculatedPerCode=True, aggregateType="limestone", considerTemperatureEffect=False)
            for iteratorVertical, age in enumerate(loading_age):
                t0 = age
                creepCoefficientTable[iteratorHorizontal][iteratorVertical] = MC2010Concrete.ϕ_compute(timeSpan,t0)
        print("           |    ","{:.0f}".format(2*area/perimetersToTest[0]),"   |   ","{:.0f}".format(2*area/perimetersToTest[1]),"   |   ","{:.0f}".format(2*area/perimetersToTest[2]),"    ")
        print("--------------------------------------------------") 
        for index, age in enumerate(loading_age):
            print("   ","{:03d}".format(age),"   |   ","{:.1f}".format(creepCoefficientTable[0][index]),"   |   ","{:.1f}".format(creepCoefficientTable[1][index]),"   |   ","{:.1f}".format(creepCoefficientTable[2][index]))

    elif whichAnalysisToPerform == 2:
        costConcrete = creepMC2010(250*250, 2*500, 'C50', '42.5N', 60, modulusCalculatedPerCode=False, considerTemperatureEffect=False, α_ThreeParModulus=36000, β_ThreeParModulus=1, 𝜏_ThreeParModulus=28)
        #Read compliance
        costCompliance = np.array(readCSVFile())

        loading_age = np.array([0.2,0.3,0.4,0.5,0.7,0.85,1,1.5,2,3,5,10,15,20,25,30,100])
        for age in loading_age:
            timeSpan = np.arange(age,16+age,0.01)
            t0 = age
            '''
            if costConcrete.modulusCalculatedPerCode is False:
                compliance = [(1/costConcrete.Eci_computeEvolution(t0))+(costConcrete.ϕ_compute(t,t0))/costConcrete.Eci_computeEvolution(t0) for t in timeSpan]
            else:
            '''
            compliance = [(1/costConcrete.Eci_computeEvolution(t0))+(costConcrete.ϕ_compute(t,t0))/costConcrete.Eci_computeEvolution(t0) for t in timeSpan]
            xaxis=[0]+[value-t0 for value in timeSpan]
            #yaxis=[0]+compliance
            yaxis=[0]+compliance
            plt.plot(xaxis,yaxis,label=t0)
            #plt.xscale("log")

        for iterator, dataSet in enumerate(costCompliance[1:]):
            plt.scatter(costCompliance[0],dataSet, label=str(loading_age[iterator])+"-COST", s=5)
        
        plt.legend()
        plt.show() 

    