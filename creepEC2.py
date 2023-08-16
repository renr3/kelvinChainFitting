import numpy as np
import matplotlib.pyplot as plt
import sys
#Module for computing compliance with MC2010


"""

"""

class creepEC2():
    """
    Define creep-related properties for creep coefficient calculation according to Eurocde 2 (2009
    """
    def __init__(self, Ac: int, u: int, strengthClass: str, cementType: str, Rh: float, modulusCalculatedPerCode:bool=True, considerTemperatureEffect:bool=True, aggregateType:str=None, α_ThreeParModulus:float=None, β_ThreeParModulus:float=None, 𝜏_ThreeParModulus:float=None):
        """
        Class constructor.
        This class defines creep-related properties for creep coefficient calculation according to Eurocde 2 (2009)

        Parameters
        ----------
        Ac (int): Cross sectional area in mm2 (3.1.4 - (5))
        u (int): Cross sectional perimeter in contact with drying, in mm (3.1.4 - (5))
        strengthClass (str): Strength class accordingly to Eurocode 2 (Table 3.1)
        cementType (str): Cement type accordingly to Eurocode 2 (check Eq. 3.2)
        Rh (float): Average relative humidity of the ambient,  in %
        modulusCalculatedPerCode (bool): Inform whether E-modulus should be calculated using EC2 formulation. If "False", parameters of a three-parameter model must be informed
        considerTemperatureEffect (bool): Inform whether temperature effects and type of cement must be taken into account for correcting t0 (age of loading, as per Eq. B.10)
        aggregateType (str): Aggregate type accordingly to types addressed by Eurocode 2. Influences EC2 formulation for E-modulus estimation. Check item (2) of 3.1.3. Options: "basalt", "quartzite", "limestone", "sandstone"
        α_ThreeParModulus (float): The alfa parameter of a three parameter model used to compute E-modulus evolution. Model must give modulus in MPa and time must be in hours. Used if modulusCalculatedPerCode is False.
        β_ThreeParModulus (float): The beta parameter of a three parameter model used to compute E-modulus evolution. Model must give modulus in MPa and time must be in hours. Used if modulusCalculatedPerCode is False.
        𝜏_ThreeParModulus (float): The tau parameter of a three parameter model used to compute E-modulus evolution. Model must give modulus in MPa and time must be in hours. Used if modulusCalculatedPerCode is False.
        """
    
        #Geometrical information about the member
        self.Ac=Ac #Cross sectional area in mm2 (3.1.4 (5))
        self.u=u #3.1.4 (5), cross sectional perimeter in contract with drying in mm
        self.h=2*self.Ac/self.u #Notional size as per 3.1.4 (5)
        
        #Validate cement type
        if cementType not in ["32.5N","32.5R","42.5N","42.5R","52.5N","52.5R"]:
            raise SystemExit("Cement type does not exist. Check Table 5.1-9")
        else:
            self.cementType = cementType

        #Concrete strengths
        #Characteristics compressive strength
        if strengthClass not in ["C12","C16","C20","C25","C30","C35","C40","C45","C50","C55","C60","C70","C80","C90","C100","C110","C120"]:
            raise SystemExit("Strength class does not exist. Check Table 3.1")
        else:
            self.fck=int(strengthClass[-2:])  
        #Mean value of compressive strength      
        self.fcm = self.fck+8 #Eq. 3.1.2 (5)
        
        #Eq. 3.2 for s-parameter
        #Types of S, N, R according to 3.1.2(6)
        if self.cementType == "32.5N":
            #Class S
            self.s=0.38 #Eq. 3.2
            self.α=-1 #Eq. B.9
        elif self.cementType == "32.5R" or self.cementType == "42.5N":
            #Class N
            self.s=0.25 #Eq. 3.2
            self.α=0 #Eq. B.9
        elif self.cementType == "42.5R" or self.cementType == "52.5N" or self.cementType == "52.5R":
            #Class R
            self.s=0.20 #Eq. 3.2
            self.α=1 #Eq. B.9
        else:
            raise SystemExit("Cement type does not exist. Check comments on Eq. 3.2")

        #Calculation of modulus of elasticity
        self.modulusCalculatedPerCode=modulusCalculatedPerCode
        if self.modulusCalculatedPerCode is True:
            #3.1.3 (2): dependency of modulus of elasticity in types of aggregates
            if aggregateType == "basalt":
                self.αE=1.2
            elif aggregateType == "quartzite":
                self.αE=1
            elif aggregateType == "limestone":
                self.αE=0.9
            elif aggregateType == "sandstone":
                self.αE=0.7
            else:
                print("Wrong or not informed type of aggregate. See 3.1.3 (2)")
            self.Ecm = self.αE*22000*(self.fcm/10)**0.3 #Table 3.1, in MPa
            self.Ec = 1.05*self.Ecm #3.1.4 (2)
        else:
            self.α_ThreeParModulus=α_ThreeParModulus
            self.β_ThreeParModulus=β_ThreeParModulus
            self.𝜏_ThreeParModulus=𝜏_ThreeParModulus
            self.Ecm = self.Ecm_computeEvolution(28) #in MPa
            self.Ec = 1.05*self.Ecm #3.1.4 (2)

        #Additional parameters of the code
        self.RH = Rh #In %
        self.considerTemperatureEffect = considerTemperatureEffect #Boolean

    #Methods to compute creep properties
    def ϕ_compute (self, t, t0, timeHistory:list=None, temperatureHistory:list=None):
        """
        Parameters:
            t: float
                Time in which ϕ is to be computed
            t0: float
                Loading age
            timeHistory: list of floats
                List with consecutive time duration, in days, associated to each isothermal temperature in temperatureHistory, until the desired age of calculation.
                Ex: if we want to compute the ϕbc at age of 30 days of a material that initially stayed 20 days under 25C and then 10 days under 30C, then timeHistory = [20,10]
            temperatureHistory: list of floats
                List with the value of isothermal temperatures, in Celsius, during each interval of time 
                Ex: In the previous example at timeHistory description, the respective timeHistory would be [25,30]
        """
        #Maturity corrections and cement type
        #Eq. B.10
        t0_adj = self.t0_adj_compute(t0, timeHistory, temperatureHistory)

        #Creep coefficient
        creepCoefficient = self.ϕ0_compute(t0_adj)*self.βc_compute(t, t0_adj)

        return creepCoefficient

    def ϕ0_compute (self, t0):
        '''
        Method to compute notional creep coefficient, accordingle to Eq. B.2

        Parameters:
            timeHistory: list of floats
                List with consecutive time duration, in days, associated to each isothermal temperature in temperatureHistory, until the desired age of calculation.
                Ex: if we want to compute the ϕbc at age of 30 days of a material that initially stayed 20 days under 25C and then 10 days under 30C, then timeHistory = [20,10]
            temperatureHistory: list of floats
                List with the value of isothermal temperatures, in Celsius, during each interval of time 
                Ex: In the previous example at timeHistory description, the respective timeHistory would be [25,30]
            t0: float
                Loading age, in days.

        Returns:
            ϕ0: float
                Notional creep coefficient Eq. B.2
        '''
        #Eq. B.3
        if self.fcm <= 35:
            ϕrh = 1+(1-self.RH/100)/(0.1*self.h**(1/3)) #Eq. B.3a
        else:
            α1 = (35/self.fcm)**0.7 #Eq. B.8c
            α2 = (35/self.fcm)**0.2 #Eq. B.8c
            ϕrh = (1+α1*(1-self.RH/100)/(0.1*self.h**(1/3)))*α2 #Eq. B.3b
        #Eq. B.4
        βfcm = 16.8/((self.fcm)**0.5)
        #Eq. B.5
        βt0 = 1/(0.1+t0**0.2)

        #Eq. B.2
        ϕ0 = ϕrh*βfcm*βt0

        return ϕ0

    def βc_compute (self, t, t0):
        '''
        Compute the parameter βc that describe the development of creep with time after loading (Eq. B.7)

        Parameters:
            t: float
                Age in which the creep function is being evaluated, in days
            t0: float
                Loading age, in days, already corrected for temperature and cement type effects (Eq. B.10), if applicable.

        Returns:
            βc: float
                βc coefficient according to Eq. B.6
        '''
        #Eq. B.8
        if self.fcm <= 35:
            βh = min(1.5*(1+(0.012*self.RH)**18)*self.h+250, 1500) #Eq. B.8a
        else:
            α3 = (35/self.fcm)**0.5 #Eq. B.8c
            βh = min(1.5*(1+(0.012*self.RH)**18)*self.h+250*α3, 1500*α3) #Eq. B.8b
        #Eq. B.7
        βc = ((t-t0)/(βh+t-t0))**0.3

        return βc

    def t0_adj_compute (self, t0, timeHistory, temperatureHistory):
        '''
        Method to compute the adjusted t0 according to Eq. B.9.
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
        #Eq. B.9
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
        #Locate t0 inside timeHistory
        currentAge=0
        #Make sure we have what it is needed for maturity correction
        if timeHistory is None:
            print("You need to inform a timeHistory when using ϕ_compute in order to allow computing the adjusted age!")
            exit()
        if temperatureHistory is None:
            print("You need to inform a temperatureHistory when using ϕ_compute in order to allow computing the adjusted age!")
            exit()
        for iteration,timeStep in enumerate(timeHistory):
            currentAge=currentAge+timeStep
            if t<currentAge:
                #We have found where t0 lies within timeHistory!
                break
        #Build a vector containing the time and temperature history of loading, and store in timeHistoryLoading and temperatureHistoryLoading
        timeHistoryNew=timeHistory[0:iteration]
        if len(timeHistoryNew) == 0:
            #If it is equal to zero, then the current t is actually within the first time emcompassed by timeHistory. So, timeHistoryNew will just be equal to t
            timeHistoryNew=[t]
        elif (t-timeHistory[iteration]>0):
            #If t is actually above the last iterated timeHistory, given by timeHistory[iteration], this would just happen if t > timeHistory[-1]
            #This means timeHistory does not encompass t, and thus calculation is impossible
            print("@adjustedAge: You need to specify a timeHistory that encompass your largest time span! Exiting the script.")
            exit()
        else:
            #Otherwise, we must add what remains until we achieve t from the last step
            timeHistoryNew.append(t-timeHistory[iteration-1])
        temperatureHistoryNew=temperatureHistory[0:iteration+1]
        adjustedAge = 0
        for Δti, T_Δti in zip (timeHistoryNew, temperatureHistoryNew):
            #Eq. B.10
            adjustedAge = adjustedAge + Δti*np.exp(13.65-4000/(273+T_Δti))
        
        return adjustedAge
 
    def fcm_computeEvolution (self, t):
        '''
        Section 3.1.2 (6): Development of strength with time 
        A method that computes fcm at a given age, in days.
        Age needs to be adjusted according to equivalent age concept
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
        
    def Ecm_computeEvolution (self, t):
        '''
        This method may either use the EC2 approach, or a three-parameter model fitted from experimental results.
        The definition is giving when constructing the class, via the parameter modulusCalculatedPerCode.

        If using the EC2 approach, the following observations are applicable:
            Section 3.1.3(3): Development of modulus of elasticity with time 
            A method that computes Eci at a given age, in days.
            Age needs to be adjusted according to equivalent age concept
            Returns the value of Eci at the requested age.

        If using the three-parameter model, the modulus is calculated at age "t" using the parameters of the model.
        Model must give modulus in MPa and time must be in hours.

        Parameters:
            age: int
                Desired age, in days. For correct modelling, it should be given in adjusted age.
        
        Returns:
            EcmAtAge: float
                Value of Eci at age in MPa
        '''
        if self.modulusCalculatedPerCode is True:
            #Compute fcm at the given age
            fcmt = self.fcm_computeEvolution(t) 
            #Compute elastic modulus evolution
            EcmAtAge = ((fcmt/self.fcm)**0.3)*self.Ecm
        else:
            #Because the three-parameter model is adjusted in terms of hours, here we must 
            #convert "t", which is given in days, to hours
            EcmAtAge = self.α_ThreeParModulus*np.exp((-self.𝜏_ThreeParModulus/(t*24))**self.β_ThreeParModulus)

        return EcmAtAge

    #Methods to compute code parameters
    def βcc_compute(self,t):
        '''
        Section 3.1.2(6): Development of strength with time
        Computes the parameter βcc, according to 3.1.2(6)

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
    whichAnalysisToPerform = 2 #0: validate with DIANA FEA values; 1: validate with COST values; 2:fix stability issues with loading age higher than 100

    if whichAnalysisToPerform == 0:
        #Read all compliances obtained in Diana from Excel file
        import pandas as pd
        creepDiana = [pd.DataFrame() for element in range(0,24)] #To store the data
        plotGroups = [[2,"Ambient temperature"], [5, "Concrete class"], [9, "Aggregate type"], [12, "Cement type"], [15, "Relative humidity"], [18, "Notional size"], [23, "Age of loading"]] #To organize the plots to be generated
        dividePlotsInGroups = [item[0] for item in plotGroups] #To allow dividing the plots
        legendName = ["20C","40C","60C","C25","C35","C45","Sandstone","Quartzite","Limestone","Basalt","N","R","S","50%","80%","100%","50 mm","150 mm","600 mm","1 day","7 day","28 day","90 day","365 day"] #To allow building legends
        #Read the Excel file
        for element in range(0,24):
            creepDiana[element]=pd.read_excel('20230803-VALIDATION_DIANA-EC2.xlsx',sheet_name=str(element+1).zfill(2), header=1, index_col=0)

        #Compute all creep curves with this code
        creepEC2concrete = [np.zeros((101)) for element in range(0,24)]
        creepEC2Configurations =   [[375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C25', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C35', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C45', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'quartzite'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'limestone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'basalt'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5R', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '32.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 50, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 100, 'sandstone'],
                                    [375*100, 3*500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500/4, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone']]
        ϕ_computeConfiguration =   [[28,[36501],[20]],
                                    [28,[36501],[40]],
                                    [28,[36501],[60]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [28,[36501],[20]],
                                    [1,[36501],[20]],
                                    [7,[36501],[20]],
                                    [28,[36501],[20]],
                                    [90,[36501],[20]],
                                    [365,[36501],[20]]]
        '''
        creepEC2Configurations =   [[375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone'],
                                    [375*100, 500, 'C30', '42.5N', 80, 'sandstone']]
        ϕ_computeConfiguration =   [[1,[36501],[20]],
                                    [7,[36501],[20]],
                                    [28,[36501],[20]],
                                    [90,[36501],[20]],
                                    [365,[36501],[20]]]
        '''
        timeSpan = creepDiana[0].iloc[:,0] #Use the same instants of Excel
        for element in range(0,24):
            #Construct the concrete
            Ac=creepEC2Configurations[element][0]
            u=creepEC2Configurations[element][1]
            concreteClass=creepEC2Configurations[element][2]
            cementType=creepEC2Configurations[element][3]
            relativeHumidity=creepEC2Configurations[element][4]
            aggregateType=creepEC2Configurations[element][5]
            EC2concrete =  creepEC2(Ac, u, concreteClass, cementType, relativeHumidity, aggregateType=aggregateType)
            #Compute the creep coeficient
            t0=ϕ_computeConfiguration[element][0]
            timeHistory=ϕ_computeConfiguration[element][1]
            temperatureHistory=ϕ_computeConfiguration[element][2]
            creepEC2concrete[element]=[EC2concrete.ϕ_compute(t+t0,t0,timeHistory,temperatureHistory) for t in timeSpan]
        
        #Build the plots
        iterator=0
        for element in range(0,24):
            plt.plot(creepDiana[element].iloc[:,0],creepDiana[element].iloc[:,1], marker="o", label=legendName[element])
            plt.plot(timeSpan,creepEC2concrete[element],marker="x",label=legendName[element]+"-custom")
            if element in dividePlotsInGroups:
                plt.title(plotGroups[iterator][1])
                plt.legend()
                plt.show()
                iterator=iterator+1
    elif whichAnalysisToPerform == 1:
        #exampleMaterial = creepEC2(400*200, 2*600, 'C30', '32.5N', 100, aggregateType='limestone')
        costConcrete = creepEC2(250*250, 2*500, 'C50', '42.5N', 60, modulusCalculatedPerCode=False, considerTemperatureEffect=False, α_ThreeParModulus=36000, β_ThreeParModulus=1, 𝜏_ThreeParModulus=28)
        
        #Read compliance
        costCompliance = np.array(readCSVFile())

        #loading_age = np.arange(10,90,5)
        initialStrain=[]
        loading_age = np.array([0.2,0.3,0.4,0.5,0.7,0.85,1,1.5,2,3,5,10,15,20,25,30,100])
        for age in loading_age:
            timeSpan = np.arange(age,16+age,0.01)
            t0 = age
            #compliance = [(1/exampleMaterial.Ecm_computeEvolution(t0))+(exampleMaterial.ϕ_compute(t,t0,[100],[20]))/exampleMaterial.Ec for t in timeSpan]
            #compliance = [(1/costConcrete.Ecm_computeEvolution(t0))+(costConcrete.ϕ_compute(t,t0,[100],[20]))/costConcrete.Ecm_computeEvolution(t0) for t in timeSpan]
            fullCompliance = [(1/costConcrete.Ecm_computeEvolution(t0))+(costConcrete.ϕ_compute(t,t0,[200],[20]))/costConcrete.Ecm_computeEvolution(t0) for t in timeSpan]
            creepCoefficient = [costConcrete.ϕ_compute(t,t0,[100],[20]) for t in timeSpan]
            initialStrain.append(1/costConcrete.Ecm_computeEvolution(t0))
            xaxis=[0]+[value-t0 for value in timeSpan]
            #yaxis=[0]+compliance
            yaxis=[0]+fullCompliance
            plt.plot(xaxis,yaxis,label=t0)
        for iterator, dataSet in enumerate(costCompliance[1:]):
            plt.scatter(costCompliance[0],dataSet, label=str(loading_age[iterator])+"-COST", s=5)

        plt.legend()
        plt.show()
    elif whichAnalysisToPerform == 2:
        #Create EC2 compliances
        costConcrete = creepEC2(250*250, 2*500, 'C50', '42.5N', 100, aggregateType='limestone')

        #loading_age = np.arange(10,90,5)
        initialStrain=[]
        loading_age = np.array([28])
        EcmEvolution=[]
        maxTime=720*500
        for age in loading_age:
            timeHistory=[1e6]
            temperatureHistory=[20]
            timeSpan = np.arange(age,maxTime+age,maxTime/1000)
            t0 = age
            fullCompliance = [(1/costConcrete.Ecm_computeEvolution(t0))+(costConcrete.ϕ_compute(t,t0,timeHistory,temperatureHistory))/costConcrete.Ecm_computeEvolution(t0) for t in timeSpan]
            initialStrain.append(1/costConcrete.Ecm_computeEvolution(t0))
            #xaxis=[0]+[value-t0 for value in timeSpan]
            xaxis=[value-t0 for value in timeSpan]
            creepCoefficient = [costConcrete.ϕ_compute(t,t0,timeHistory,temperatureHistory) for t in timeSpan]
            #yaxis=[0]+compliance
            #yaxis=[0]+compliance
            yaxis=creepCoefficient
            EcmEvolution.append(costConcrete.Ecm_computeEvolution(t0))
            plt.plot(xaxis,yaxis,label=t0)
        plt.legend()
        plt.xscale('log')
        plt.xlabel('Days after loading')
        plt.ylabel("Creep coeficient (-)")
        plt.show()

        #Analyse Ecm evolution with loading age
        '''
        plt.plot(loading_age,EcmEvolution, label="Ecm evolution with loading age")
        plt.xlabel("loading age (days)")
        plt.ylabel("Ecm value (MPa)")
        plt.show()
        '''
    print("Fim")