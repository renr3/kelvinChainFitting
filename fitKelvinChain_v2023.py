import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#Bazant and Osman DPL strategy

class kelvinChainModel ():
    #Model to define an aging Kelvin Chain Model

    def __init__(self, ageingComplianceSeries:list):
        """
        Class constructor.
        This class defines a Kelvin Chain model.
        It is able to handle aging and non-aging models.

        Parameters
        ----------
        complianceSeries (list): List of "m" ageing compliances, comprised by "n" datapoints, formated in the following way: 
        [[t'0, list[[t10,t20,t30,...,tn0],[J10,J20,J30,...,Jn0]],
         [t'1, list[[t11,t21,t31,...,tn1],[J11,J21,J31,...,Jn1]],
         [t'2, list[[t12,t22,t32,...,tn2],[J12,J22,J32,...,Jn2]],
         ...
         [t'm, list[[t1m,t2m,t3m,...,tnm],[J1m,J2m,J3m,...,Jnm]]]
        Not all compliances must have the same number of data points.
        The user can also give just one compliance, in which case a non-aging Kelvin chain model will be fitted.
        Units:
            times must be in days
            compliances must be in Pascals
        """
        self.complianceAges=[complianceEntry[0] for complianceEntry in ageingComplianceSeries]
        self.complianceSeries=[complianceEntry[1] for complianceEntry in ageingComplianceSeries]

        #Additional paramters that may be computed along the various Kelvin Chain fitting strategies available
        self.doublePowerLawModel=None
        self.modifiedDoublePowerLawModel=None
        self.dirichletSeriesModel=None

    #Smoothing curves method
    def fitDoublePowerLaw(self, typeOfDPL:str='classical', graphVisualization:bool=False, mrkSz:int=5, plotSlicingFactor:int=1, axisObjectToPlot:object=None):
        """
        This is the method for fitting a DPL, which can be of the form 'classical' form:

        J(t,t')=(1/E0)+(ϕ1/E0)*(t'^(-m))*((t-t')^n)

        Or the 'modified' form, in which the E-modulus evolution is included in the creep law:

        J(t,t')=(1/E0*exp(p/(t**q)))+(ϕ1/E0)*(t'^(-m))*((t-t')^n)

        Many guidelines of the procedure herein implemented are based on the paper "Double power law for basic creep of concrete", by Bazant and Osman from 1976.

        Parameters
        ----------
        typeOfDPL (str): defines which form of the DPL is to be fitted, the 'classical' or the 'modified' form
        graphVisualization (bool): defines if a graphic should be plotted with the results of the fitting
        mrkSz (int): if graph visualization is true, defines the marker size to be used in the curves
        plotSlicingFactor (int): if graph visualization is true, defines the density of points to be considered, equal to one every plotSlicingFactor
        axisObjectToPlot (axis object): if graph visualization is true, you can pass a matplotlib axis object to be used in the plot (useful when you want to plot over an already existing plot)
        """
        #TODO: Check for non-aging compliance

        #INITIAL VERIFICATION
        #Check time length of creep time and ages
        #Explanation: When only short-time creep data with up to one month duration and up to one month age at loading are available, all four material parameters cannot be determined from creep data. In such a case it is necessary to assume exponents "m" and "n", whereupon determination of Eo and ϕ1 from the short-time data is normally possible.
        maxTime=0
        for index,compliance in enumerate(self.complianceSeries):
            if maxTime < max(compliance[0])-self.complianceAges[index]:
                maxTime = max(compliance[0])-self.complianceAges[index]
        if round(maxTime) <= 30:
            print("All creep data series are shorter than 1 month, which may impossibilitate correct determination of all four parameters of DPL by direct fitting")
            print("Consider assuming exponents m and n as equal to 1/3 and 1/8, respectively")
        maxAge=0
        for complianceAge in self.complianceAges:
            if maxAge < complianceAge:
                maxAge = complianceAge
        if round(maxAge) <= 30:
            print("The ageing curves encompass a time span shorter than 1 month, which may impossibilitate correct determination of all four parameters of DPL by direct fitting")
            print("Consider assuming exponents m and n as equal to 1/3 and 1/8, respectively")
        
        #PERFORM FITTING
        #Build the vector data
        t_line = np.array([[self.complianceAges[index] for element in complianceSeries[0]] for index, complianceSeries in enumerate(self.complianceSeries)])
        t=np.array([complianceSeries[0] for complianceSeries in self.complianceSeries])
        ydata=np.array([complianceSeries[1] for complianceSeries in self.complianceSeries])
        
        #Store the conversion factors
        conversionFactor_compliance=ydata.max()
        
        #Normalize the data to try improving fitting
        t_line_norm = t_line.flatten()
        t_norm=t.flatten()
        ydata_norm=(ydata.flatten())/conversionFactor_compliance

        #Fit the correct type of DPL
        if typeOfDPL == 'classical':
            #Perform the fitting with initial guess given by Bazant and Osman
            popt, pcov = curve_fit(self.DPLfunction, (t_norm, t_line_norm), ydata_norm, p0=[0.5, 0.5, 1/3, 1/8], sigma=ydata_norm**10)
            #Store optimizaton results
            self.doublePowerLawModel={'model': 'DPL','E0':popt[0],'ϕ_1':popt[1],'m':popt[2],'n':popt[3],'conversionFactor_compliance': conversionFactor_compliance, 'pcov':pcov}
        elif typeOfDPL == 'modified':
            popt, pcov = curve_fit(self.modifiedDPLfunction, (t_norm, t_line_norm), ydata_norm, p0=[0.5, 0.5, 1/3, 1/8, 1, 1], sigma=ydata_norm**10)
            #Store optimizaton results
            self.modifiedDoublePowerLawModel={'model': 'modifiedDPL','E0':popt[0],'ϕ_1':popt[1],'m':popt[2],'n':popt[3],'p':popt[4],'q':popt[5],'conversionFactor_compliance': conversionFactor_compliance, 'pcov':pcov}

        #Show how good the optimization was
        if graphVisualization is True:
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(self.complianceSeries))))
            for enum, complianceSeries in enumerate(self.complianceSeries):
                c = next(color)
                #plt.plot(complianceSeries[0][::plotSlicingFactor]-self.complianceAges[enum], np.array(complianceSeries[1][::plotSlicingFactor])*(1e6), label="EC2-"+"{:.1f}".format(self.complianceAges[enum])+" days", c=c, marker='o',markersize=mrkSz)
                if typeOfDPL == 'classical':
                    if axisObjectToPlot is None:
                        plt.plot(complianceSeries[0][::plotSlicingFactor]-self.complianceAges[enum],(np.array([self.DPLfunction((t,self.complianceAges[enum]),popt[0],popt[1],popt[2],popt[3]) for t in complianceSeries[0][::plotSlicingFactor]]))*conversionFactor_compliance*(1e6), c=c, marker='x',markersize=mrkSz,linestyle='--', label="DPL-"+"{:.1f}".format(self.complianceAges[enum])+"d")
                    else: 
                        axisObjectToPlot.plot(complianceSeries[0][::plotSlicingFactor]-self.complianceAges[enum],(np.array([self.DPLfunction((t,self.complianceAges[enum]),popt[0],popt[1],popt[2],popt[3]) for t in complianceSeries[0][::plotSlicingFactor]]))*conversionFactor_compliance*(1e6), c=c, marker='x',markersize=mrkSz,linestyle='--', label="DPL-"+"{:.1f}".format(self.complianceAges[enum])+"d")
                elif typeOfDPL == 'modified':
                    if axisObjectToPlot is None:
                        plt.plot(complianceSeries[0][::plotSlicingFactor]-self.complianceAges[enum],(np.array([self.modifiedDPLfunction((t,self.complianceAges[enum]),popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) for t in complianceSeries[0][::plotSlicingFactor]]))*conversionFactor_compliance*(1e6), c=c, marker='x',markersize=mrkSz,linestyle='-', label="DPL-$t_{0}$="+"{:.1f}".format(self.complianceAges[enum])+"d")
                    else: 
                        axisObjectToPlot.plot(complianceSeries[0][::plotSlicingFactor]-self.complianceAges[enum],(np.array([self.modifiedDPLfunction((t,self.complianceAges[enum]),popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) for t in complianceSeries[0][::plotSlicingFactor]]))*conversionFactor_compliance*(1e6), c=c, marker='x',markersize=mrkSz,linestyle='-', label="DPL-$t_{0}$="+"{:.1f}".format(self.complianceAges[enum])+"d")
            if axisObjectToPlot is None:
                plt.xscale('log')
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.xlabel("Days since loading")
                plt.ylabel("Compliance ($10^{-6}$$MPa^{-1}$)")
                plt.grid(which='minor', alpha=0.2)
                plt.grid(which='major', alpha=0.5)
                plt.tight_layout()
                plt.show()
            else:
                #Do Nothing!
                0

    def DPLfunction(self, T, E_0, ϕ_1, m, n):
        """
        This method is to evaluate a classic DPL, given the parameters informed
        """
        t,t_line=T        
        return ((1/(E_0)))+(ϕ_1/E_0)*((np.array(t_line))**(-m))*((np.array(t-t_line))**(n))

    def modifiedDPLfunction(self, T, E_0, ϕ_1, m, n, p, q):
        """
        This is a method is to evaluate a modified DPL, with the constant term possessing an age dependency, given the parameters informed
        """
        t,t_line=T        
        #return (1/(E_0*np.exp(-q/(np.array(t_line))**(p))))+(ϕ_1/(E_0*np.exp(-q/(np.array(t_line))**(p))))*((np.array(t_line))**(-m))*((np.array(t-t_line))**(n))
        return (1/(E_0*np.exp(-q/(np.array(t_line))**(p))))+(ϕ_1/(E_0))*((np.array(t_line))**(-m))*((np.array(t-t_line))**(n))

    def computeDoublePowerLaw(self, creepTimesList: list, loadingAge:float):
        """
        A method to compute the values of the fitted DPL stored in self.doublePowerLawModel. It returns a list containing the points of the DPL in each creep time given in creepTimesList, for a given loadingAge.

        Parameters
        ----------
        creepTimesList (list): A list containing all creep times in which the DPL is to be evluated.
        loadingAge (float): a loading age.
        """
        
        J=np.array([self.DPLfunction((t,loadingAge),self.doublePowerLawModel['E0'],self.doublePowerLawModel['ϕ_1'],self.doublePowerLawModel['m'],self.doublePowerLawModel['n']) for t in creepTimesList])*self.doublePowerLawModel['conversionFactor_compliance']
        
        return J

    def computeModifiedDoublePowerLaw(self, creepTimesList: list, loadingAge:float):
        """
        A method to compute the values of the fitted DPL stored in self.doublePowerLawModel. It returns a list containing the points of the DPL in each creep time given in creepTimesList, for a given loadingAge.

        Parameters
        ----------
        creepTimesList (list): A list containing all creep times in which the DPL is to be evluated.
        loadingAge (float): a loading age.
        """
        
        J=np.array([self.modifiedDPLfunction((t,loadingAge),self.modifiedDoublePowerLawModel['E0'],self.modifiedDoublePowerLawModel['ϕ_1'],self.modifiedDoublePowerLawModel['m'],self.modifiedDoublePowerLawModel['n'],self.modifiedDoublePowerLawModel['p'],self.modifiedDoublePowerLawModel['q']) for t in creepTimesList])*self.modifiedDoublePowerLawModel['conversionFactor_compliance']
        
        return J
    
    #Kelvin chain (Dirichlet series) fitting
    def fitDirichletSeries(self, fittingData: dict, loadingAgesInterval: list=None, creepTimesInterval: list=None, qμ:float=1, retardationTimesRange:list=None, retardationTimesFactor:int=None, graphVisualization:bool=False, mrkSz:int=5, plotSlicingFactor:int=1, axisObjectToPlot:object=None):
        """
        This is the method for fitting a Dirichlet series of the form:

        J(t,t')=(1/E(t'))+Σ{(1/Eμ(t'))*(1-exp(-(t-t')/𝜏μ))}

        Parameters
        ----------
        fittingData (dict): The data in which the fitting will be applied. It is to be given in the form of a dictionary. 
        It can be either:
            -one of the implemented models in this class (DPL, etc) (in such case, the data will be generated from the model)
            -or raw data that shall be formatted in the dictionary as: {'model': 'raw', 't_line': [t'1, t'2, t'3, ...], 't': [[t11,t12,t13,t14,...],[t21,t22,t23,t24,...],...], 'J': [[J11,J12,J13,...], [J21,J22,J23,...],...]}, in which each t'n in 't_line' is related to a list of creep times [tn1,tn2,...] in 't' and a list of compliance values [Jn1,Jn2,Jn3,...] in 'J'. In this last case, loadingAgesInterval and creepTimesInterval are disregarded.
        loadingAgesInterval (list): A list [t_line_min, t_line_max, number], in which t_line_min is the minimum loading age, t_line_max is the maximum loading age, and number is the number of loading ages to be sampled from that interval. Mandatory if the fitting is performed based on a previously fitted model for data smoothing.
        creepTimesInterval (list): A list [t_min, t_max, number], in which t_min is the minimum creep time, t_max is the maximum creep time, and number is the number of creep times to be sampled from that interval. Only used if the fitting is performed based on a previously fitted model for data smoothing.
        qμ (float): The exponent qμ in the time factor of the Dirichlet series, such as (t/𝜏μ)^qμ. Traditionally, qμ=1, but the book "Mathematical Modelling of Creep and Shrinkage" indicates qμ=2/3 may better suit concrete data
        retardationTimesRange (list): a list [𝜏1, 𝜏N] in a way that 𝜏1 is the smallest 𝜏 to consider, and 𝜏N is the largest 𝜏 to consider
        retardationTimesFactor (int): an optional factor that will change the traditional 10^(1/qu) factor for distributing the retardation times to (retardationTimesFactor). Useful if a coarser or finer mesh of retardation times is desired.
        graphVisualization (bool): View the result of fitting graphically.
        mrkSz (int): the marker size in the plot.
        plotSlicingFactor (int): the factor that will be used in the graphVisualization to make the graphs less dense.

        Observation:
        ----------
        Some observations regarding this implementation is given:
        - The paper "Dirichlet series creep function for aging concrete", by Bazant and Wu, is the main guidance and referenced by many other sources when it comes to an algorithm for fitting Dirichlet series. It outlines the mechanics behind how to fit a Dirichlet function to an ageing viscoelastic compliance
        - It saus the solution of Dirichlet fitting series is unique and not oversensitive to inaccuracy of data only if the 𝜏μ-values are specified and choice of 𝜏μ is such that the difference between any two 𝜏μ-values is not too small.
        - 𝜏μ should also not be too large for good accuracy of the representaition.
        - For practicity, the m-values of 𝜏μ (μ=1,...,m) are best chosen as evenly distributed in log 𝜏μ-scale.
            . 𝜏μ=10^{μ-1}*𝜏1 (μ=1,...,m)
            . 𝜏1 is the point where the creep curve when plotted in log (t-t') scale begins to rise
        - The book "Applied Analysis" by Lanczos, indicated also by the above paper, explains in detail about the pitfals of the problem of fitting exponential series to real world data (pages 272 to 280). The main problem of fitting exponential series to experimental measurements is that the parameters of an exponential series are highly sensitive to changes in the data (ill-conditioned), such  (𝜏μthat a single set of data may be sucessfuly represented by various largely different solutions (uniqueness of solution is difficult).
        - In the book "Mathematical Modelling of Creep and Shrinkage" by Bazant, there are the following recommendations:
            . 𝜏1 should be taken as a very small number so to represent the instantaneous part of the compliance function. A recomendation is 𝜏1=10^-9 day.
            . 𝜏μ's must not be spaced too sparsely in the log(t-t') scale and must cover the entire time range of interest
            . the smallest 𝜏μ must be such that 𝜏2 <= 3*𝜏min, in which 𝜏min is the smallest time delay after load application for which the response is of interest
            . also, the smallest 𝜏μ must be smaller than the age of concrete t0 when the structure is first loaded, such as that 𝜏2 <= 0.1t0
            . the largest 𝜏μ, identified as 𝜏N, must be such that 𝜏N >= 0.5*𝜏max
            . while the recomendation is to use the Dirichlet series with the exponent with a ratio of (t/𝜏μ)^qμ, with qμ = 1, in "Mathematical..." it is recommended the use of a qμ=2/3, in which case 𝜏μ=10^(1/qμ)*𝜏_(μ-1)
        - Because of the aforementioned criteria over the retardation times values, the number of elements to be included in the chain is defined by the number of retardation times necessary to cover the entire time range
        """
        #Identify what kind of data is passed and get the fitting data:
        if fittingData['model']=='DPL':
            t_lines_real=np.geomspace(loadingAgesInterval[0],loadingAgesInterval[1],num=loadingAgesInterval[2])
            t_real=[np.geomspace(creepTimesInterval[0],creepTimesInterval[1],num=creepTimesInterval[2])+t_line for t_line in t_lines_real]
            #J_real=np.array([self.computeDoublePowerLaw(t,t_line) for t_line, t in zip(t_lines_real, t_real)])/self.doublePowerLawModel['conversionFactor_compliance']
            conversionFactor=self.doublePowerLawModel['conversionFactor_compliance']
            #conversionFactor=1
            J_real=np.array([self.computeDoublePowerLaw(t,t_line) for t_line, t in zip(t_lines_real, t_real)])/conversionFactor
            #modulusGuess=self.doublePowerLawModel['E0']/conversionFactor
            modulusGuess=self.doublePowerLawModel['E0']
        elif fittingData['model']=='modifiedDPL':
            t_lines_real=np.geomspace(loadingAgesInterval[0],loadingAgesInterval[1],num=loadingAgesInterval[2])
            t_real=[np.geomspace(creepTimesInterval[0],creepTimesInterval[1],num=creepTimesInterval[2])+t_line for t_line in t_lines_real]
            conversionFactor=self.modifiedDoublePowerLawModel['conversionFactor_compliance']
            #conversionFactor=1
            J_real=np.array([self.computeModifiedDoublePowerLaw(t,t_line) for t_line, t in zip(t_lines_real, t_real)])/conversionFactor
            #J_real=np.array([self.computeModifiedDoublePowerLaw(t,t_line) for t_line, t in zip(t_lines_real, t_real)])
            modulusGuess=self.modifiedDoublePowerLawModel['E0']
        elif fittingData['model']=='raw':
            #Raw data has been given
            t_lines_real=fittingData['t_line']
            #Populate loadingAgesInterval variable
            loadingAgesInterval=[min(t_lines_real), max(t_lines_real)]
            t_real=fittingData['t']
            #Populate creepTimesInterval variable
            creepTimesInterval=[min([min(creepTimeSeries-t_lines_real[enum]) for enum, creepTimeSeries in enumerate(t_real)]), max([max(creepTimeSeries-t_lines_real[enum]) for enum, creepTimeSeries in enumerate(t_real)])]
            J_real=fittingData['J']
            conversionFactor=max([max(J_series) for J_series in J_real])
            J_real=np.array(J_real)/conversionFactor
            modulusGuess=1/J_real[0][0]
        else:
            print("The data given for fitting a Dirichlet series is not correctly formatted. See the comments of the method fitDirichletSeries for more details.")
            exit()

        #Now, layout the optimization algorithm (check if weights are needed, but most probably they are not)

        #Also, probably one of the inputs of this function will be the retardation times (they are chosen a priori and we only optimize for the moduluses)
        #Define the retardation times based on the creep time span of interest
        if retardationTimesFactor is None:
            factorRetardationTimes=10**(1/qμ)
        else:
            factorRetardationTimes=retardationTimesFactor
        if retardationTimesRange is None:
            𝜏1=10**(-9) #the first retardation time is to model the instantaneous response
            𝜏2=min([3*creepTimesInterval[0], 0.1*loadingAgesInterval[0]]) #according to recommendations of "Mathematical..." book
            𝜏N_proposed=0.5*creepTimesInterval[1]
            retardationTimes=[𝜏1, 𝜏2]
        else:
            𝜏1=retardationTimesRange[0]
            𝜏2=retardationTimesRange[1]
            𝜏N_proposed=retardationTimesRange[2]
            retardationTimes=[𝜏1,𝜏2]
        #Start building the retardation time array
        
        #Now, build the vector retardationTimes. We will stop when the last retardation time is equal or higher than 𝜏N_proposed
        while retardationTimes[-1] <= 𝜏N_proposed:
            retardationTimes.append(retardationTimes[-1]*factorRetardationTimes)

        #Perform the fitting with initial guess given by Bazant and Osman
        modulusList=[]
        pcovList=[]
        for enum, (t_line, t_creep, J) in enumerate(zip(t_lines_real, t_real, J_real)):
            #Select proper bounds
            if enum == 0:
                selectedBounds=(0,1e15)
                p0_guess=[modulusGuess for 𝜏 in retardationTimes]
            else:
                #selectedBounds=(modulusList[-1],[np.inf for elements in modulusList[-1]])
                selectedBounds=(0,1e15)
                #p0_guess=[modulus for modulus in modulusList[-1]]
                p0_guess=[modulusGuess for 𝜏 in retardationTimes]
                #print(selectedBounds)
            popt, pcov = curve_fit(lambda t, *modulus: self.dirichletFunction(retardationTimes, qμ, t_line, t, modulus), t_creep, J, p0=p0_guess, bounds=selectedBounds, maxfev = 5000)
            modulusList.append(popt)
            pcovList.append(pcov)

        self.dirichletSeriesModel={'model': 'dirichlet','agingModulus': modulusList,'retardationTimes': retardationTimes,'qu': qμ,'t_line':t_lines_real, 'conversionFactor_compliance': conversionFactor, 'pcov': pcovList}

        #Show how good the optimization was
        if graphVisualization is True:
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(t_lines_real))))
            for enum, t_lines in enumerate(t_lines_real):
                c = next(color)
                #plt.plot(t_real[enum]-t_lines, J_real[enum]*1e6, label="data-"+"{:.1f}".format(t_lines)+" days", c=c, marker='o', linestyle='none')
                J_dirichletSeries = self.dirichletFunction(retardationTimes, qμ, t_lines, t_real[enum], modulusList[enum])
                if axisObjectToPlot is None:
                    plt.plot(t_real[enum][::plotSlicingFactor]-t_lines, J_dirichletSeries[::plotSlicingFactor]*conversionFactor*1e6, label="KC-$t_{0}$="+"{:.1f}".format(t_lines)+"d", c=c, marker='v',markersize=mrkSz, linestyle='--')
                else:
                    axisObjectToPlot.plot(t_real[enum][::plotSlicingFactor]-t_lines, J_dirichletSeries[::plotSlicingFactor]*conversionFactor*1e6, label="KC-$t_{0}$="+"{:.1f}".format(t_lines)+"d", c=c, marker='v',markersize=mrkSz, linestyle='--')
                #plt.plot(complianceSeries[0][::plotSlicingFactor]-self.complianceAges[enum],(np.array([self.DPLfunction((t,self.complianceAges[enum]),popt[0],popt[1],popt[2],popt[3]) for t in complianceSeries[0][::plotSlicingFactor]]))*conversionFactor_compliance*(1e6), c=c, marker='x',markersize=mrkSz,linestyle='none', label="DPL-"+"{:.1f}".format(self.complianceAges[enum])+" days")
            if axisObjectToPlot is None:
                plt.xscale('log')
                plt.yscale('log')
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.xlabel("Days since loading")
                plt.ylabel("Compliance ($10^{-6}$$MPa^{-1}$)")
                plt.tight_layout()
                plt.grid(which='minor', alpha=0.2)
                plt.grid(which='major', alpha=0.5)
                plt.show()
            else:
                #Do nothing!
                0
        
    def dirichletFunction(self, retardationTimes, qu, t_line, t, *modulus):
        """
        A method to build a Dirichlet function of the form:

        J(t,t')=(1/E(t'))+Σ{(1/Eμ(t'))*(1-exp([(t'/𝜏μ)^(qu)]-[(t/𝜏μ)^(qu)]))}

        It will online compute J(t,t') for one t' at a time.
        """
        J=0
        for Eμ,𝜏μ in zip(modulus[0],retardationTimes):
            J = J + (1/Eμ)*(1-np.exp(-(((t)**qu)-((t_line)**qu))/𝜏μ))
        
        return J