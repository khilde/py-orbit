import random
import math
from orbit_utils import Function

class probabilityStrippingWithChicane:

    def __init__ (self,functionx,functiony,n,maxXValue,gamma,beta,computeStripping=True):
        self.name="probabilityStrippingWithChicane"
        self.magneticFieldx=functionx
        self.magneticFieldy=functiony
        self.maxXValue=maxXValue
        self.n=n
        self.gamma=gamma
        self.beta=beta
        self.accumlatedSum= None
        self.CDF= None
        self.computeStripping=computeStripping
        
        #tracking functions for a charged particle that is stripped to neutral, they are implicitly scaled by rigidity (ie values *rigidity which means to get value need to divide by rigidity)
        self.deltaxp_rigidity= None  
        self.deltax_rigidity= None 
        self.deltayp_rigidity= None  
        self.deltay_rigidity= None 
        
        #tracking functions for a neutral particle that is stripped to have charge, they are implicitly scaled by rigidity (ie values *rigidity which means to get value need to divide by rigidity)
        self.deltaxp_m_rigidity= None  
        self.deltax_m_rigidity= None          
        self.deltayp_m_rigidity= None  
        self.deltay_m_rigidity= None   
        
        self.InverseFunction= None     
    def tau(self,x):
    	    A1=2.47e-6
    	    A2=4.49e9
    	    c=299792458   
    	    if self.magneticFieldx.getY(x) ==0 and self.magneticFieldy.getY(x) ==0:
    	    	    #tau is infinite. ie it cant be stripped
    	    	    return -1
    	    else : 
    	    	    totalMagneticField=math.sqrt(self.magneticFieldx.getY(x)*self.magneticFieldx.getY(x)+self.magneticFieldy.getY(x)*self.magneticFieldy.getY(x))
    	    	    if A2/self.gamma/self.beta/c/totalMagneticField >200:
    	    	    	    #tau is HUGE so assume it is infinite 
    	    	    	    #forumla is tau=A1/C*exp(A2/C), C=gamma*beta*c*B , as B->0 tau->inf, and is dominated by exp(A2/C). 
    	    	    	    #A2/C>200 implies exp(A2/C)>exp(200) and A1/C>(A1/A2*200) which implies A1/C*exp(A2/C) >A1/A2*200exp(200)=7.95e73 (which is way older than the age of the universe so it should be save to say tau is infinite
    	    	    	    return -1
    	    	    else:    	    	   
			    #print "totalMagneticField= ",totalMagneticField
			    #print "self.magneticFieldx.getY(x)= ",self.magneticFieldx.getY(x)
			    #print "self.magneticFieldy.getY(x)= ",self.magneticFieldy.getY(x)
			    return A1/self.gamma/self.beta/c/totalMagneticField*math.exp(A2/self.gamma/self.beta/c/totalMagneticField)   	    
    def computeFunctions(self):
    	    c=299792458  
    	    stepSize=self.maxXValue/self.n
    	    theSum=0;
    	    theSumDeltaxp=0
    	    theSumDeltax=0
    	    theSumDeltaxp_m=0
    	    theSumDeltax_m=0    
	    
    	    theSumDeltayp=0
    	    theSumDeltay=0
    	    theSumDeltayp_m=0
    	    theSumDeltay_m=0  
	    
    	    self.accumlatedSum= Function()
    	    self.CDF= Function()    
    	    self.deltaxp_rigidity=Function()
    	    self.deltax_rigidity=Function()
    	    self.deltaxp_m_rigidity=Function()
    	    self.deltax_m_rigidity=Function()    
	    
    	    self.deltayp_rigidity=Function()
    	    self.deltay_rigidity=Function()
    	    self.deltayp_m_rigidity=Function()
    	    self.deltay_m_rigidity=Function()    	    
	    
    	    self.InverseFunction=Function()
   	    #normalizeValue=1-self.probabilityOfSurvival(self.maxXValue)
   	    #print "normalizeValue= %d" %normalizeValue
   	    canBeStripped=False
    	    for i in range(self.n):
    	    	    x = stepSize*i
    	    	    y=1
    	    	    if self.computeStripping:
    	    	    	    y = self.tau(x)
    	    	    #if tau is not infinite
		    if y>0 and self.computeStripping:
		    	    theSum=theSum+stepSize/self.gamma/self.beta/c/y
		    temp_theSumDeltaxp=theSumDeltaxp
		    theSumDeltaxp=theSumDeltaxp+self.magneticFieldy.getY(x)*stepSize 
		    theSumDeltax=theSumDeltax+(theSumDeltaxp+temp_theSumDeltaxp)/2.*stepSize
		    
		    self.deltaxp_rigidity.add(x,theSumDeltaxp)
		    self.deltax_rigidity.add(x,theSumDeltax)
		    
		    theSumDeltaxp_m=0
		    theSumDeltax_m=0
		    if self.computeStripping:
			    for j in range(i,self.n):
				    x_m = stepSize*j
				    temp_theSumDeltaxp_m=theSumDeltaxp_m
				    theSumDeltaxp_m=theSumDeltaxp_m+self.magneticFieldy.getY(x_m)*stepSize
				    theSumDeltax_m=theSumDeltax_m+(theSumDeltaxp_m+temp_theSumDeltaxp_m)/2.*stepSize
			    self.deltaxp_m_rigidity.add(x,theSumDeltaxp_m)
			    self.deltax_m_rigidity.add(x,theSumDeltax_m)
		    
		    temp_theSumDeltayp=theSumDeltayp
		    theSumDeltayp=theSumDeltayp+self.magneticFieldx.getY(x)*stepSize
		    theSumDeltay=theSumDeltay+(theSumDeltayp+temp_theSumDeltayp)/2.*stepSize
		    
		    self.deltayp_rigidity.add(x,theSumDeltayp)
		    self.deltay_rigidity.add(x,theSumDeltay)
		    
		    theSumDeltayp_m=0
		    theSumDeltay_m=0  
		    if self.computeStripping:
			    for j in range(i,self.n):
				    y_m = stepSize*j
				    temp_theSumDeltayp_m=theSumDeltayp_m
				    theSumDeltayp_m=theSumDeltayp_m+self.magneticFieldx.getY(y_m)*stepSize
				    theSumDeltay_m=theSumDeltay_m+(theSumDeltayp_m+temp_theSumDeltayp_m)/2.*stepSize
			    self.deltayp_m_rigidity.add(x,theSumDeltayp_m)
			    self.deltay_m_rigidity.add(x,theSumDeltay_m)
		    if self.computeStripping and y>0:
		    	    #y<0 means sum hasnt change and inverting with fail
			    self.accumlatedSum.add(x,theSum)
			    #3e-7 gives 1-exp(-15)=exp(-3e-7)
			    if theSum<15 and theSum>3e-7: 
			    	    canBeStripped=True
				    self.CDF.add(x,(1-math.exp(-theSum)))
	    if canBeStripped:			    
		    wasSuccess=False
		    if self.computeStripping:
			    wasSuccess=self.CDF.setInverse(self.InverseFunction)
			    if wasSuccess is 1:
				    print "successfully Inverted"
			    else:
				    print "failed to Invert"
	    else:
	    	    #it cant be stripped so set CDF at end to something bigger than 10 so later when check to see if a particle is stripped it will fail (0<=random1<1 will never be bigger than 10)
	    	    self.CDF.add(self.maxXValue,10)				    
    def getaccumlatedSum(self):
    	    return self.accumlatedSum
    def getCDF(self):
    	    return self.CDF
    def getdeltaxp_rigidity(self):
    	    return self.deltaxp_rigidity
    def getdeltax_rigidity(self):
    	    return self.deltax_rigidity    	
    def getdeltaxp_m_rigidity(self):
    	    return self.deltaxp_m_rigidity
    def getdeltax_m_rigidity(self):
    	    return self.deltax_m_rigidity    
 	    
    def getdeltayp_rigidity(self):
    	    return self.deltayp_rigidity
    def getdeltay_rigidity(self):
    	    return self.deltay_rigidity    	
    def getdeltayp_m_rigidity(self):
    	    return self.deltayp_m_rigidity
    def getdeltay_m_rigidity(self):
    	    return self.deltay_m_rigidity     	    
    def getInverseFunction(self):
    	    return self.InverseFunction
    def getLength(self,atX):
        return (self.NotNormalizedFunction.getY(atX))
        
