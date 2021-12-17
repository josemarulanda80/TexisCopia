# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 15:31:55 2021

@author: josem
"""
import numpy as np
import matplotlib.pyplot as plt


class InterpolationFBR:
    def __init__(self,m,n,lowerLimit,upperLimit):
        self.m=m
        self.n=n
    def generatePoints(self,lowerLimit,upperLimit):
        x=np.linspace(lowerLimit,upperLimit, num=self.m)   # generar el intervalo
        xc=np.linspace(lowerLimit,upperLimit, num=self.m)   
        return x,xc
    def functionPointsMonosCenters(self,x,xc):
         pointsMonosCenters=np.zeros((self.n,self.m))
         for i in range (0,self.n):
             pointsMonosCenters[i] =x[i]-xc
         return pointsMonosCenters
    def functionR(self,x,xc):
        #Dimensiono la matriz r
        r=np.zeros((self.n,self.m))
        #recorro las columnas
        for i in range (0,self.m):
            r[i,:]=abs(x[i]-xc)
            # vuelvo diagonal 1 para evitar el logaritmo natural cero
        return r
    def functionMulticuadricParameter(self,x,b):
        distanceMoreSmall =np.zeros(self.m)
        a= np.zeros(self.m)
        MQparameter= np.zeros(self.m) #parametro de la multicuadric
        aux=np.zeros(self.m-1)
        aux2=np.zeros(self.m)
        aux=x[1:self.m].copy()
        aux2=np.append(aux,x[self.m-2])
        distanceMoreSmall=abs(x-aux2)
        a=distanceMoreSmall*b
        MQparameter=a.copy()**2
        return MQparameter

class ProcessInterpolation:
    def __init__(self,r,m,q,parameterMQ,typeFunction,x):
        self.r=r
        self.m=m
        self.q=q
        self.parameterMQ=parameterMQ
        self.typeFunction=typeFunction
        self.x=x

    def matrixA(self):
        matrix=np.zeros((self.m,self.m))
        if self.typeFunction == 1:
            matrix= np.sqrt((self.r**2+self.parameterMQ))
        else:
            matrix= (self.r)**(2*self.q)*np.log(self.r)
            matrix =np.where(self.r==0,0,matrix)
        return matrix
    def getRealFunctions(self):
        realPolynomialFunction = self.x**4+self.x+0.5
        derivativeFunctionPolynomialReal=4*self.x**3+1
        secondDerivativeFunctionPolynomicsReal=12*self.x**2
        thirdDerivativeFunctionPolynomicReal=24*self.x
        quarterDerivativedFunctionPolynomicReal=24
        realOscillatingFunction = 0.02*(12+3*self.x-3.5*self.x**2+7.2*self.x**3)*(1+np.cos(4*np.pi*self.x))*(1+0.8*np.sin(3*np.pi*self.x))
        derivativeOscillatingFunctionReal= 0.02*(21.6*self.x**2 - 7*self.x + 3)*(1 + 0.8*np.sin(3 *np.pi*self.x))*(1 + np.cos(4 *np.pi*self.x)) - 0.251327*(7.2*self.x**3 - 3.5 *self.x**2 + 3*self.x + 12)*(1 + 0.8*np.sin(3 *np.pi* self.x))*np.sin(4 *np.pi*self.x) + 0.150796*(7.2 *self.x**3 - 3.5* self.x**2 + 3*self.x + 12)*np.cos(3 *np.pi*self.x)*(1 + np.cos(4 *np.pi* self.x))
        secondDerivativeOscillatingFunctionActual=0.02*((43.2*self.x-7)*(1+0.8*np.sin(9.42477*self.x))*(1+np.cos(12.56637*self.x))+(7.53982*np.cos(9.42477*self.x)*(1+np.cos(12.56637*self.x))-12.56637*np.sin(12.56637*self.x)*(1+0.8*np.sin(9.42477*self.x)))*(21.6*self.x**2-7*self.x+3))-0.251327*(np.cos(12.56637*self.x)*12.56637*(7.2*self.x**3-3.5*self.x**2+3*self.x+12)*(1+0.8*np.sin(9.42477*self.x))+((21.6*self.x**2-7*self.x+3)*(1+0.8*np.sin(9.42477*self.x))+7.53982*np.cos(9.42477*self.x)*(7.2*self.x**3-3.5*self.x**2+3*self.x+12))*np.sin(12.56637*self.x))+0.150796*(-9.42477*np.sin(9.42477*self.x)*(7.2*self.x**3-3.5*self.x**2+3*self.x+12)*(1+np.cos(12.56637*self.x))+np.cos(9.42477*self.x)*((21.6*self.x**2-7*self.x+3)*(1+np.cos(12.56637*self.x))-12.56637*np.sin(12.56637*self.x)*(7.2*self.x**3-3.5*self.x**2+3*self.x+12))) 
        thirdDerivativeFunctionOscilant=0.06*((21.6*self.x**2)-7*self.x+3)*((-16*np.pi**2)*(1+0.8*np.sin(3*np.pi*self.x))*np.cos(4*np.pi*self.x)-71.0612*np.sin(3*np.pi*self.x)*(1+np.cos(4*np.pi*self.x))-189.496*np.sin(4*np.pi*self.x)*np.cos(3*np.pi*self.x))+0.02*((7.2*self.x**3)-(3.5*self.x**2)+3*self.x+12)*((64*np.pi**3)*(1+0.8*np.sin(3*np.pi*self.x))*np.sin(4*np.pi*self.x)+2678.94*np.sin(3*np.pi*self.x)*np.sin(4*np.pi*self.x)-3571.92*np.cos(3*np.pi*self.x)*np.cos(4*np.pi*self.x)-669.736*np.cos(3*np.pi*self.x)*(1+np.cos(4*np.pi*self.x)))+ 0.864*(1+0.8*np.sin(3*np.pi*self.x))*(1+np.cos(4*np.pi*self.x))+0.06*(43.2*self.x-7)*(7.53982*np.cos(3*np.pi*self.x)*(1+np.cos(4*np.pi*self.x))-4*np.pi*(1+0.8*np.sin(3*np.pi*self.x))*np.sin(4*np.pi*self.x))
        quarterDerivativedFunctionOscilant=0.08*((21.6*self.x**2)-7*self.x+3)*((64*np.pi**3)*(0.8*np.sin(3*np.pi*self.x)+1)*np.sin(4*np.pi*self.x)+2678.94*np.sin(3*np.pi*self.x)*np.sin(4*np.pi*self.x)-3571.92*np.cos(3*np.pi*self.x)*np.cos(4*np.pi*self.x)-669.736*np.cos(3*np.pi*self.x)*(np.cos(4*np.pi*self.x)+1))+0.02*(12+3*self.x-3.5*(self.x**2)+7.2*(self.x**3))*((256*np.pi**4)*(0.8*np.sin(3*np.pi*self.x)+1)*np.cos(4*np.pi*self.x)+67329.2*np.sin(3*np.pi*self.x)*np.cos(4*np.pi*self.x)+6312.11*np.sin(3*np.pi*self.x)*(np.cos(4*np.pi*self.x)+1)+93512.7*np.sin(4*np.pi*self.x)*np.cos(3*np.pi*self.x))+0.12*(43.2*self.x-7)*((-16*np.pi**2)*(0.8*np.sin(3*np.pi*self.x)+1)*np.cos(4*np.pi*self.x)-71.0612*np.sin(3*np.pi*self.x)*(np.cos(4*np.pi*self.x)+1)-189.496*np.sin(4*np.pi*self.x)*np.cos(3*np.pi*self.x))+3.456*(7.53982*np.cos(3*np.pi*self.x)*(np.cos(4*np.pi*self.x)+1)-4*np.pi*(0.8*np.sin(3*np.pi*self.x)+1)*np.sin(4*np.pi*self.x))
        return  realPolynomialFunction,derivativeFunctionPolynomialReal,secondDerivativeFunctionPolynomicsReal,realOscillatingFunction,derivativeOscillatingFunctionReal,secondDerivativeOscillatingFunctionActual,thirdDerivativeFunctionPolynomicReal,thirdDerivativeFunctionOscilant,quarterDerivativedFunctionOscilant,quarterDerivativedFunctionPolynomicReal


    def functionG(self,pointsMonosCenters,typeFunctionRadial):
        if typeFunctionRadial ==1:
            g = np.sqrt((self.r**2+self.parameterMQ))
            derivateG=(((self.r**2)+self.parameterMQ)**(-0.5))*(pointsMonosCenters)
            secondDerivateG=((self.r**2)+(self.parameterMQ)-pointsMonosCenters**2)/((self.r**2+self.parameterMQ)**1.5)
            thirdDerivatedG=0
            quarter=0
        elif  typeFunctionRadial ==2:
            g = (self.r)**(2*self.q)*np.log(self.r)
            g =np.where(self.r==0,0,g)
            derivateG =((self.r)**(2*self.q-1))*(2*self.q*np.log(self.r)+1)*((pointsMonosCenters)/self.r)
            derivateG =np.where(self.r==0,0,derivateG)
            secondDerivateG = self.r**(2.0*(self.q-1))*(4*self.q**2.0*np.log(self.r)+4.0*self.q-2*self.q*np.log(self.r)-1)   
            secondDerivateG=np.where(self.r==0,0,secondDerivateG)     
            thirdDerivatedG= ((2*self.r**(2*self.q-3))*(4*(self.q**3)*np.log(self.r)+6*(self.q**2)*(1-np.log(self.r))+2*self.q*np.log(self.r)-6*self.q+1))*((pointsMonosCenters)/self.r)**3
            thirdDerivatedG=np.where(self.r==0,0,thirdDerivatedG)   
            quarter=(self.r**(2.0*self.q-4))*((32*self.q**3)-(72*self.q**2)+4*(4*self.q**3-12*self.q**2+11*self.q-3)*self.q*np.log(self.r)+44*self.q-6)
            quarter=np.where(self.r==0,0,quarter)   
        else:
            print("Error")
        return g, derivateG,secondDerivateG,thirdDerivatedG,quarter

    def pesosW(self,matrix,y):
        matrixInverse = np.linalg.inv(matrix)
        pesosW =  np.matmul(matrixInverse, y)
        return pesosW
    def numericalApproximations(self,gMQ,derivadagMQ,segundaDerivadagMQ,pesosW,thirdDerivatedMg,quarterDerivatived):
        functionNumerica=np.matmul(gMQ,pesosW)
        derivateNumerica=np.matmul(derivadagMQ,pesosW)
        secondDerivateNumerica=np.matmul(segundaDerivadagMQ,pesosW)
        thirdDerivateNumerica=np.matmul(thirdDerivatedMg,pesosW)
        quarterDerivatived=np.matmul(quarterDerivatived,pesosW)
        return functionNumerica,derivateNumerica,secondDerivateNumerica,thirdDerivateNumerica,quarterDerivatived
    
    def compareGraphicInterpolation(self,functionNumerica,lowerLimit,upperLimit,left,real,interpolated,title,polynomialfunction):
       plt.figure()
       plt.xlabel('Eje X')
       plt.ylabel('Eje Y')
       x1= np.linspace(lowerLimit,upperLimit, num=len(self.x))
       plt.plot(x1,polynomialfunction,color="blue", linewidth=2.5, linestyle="-", label=real)
       plt.plot(self.x,functionNumerica,color="red",  linewidth=2.5, label=interpolated)
       plt.title(title)
       plt.legend(loc= left)
       plt.show()
       return plt.show()
   
    def dataPresentation(self,functionNumerica,lowerLimit,upperLimit,realPolynomialFunction,derivateNumerica,derivativeFunctionPolynomialReal,secondDerivateNumerica,secondDerivativeFunctionPolynomicsReal,functionNumerica1,realOscillatingFunction,derivateNumerica1,derivativeOscillatingFunctionReal,secondDerivateNumerica1,secondDerivativeOscillatingFunctionActual,thirdDerivater,thirdDerivativeFunctionPolynomicReal,thirdDerivater1,thirdDerivativeFunctionOscilant,quarterDerivatived1,quarterDerivativedFunctionOscilant):
        if self.typeFunction==1:
            self.compareGraphicInterpolation(functionNumerica,lowerLimit,upperLimit,'upper left','Función real', 'Función interpolada MQ','Grafica funciones',realPolynomialFunction)
            self.compareGraphicInterpolation(derivateNumerica,lowerLimit,upperLimit,'upper right','Derivada real', 'Derivada interpolada MQ','Grafica derivadas',derivativeFunctionPolynomialReal)
            self.compareGraphicInterpolation(secondDerivateNumerica,lowerLimit,upperLimit,'upper left','Segunda derivada real', 'Segunda derivada interpolada MQ','Grafica segundas derivadas',secondDerivativeFunctionPolynomicsReal)
            self.compareGraphicInterpolation(functionNumerica1,lowerLimit,upperLimit,'upper left','Función real', 'Función interpolada MQ','Grafica funciones',realOscillatingFunction)
            self.compareGraphicInterpolation(derivateNumerica1,lowerLimit,upperLimit,'upper right','Derivada real', 'Derivada interpolada MQ','Grafica derivadas',derivativeOscillatingFunctionReal)
            self.compareGraphicInterpolation(secondDerivateNumerica1,lowerLimit,upperLimit,'upper left','Segunda derivada real', 'Segunda derivada interpolada MQ','Grafica segundas derivadas',secondDerivativeOscillatingFunctionActual)
        elif self.typeFunction ==2:
            self.compareGraphicInterpolation(functionNumerica,lowerLimit,upperLimit,'upper left','Función real', 'Función interpolada TPS','Grafica funciones',realPolynomialFunction)
            self.compareGraphicInterpolation(derivateNumerica,lowerLimit,upperLimit,'upper right','Derivada real', 'Derivada interpolada TPS','Grafica derivadas',derivativeFunctionPolynomialReal)
            self.compareGraphicInterpolation(secondDerivateNumerica,lowerLimit,upperLimit,'upper left','Segunda derivada real', 'Segunda derivada interpolada TPS','Grafica segundas derivadas',secondDerivativeFunctionPolynomicsReal)
            self.compareGraphicInterpolation(thirdDerivater,lowerLimit,upperLimit,'upper left','Tercera derivada real', 'Tercera derivada interpolada TPS','Grafica terceras derivadas',thirdDerivativeFunctionPolynomicReal)
            self.compareGraphicInterpolation(functionNumerica1,lowerLimit,upperLimit,'upper left','Función real', 'Función interpolada TPS','Grafica funciones',realOscillatingFunction)
            self.compareGraphicInterpolation(derivateNumerica1,lowerLimit,upperLimit,'upper right','Derivada real', 'Derivada interpolada TPS','Grafica derivadas',derivativeOscillatingFunctionReal)
            self.compareGraphicInterpolation(secondDerivateNumerica1,lowerLimit,upperLimit,'upper left','Segunda derivada real', 'Segunda derivada interpolada TPS','Grafica segundas derivadas',secondDerivativeOscillatingFunctionActual)
            self.compareGraphicInterpolation(thirdDerivater1,lowerLimit,upperLimit,'upper left','Tercera derivada real', 'Tercera derivada interpolada TPS','Grafica Terceras derivadas',thirdDerivativeFunctionOscilant)
            self.compareGraphicInterpolation(quarterDerivatived1,lowerLimit,upperLimit,'upper left','Cuarta derivada real', 'Cuarta derivada interpolada TPS','Grafica cuartas derivadas',quarterDerivativedFunctionOscilant)
        else:
            print("Error")