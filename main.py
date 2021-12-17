# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 15:05:16 2021

@author: josem
"""
from dataReading import ReadingData 
from classSelect import switch
from interpolation import InterpolationFBR,ProcessInterpolation
from flowbetweenplates import FlowPlates,ProcessPlates
from cylinder import FlowCylinder,ProcessCylinder
from plateTemperature import Plate2d
import matplotlib.pyplot as plt
import sys
from time import time
tiempo_inicial = time() 
datas = ReadingData('reading.txt')
typeFunctionRadial,typeProblem,m,n=datas.informationBasic()
q =datas.informationTps()
with switch(typeProblem) as sw:
     if sw.case(1,True):
         beta=datas.informationMulticuadric()
         lowerLimit,upperLimit=datas.informationPlates()
         interpolation=InterpolationFBR(m,n,lowerLimit,upperLimit)
         x,xc=interpolation.generatePoints(lowerLimit,upperLimit)
         pointsMonosCenters=interpolation.functionPointsMonosCenters(x,xc)
         r=interpolation.functionR(x,xc)
         MQparameter=interpolation.functionMulticuadricParameter(x,beta)
         process=ProcessInterpolation(r,m,q,MQparameter,typeFunctionRadial,x)
         matrix=process.matrixA()
         realPolynomialFunction,derivativeFunctionPolynomialReal,secondDerivativeFunctionPolynomicsReal,realOscillatingFunction,derivativeOscillatingFunctionReal,secondDerivativeOscillatingFunctionActual,thirdDerivativeFunctionPolynomicReal,thirdDerivativeFunctionOscilant,quarterDerivativedFunctionOscilant,quarterDerivativedFunctionPolynomicReal=process.getRealFunctions()
         g, derivateG,secondDerivateG,thirdDerivatedG,quarter=process.functionG(pointsMonosCenters,typeFunctionRadial)
         pesosWPolynomial = process.pesosW(matrix,realPolynomialFunction)
         pesosWOscilants=process.pesosW(matrix,realOscillatingFunction)
         functionNumerica,derivateNumerica,secondDerivateNumerica,thirdDerivater,quarterDerivatived=process.numericalApproximations(g,derivateG,secondDerivateG,pesosWPolynomial,thirdDerivatedG,quarter)
         functionNumerica1,derivateNumerica1,secondDerivateNumerica1,thirdDerivater1,quarterDerivatived1=process.numericalApproximations(g,derivateG,secondDerivateG,pesosWOscilants,thirdDerivatedG,quarter)
         process.dataPresentation(functionNumerica,lowerLimit,upperLimit,realPolynomialFunction,derivateNumerica,derivativeFunctionPolynomialReal,secondDerivateNumerica,secondDerivativeFunctionPolynomicsReal,functionNumerica1,realOscillatingFunction,derivateNumerica1,derivativeOscillatingFunctionReal,secondDerivateNumerica1,secondDerivativeOscillatingFunctionActual,thirdDerivater,thirdDerivativeFunctionPolynomicReal,thirdDerivater1,thirdDerivativeFunctionOscilant,quarterDerivatived1,quarterDerivativedFunctionOscilant)
     if sw.case(2,True):
         beta=datas.informationMulticuadric()
         lowerLimit,upperLimit=datas.informationPlates()
         plates=FlowPlates(m,n,lowerLimit,upperLimit)
         x,xc=plates.generatePoints(lowerLimit,upperLimit)
         pointsMonosCenters=plates.functionPointsMonosCenters(x,xc)
         r=plates.functionR(x,xc)
         MQparameter=plates.functionMulticuadricParameter(x,beta)
         upperSpeed,lowerSpeed,pressureGradient,u,radioCylinder=datas.informationCylinderOrPlates()
         process=ProcessPlates(r,m,q,MQparameter,typeFunctionRadial,x,upperSpeed,lowerSpeed,pressureGradient,u)
         matrix=process.matrixA(pointsMonosCenters)
         upperSpeed,lowerSpeed,pressureGradient,u,radioCylinder=datas.informationCylinderOrPlates()
         y=process.vectorY()
         g, derivateG,secondDerivateG,thirdDerivatedG,quarter=process.functionG(pointsMonosCenters,typeFunctionRadial)
         pesosW = process.pesosW(matrix,y)
         functionNumerica,derivateNumerica,secondDerivateNumerica,thirdDerivate,quarterDerivatived=process.numericalApproximations(g,derivateG,secondDerivateG,pesosW,thirdDerivatedG,thirdDerivatedG)
         k=x
         vx,k =process.function(y,u,upperLimit)
         process.showDatas(functionNumerica,k,vx)
     if sw.case(3,True):
         beta=datas.informationMulticuadric()
         upperSpeed,lowerSpeed,pressureGradient,u,ra5dioCylinder=datas.informationCylinderOrPlates()
         cylinder=FlowCylinder(m,n,-radioCylinder,radioCylinder)
         x,xc=cylinder.generatePoints(-radioCylinder,radioCylinder)
         pointsMonosCenters=cylinder.functionPointsMonosCenters(x,xc)
         r=cylinder.functionR(x,xc)
         MQparameter=cylinder.functionMulticuadricParameter(x,beta)
         process=  ProcessCylinder(r, m, q, MQparameter, typeFunctionRadial, x, upperSpeed, lowerSpeed, pressureGradient, u,radioCylinder,n,xc,b=1)
         matrix=process.matrixA(pointsMonosCenters,r)
         g, derivateG,secondDerivateG,thirdDerivatedG,quarter=process.functionG(pointsMonosCenters,typeFunctionRadial)
         y=process.vectorY()
         pesosW=process.pesosW(matrix,y)
         functionNumerica,derivateNumerica,secondDerivateNumerica,thirdDerivate,quarterDerivatived=process.numericalApproximations(g,derivateG,secondDerivateG,pesosW,thirdDerivatedG,thirdDerivatedG)
         #k=x
         vx,k =process.function()
         process.showDatas(functionNumerica,k,vx)
         process.showCylinder(functionNumerica,k,vx)
     if sw.case(4,True):
         temp,lx,ly,b=datas.dataTemperatureYplate()
         plate = Plate2d(m, n, lx, ly, temp, b ,q)
         x,xc,y,yc,arrayNodos,format1=plate.generatorPointsPlate()
         matrix,temperature=plate.matrixA(arrayNodos, format1)
         pesosW=plate.pesosW(matrix, temperature)
         g, derivateG,secondDerivateG,erre=plate.functionGTPS(arrayNodos)
         functionNumerica,derivadaNumerica,segundaDerivadaNumerica=plate.numericalApproximations( g, derivateG,secondDerivateG,pesosW)
         plate.postProceso(arrayNodos[:,0],arrayNodos[:,1],functionNumerica,'Mapa de calor aproximado (°C)')
         arrayTemperature=plate.analyticalTemperature(arrayNodos,temp['TemperaturaInferior'])
         print(temp['TemperaturaInferior'])
         plt.figure()
         plate.postProceso(arrayNodos[:,0],arrayNodos[:,1],arrayTemperature,'Mapa de calor analitico (°C)')
tiempo_final = time() 
tiempo_ejecucion = tiempo_final - tiempo_inicial
print('El tiempo de ejecucion fue:',tiempo_ejecucion)
sys.exit()