# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 10:22:01 2021
@author: josem
"""

class ReadingData():
    def __init__(self,l):
        self.l=l
    def informationBasic(self):
        f = open(self.l,'r')
        lines=f.readlines()
        typeFunctionRadial =int(lines[1])
    #    print(typeFunctionRadial)
        typeProblem = int(lines[3])
    #    print(typeProblem)
        m = int(lines[5])
    #    print(m)
        n=int(lines[7])
    #    print(n)
        return typeFunctionRadial,typeProblem,m,n
    #tipoFuncionBaseRadial,tipoDeProblema,m,n = informationBasic('reading.txt')
    def informationMulticuadric(self):
        f = open(self.l,'r')
        lines=f.readlines()
        beta=float(lines[9])
        return beta
    #bera=informationMulticuadric('reading.txt')
    #print(bera)
    def informationTps(self):
        f = open(self.l,'r')
        lines=f.readlines()
        a=float(lines[11])
        float(a)
        return a
    #are =informationTps('reading.txt')
    #print(are)
    def informationPlates(self):
        f = open(self.l,'r')
        lines=f.readlines()
        lowerLimit = float(lines[13])
    #    print('Limite inferior',lowerLimit)
        upperLimit =float(lines[15])
    #    print('limite Superior', upperLimit)
        return lowerLimit,upperLimit
    #d,dd =informationPlates('reading.txt')
    def informationCylinderOrPlates(self):
        f = open(self.l,'r')
        lines=f.readlines()
        upperSpeed=float(lines[17])
    #    print('velocidadalara: ',upperSpeed)
        lowerSpeed=float(lines[19])
    #    print('velocidadbaja: ', lowerSpeed)
        pressureGradient=float(lines[21])
        u =float(lines[23])
    #    print("viscocidad",u)
        radioCylinder=float(lines[25])
    #    print(radioCylinder)
        return upperSpeed,lowerSpeed,pressureGradient,u,radioCylinder
    #upperSpeed,lowerSpeed,pressureGradient,u,radioCylinder=informationCylinderOrPlates('reading.txt')
    def dataTemperatureYplate(self):
        f = open(self.l,'r')
        lines=f.readlines() 
        b = float(lines[27])
        temperature={}
        temperature['TemperaturaInferior']=float(lines[29])
        temperature['TemperaturaSuperior']=float(lines[31])
        temperature['TemperaturaIzquierda']=float(lines[33])
        temperature['TemperaturaDerecha']=float(lines[35])
        lx=int(lines[37])
        ly=int(lines[39])
        return temperature,lx,ly,b
  #  temp,lx,ly,b= dataTemperatureYplate('reading.txt')
   # print(temp['TemperaturaInferior'])
    def dateInformationSphere(self):
        f= open(self.l,'r')
        lines=f.readlines()
        radioSphere=float(lines[41])
        radioCylinder=float(lines[25])
        largeCylinder= float(lines[43])
        return  radioSphere,radioCylinder, largeCylinder
    #R,b,a = dateInformationSphere('reading.txt')
    def datePrueba(self):
        f= open(self.l,'r')
        lines=f.readlines() 
        lx=int(lines[45])
        ly=int(lines[47])
        deltaPx= float(lines[49])
        deltaPy=float(lines[51])
        velocidad={}
        velocidad['vxSuperior']=float(lines[53])
        velocidad['vxInferior']=float(lines[55])
        velocidad['vySuperior']=float(lines[57])
        velocidad['vyInferior']=float(lines[59])
        b= float(lines[27])
        return lx,ly,deltaPx,deltaPy,b,velocidad
l=ReadingData('reading.txt')

        