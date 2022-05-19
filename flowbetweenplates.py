# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 19:20:32 2021

@author: josem
"""
from interpolation import InterpolationFBR,ProcessInterpolation
import matplotlib.pyplot as plt
import numpy as np
# Herencia simple
# class FlowPlates(ProcessInterpolation):
#    pass
class FlowPlates(InterpolationFBR):
    def __init__(self,m,n,lowerLimit,upperLimit):
        InterpolationFBR.__init__(self,m,n,lowerLimit,upperLimit)
class ProcessPlates(ProcessInterpolation):
    def __init__(self,r,m,q,parameterMQ,typeFunction,x,upperSpeed,lowerSpeed,pressureGradient,u):
        ProcessInterpolation.__init__(self,r,m,q,parameterMQ,typeFunction,x)
        self.upperSpeed=upperSpeed
        self.lowerSpeed=lowerSpeed
        self.pressureGradient=pressureGradient
        self.u=u        
    def matrixA(self,pointsMonosCenters):
        matrix=np.zeros((self.m,self.m))
        if self.typeFunction == 1:
            matrix=((self.r**2)+(self.MQparameter)-pointsMonosCenters**2)/((self.r**2+self.MQparameter)**1.5)
            matrix[0,:]= np.sqrt((self.r[0,:]**2+self.parameterMQ))
            matrix[self.m-1,:]= np.sqrt((self.r[self.m-1,:]**2+self.parameterMQ))
        elif self.typeFunction == 2:
            matrix=np.zeros((self.m,self.m))
            matrix[self.m-1,:]= (self.r[self.m-1,:])**(2*self.q)*np.log(self.r[self.m-1,:])
            #matriz[1:m-1,:]= (r[1:m-1])**(2*GradoTps-2)*(2*GradoTps-1)*((2*GradoTps*np.log(r[1:m-1])+1)+2*GradoTps)    
            matrix[1:self.m-1,:]= self.r[1:self.m-1,:]**(2.0*(self.q-1))*(4*self.q**2.0*np.log(self.r[1:self.m-1,:])+4.0*self.q-2*self.q*np.log(self.r[1:self.m-1,:])-1)
            matrix[0,:]= (self.r[0,:])**(2*self.q)*np.log(self.r[0,:])
            matrix=np.where(self.r==0,0,matrix)
        return matrix
    def vectorY(self):
        y =np.zeros(self.m)
        y[self.m-1]= self.upperSpeed
        y[0]=self.lowerSpeed
        y[1:self.m-1]=self.pressureGradient*(1/self.u)
        return y
    def function(self,y,u,upperLimit):
        x1= np.linspace(0,upperLimit, num=len(y))
        a=y[self.m-1].copy()
        b=y[0].copy()
        vx =self.pressureGradient*(1/(2*u))*(x1**2)+(((a-b)/upperLimit)-self.pressureGradient*(1/(2*u))*upperLimit)*x1+b
        return vx,x1
    def graphFunctionPlate(self,y_pos,x_direct,title,ejex,ejey):
        y_direct = np.zeros(len(y_pos))
        fig1, ax= plt.subplots()
        x_pos = np.zeros(len(x_direct))
        #con esto dibujo las flechas
        ax.quiver(x_pos,y_pos,x_direct,y_direct,angles='xy', scale_units='xy', scale=1, color='red')
        plt.title(title, 
        fontdict={'family': 'serif', 
                            'color' : 'black',
                            'weight': 'bold',
                            'size': 14})
        plt.xlabel(ejex)
        plt.ylabel(ejey)
        ax.plot(x_direct , y_pos,color='green',linewidth=2.5)
        return plt.show()    
    def compareGraph(self,funcionNumerica,k,vx,ejex,ejey,loc1):
       plt.figure()
       plt.xlabel('Eje Vx')
       plt.ylabel('Eje y')
       plt.plot(vx,k,color="blue", linewidth=2.5, linestyle="-", label="Función Real")
       plt.plot(funcionNumerica,self.x,'o',color="red",  linewidth=2.5, label="Función interpolada")
       plt.title("Grafica función")
       plt.legend(loc=loc1)
       plt.show()
       return plt.show()
    def showDatas(self,functionNumerica,k,vx):
        if self.typeFunction ==1:
            self.graphFunctionPlate(self.x,functionNumerica,"Aproximación perfil de velocidad MQ",'Vx','Y')
            self.graphFunctionPlate(k,vx,self.m,"Perfil de velocidad Real",'Vx','Y',)
            self.compareGraph(self.x,functionNumerica,k,vx,'Vx','Y','upper left')
        elif self .typeFunction==2:
            self.graphFunctionPlate(self.x,functionNumerica,"Aproximación perfil de velocidad TPS",'Vx','Y')
            self.graphFunctionPlate(k,vx,"Perfil de velocidad Real",'Vx','Y',)
            self.compareGraph(functionNumerica,k,vx,'Vx','Y','upper left')
        