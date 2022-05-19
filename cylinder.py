# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 16:37:16 2021

@author: josem
"""
import numpy as np
from interpolation import InterpolationFBR
from flowbetweenplates import ProcessPlates
import matplotlib.pyplot as plt
class FlowCylinder(InterpolationFBR):
      def __init__(self,m,n,lowerLimit,upperLimit):
        InterpolationFBR.__init__(self,m,n,lowerLimit,upperLimit)
class FunctiosG:
    def functionTps(self,x,xc,q,b):
    	 r=abs(x-xc)
    	 if r==0.0:
    	  	frbf=0.0
    	 else: 
    	  	frbf=(r*b)**(2.0*q)*np.log(b*r)
    	 return frbf
    def derivateTPS(self,x,xc,q,b):
    	r=abs(x-xc)
    	if r==0:
    		f_dev=0.0
    	else:
    		f_dev=(b)**(2.0*q)*r**(2.0*q-1.0)*(2.0*q*np.log(b*r)+b)
    	return f_dev
    def secondDerivateTPS(self,x,xc,q,b):
    	r=abs(x-xc)
    	if r==0:
    		f_seg_dev=0.0
    	else:
    		f_seg_dev=(b**(2.0*q))*r**(2.0*(q-1))*(4*q**2.0*np.log(b*r)+4*q*b-2*q*np.log(b*r)-b)
    	return f_seg_dev
    def parameter1(self,x,xc):
    	r=abs(x-xc)
    	if r==0.0:
    		par=0.0
    	else:
    		par= (x-xc)/r
    	return par
class ProcessCylinder(ProcessPlates,FunctiosG):
    def __init__(self, r, m, q, parameterMQ, typeFunction, x, upperSpeed, lowerSpeed, pressureGradient, u,radioCylinder,n,xc,b=1):
        ProcessPlates.__init__(self, r, m, q,parameterMQ ,typeFunction, x, upperSpeed, lowerSpeed, pressureGradient, u)
        FunctiosG.__init__(self)
        self.radioCylinder =radioCylinder
        self.n=n
        self.xc=xc
        self.b=b
    def matrixA(self,pointsMonosCenters,r):
        A = np.zeros((self.m,self.n))
        if self.typeFunction==2:
            for i in range(0,self.m):
                for j in range(0,self.n):
                    if i== 0:
                        A[i,j]=self.functionTps(self.x[i],self.xc[j],self.q,self.b)
                    elif i == (self.m-1):
                        A[i,j]=self.functionTps(self.x[i],self.xc[j],self.q,self.b)
                    else:
                        A[i,j]=(self.x[i])*self.secondDerivateTPS(self.x[i],self.xc[j],self.q,self.b)+self.derivateTPS(self.x[i],self.xc[j],self.q,self.b)* self.parameter1(self.x[i],self.xc[j])
        elif self.typeFunction==1:
            A=self.x*(((self.r**2)+(self.parameterMQ)-pointsMonosCenters**2)/((r**2+self.parameterMQ)**1.5))+(((r**2)+self.parameterMQ)**(-0.5))*(pointsMonosCenters)
            A[0,:]= np.sqrt((self.r[0,:]**2+self.parameterMQ))
            A[self.m-1,:]= np.sqrt((self.r[self.m-1,:]**2+self.parameterMQ)) 
        else:
            print("Error")
        return A
    def vectorY(self):
        y =np.zeros(self.m)
        y[self.m-1]= self.upperSpeed
        y[0]=self.lowerSpeed
        y[1:self.m-1]=self.pressureGradient*(1/self.u)
        return y*self.x
    def function(self):
        x1= np.linspace(-self.radioCylinder,self.radioCylinder, num=len(self.x))
        vx = ((x1**2-self.radioCylinder**2)/(4*self.u))*(self.pressureGradient)
        return vx,x1
    def graphFunction(self,y_pos,x_direct,title,ejex,ejey,space=40,large=350):
        y_direct = np.zeros(len(y_pos))
        fig1, ax= plt.subplots()
        x_pos = np.zeros(len(x_direct))
        #con esot dibujo las flechas
        for i in range(0,large,space):
            ax.quiver(x_pos+i,y_pos,x_direct,y_direct,angles='xy', scale_units='xy', scale=1, color='red')
            plt.title(title, 
                  fontdict={'family': 'serif', 
                            'color' : 'black',
                            'weight': 'bold',
                            'size': 14})
            plt.xlabel(ejex)
            plt.ylabel(ejey)
            ax.plot(x_direct + i, y_pos,color='green',linewidth=2.5)
        # linea  vertical     
        ax.axhline(-y_pos[0], -10, large+20,color = 'black')
        ax.axhline(y_pos[0], -10, large+20,color = 'black')
        # ax.quiver(x_pos,y_pos,x_direct,y_direct,angles='xy', scale_units='xy', scale=1, color='red')
        # plt.title(title, 
        # fontdict={'family': 'serif', 
        #                     'color' : 'black',
        #                     'weight': 'bold',
        #                     'size': 14})
        # plt.xlabel(ejex)
        # plt.ylabel(ejey)
        #ax.plot(x_direct , y_pos,color='green',linewidth=2.5)
        return plt.show()
    def showCylinder(self,functionNumerica,k,vx):
        if self.typeFunction ==1:
            self.graphFunction(self.x,functionNumerica,"Aproximación perfil de velocidad MQ",'Vz','r',40,350)
            self.graphFunction(k,vx,"Perfil de velocidad real",'Vz','r',40,350)
        elif self .typeFunction==2:
            self.graphFunction(self.x,functionNumerica,"Aproximación perfil de velocidad TPS",'Vz','r',40,350)
            self.graphFunction(k,vx,"Perfil de velocidad real",'Vz','r',40,350)
            
        
