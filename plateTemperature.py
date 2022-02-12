# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 13:53:17 2021

@author: josem


"""
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from numpy.random import seed
class FunctionsG2d:
    def functionRadial(self,r):
        if r==0.0:
            function=0
        else:
            function=(r*self.b)**(2.0*self.q)*np.log(self.b*r)
        return function
    
    def quarterCoordinate(self,x,xc,r):
        coordinate4= ((2)/((x-xc)**3))-((2)/((x-xc)**3))
        return coordinate4
    
    def thirdCoordinate(self,x,xc,r):
        coordinate3=(1/((r**2)))-(1/((r**2)))
        return coordinate3
    #derivadda de r respecto  x
    
    def secondCoordinate(self,x,xc,r):
        
        if r==0.0:
            coordinate2=0.0;
        else:
            coordinate2= (1/r)-(((x-xc)**2.0)/r**3.0)
        
        
        return coordinate2
    # puntos menos centro
    def firstCoordinate(self,x,xc,r):
        if r==0.0:
            coordinate1=0.0
        else:
            coordinate1= ((x-xc)/r)**2
        return coordinate1
    def quarterDerivateRadial(self,x,xc,r):
        if r==0:
            quarter=0
        else:
            quarter=(self.r**(2.0*self.q-4))*((32*self.q**3)-(72*self.q**2)+4*(4*self.q**3-12*self.q**2+11*self.q-3)*self.q*np.log(self.r)+44*self.q-6)
        return quarter
    def thirdDerivativeRadial(self,x,xc,r):
        if r ==0:
            third = 0
        else :
             third=((2*self.r**(2*self.q-3))*(4*(self.q**3)*np.log(self.r)+6*(self.q**2)*(1-np.log(self.r))+2*self.q*np.log(self.r)-6*self.q+1))*((x-xc)/self.r)**3
        return third

    def secondDerivativeRadial(self,r):
        if r==0:
            second=0
        else:
            second= (self.b**(2.0*self.q))*r**(2.0*(self.q-1))*(4*self.q**2.0*np.log(self.b*r)+4.0*self.q*self.b-2*self.q*np.log(self.b*r)-self.b)
        return second
    def derivativeRadial(self,r):
        if r==0:
            derivative=0
        else:
        
            derivative=(self.b)**(2.0*self.q)*r**(2.0*self.q-1.0)*(2.0*self.q*np.log(self.b*r)+self.b)
        return derivative
    
class Plate2d(FunctionsG2d):
    def __init__(self,m,n,lx,ly,temp,b,q):
        FunctionsG2d.__init__(self)
        self.m=m
        self.n=n
        self.lx=lx
        self.ly=ly
        self.temp=temp
        self.b=b
        self.q=q
    def generatorPointsPlate(self):
        x=xc=np.linspace(0,self.lx,self.m)
        y=yc=np.linspace(0,self.ly,self.n)
        #nodos=[]
        nodos=np.empty((0,4))
        formato={}
        cont=0
        for i in range (0,self.m):
            for j in range (0,self.n):
                cont = cont +1
                if x[i] ==0 and y[j] != 0 and y[j] != self.ly:
                    #nodos.append([x[i],y[j],1,self.temp['TemperaturaIzquierda']])
                    nodos=np.append(nodos,np.array([[x[i],y[j],1,self.temp['TemperaturaIzquierda']]]),axis=0)
                    formato['TemperaturaIzquierda']=int(1)
                elif x[i] ==self.lx and y[j] != 0 and y[j] != self.ly:
                    #nodos.append([x[i],y[j],2,self.temp['TemperaturaDerecha']])
                    nodos=np.append(nodos,np.array([[x[i],y[j],2,self.temp['TemperaturaDerecha']]]),axis=0)
                    formato['TemperaturaDerecha']=int(2)
                elif y[j]==0:
                    #nodos.append([x[i],y[j],3,self.temp['TemperaturaInferior']])
                    nodos=np.append(nodos,np.array([[x[i],y[j],3,self.temp['TemperaturaInferior']]]),axis=0)
                    
                    formato['TemperaturaInferior']=int(3)
                elif y[j]==self.ly:
                    #nodos.append([x[i],y[j],4,self.temp['TemperaturaSuperior']])
                    nodos=np.append(nodos,np.array([[x[i],y[j],4,self.temp['TemperaturaSuperior']]]),axis=0)
                    
                    formato['TemperaturaSuperior']=int(4)
                else:
                    #nodos.append([x[i],y[j],5,0])
                    nodos=np.append(nodos,np.array([[x[i],y[j],5,0]]),axis=0)
                    
                arrayNodos= nodos
        return x,xc,y,yc,arrayNodos,formato
   
    def matrixA(self,nodos,e):
        matriz=np.zeros((self.m*self.m,self.n*self.n))
        temp=np.zeros(self.m*self.m)
        for i in range(0,self.n*self.m):
    
            if nodos[i][2]==e['TemperaturaIzquierda']:
    
                for j in range(0,self.n*self.m):
                    r=((nodos[i][0]-nodos[j][0])**2 + (nodos[i][1]-nodos[j][1])**2)**0.5
                    if r ==0:
                        matriz[i,j]=0
                    else:
                        matriz[i,j]= self.functionRadial(r)
                    temp[i]=nodos[i][3]
                
          
            elif nodos[i][2]==e['TemperaturaDerecha']:
    
                for j in range(0,self.n*self.m):
                    r=((nodos[i][0]-nodos[j][0])**2 + (nodos[i][1]-nodos[j][1])**2)**0.5
                    if r ==0:
                        matriz[i,j]=0
                    else:
                        matriz[i,j]= self.functionRadial(r)
                    temp[i]=nodos[i][3]
            elif nodos[i][2]==e['TemperaturaInferior']:          
    
                for j in range(0,self.n*self.m):
                    r=((nodos[i][0]-nodos[j][0])**2 + (nodos[i][1]-nodos[j][1])**2)**0.5
                    if r ==0:
                        matriz[i,j]=0
                    else:
                        matriz[i,j]= self.functionRadial(r)
                    temp[i]=nodos[i][3]
                   
            elif nodos[i][2]==e['TemperaturaSuperior']:
    
                for j in range(0,self.n*self.m):
                    r=((nodos[i][0]-nodos[j][0])**2 + (nodos[i][1]-nodos[j][1])**2)**0.5
                    if r ==0:
                        matriz[i,j]=0
                    else:
                        matriz[i,j]= self.functionRadial(r)
                    temp[i]=nodos[i][3]
                 
            else :
                for j in range(0,self.n*self.m):
                    r=((nodos[i][0]-nodos[j][0])**2 + (nodos[i][1]-nodos[j][1])**2)**0.5
                    if r==0:
                        matriz[i,j]=0
                    else:
                        matriz[i,j]= self.secondDerivativeRadial(r)* self.firstCoordinate(nodos[i][0],nodos[j][0],r)+self.derivativeRadial(r)*self.secondCoordinate(nodos[i][0],nodos[j][0],r) + self.secondDerivativeRadial(r)* self.firstCoordinate(nodos[i][1],nodos[j][1],r)+self.derivativeRadial(r)*self.secondCoordinate(nodos[i][1],nodos[j][1],r)
                    temp[i]=nodos[i][3]
        return matriz, temp
    def pesosW(self,matrix,y):
        matrixInverse = np.linalg.inv(matrix)
        pesosW =  np.matmul(matrixInverse, y)
        return pesosW    
    def functionGTPS(self,nodos):
        r =np.zeros((self.m*self.m,self.n*self.n))
        for i in range(0,self.n*self.m):
            for j in range(0,self.n*self.m):
                 r[i,j]=((nodos[i][0]-nodos[j][0])**2 + (nodos[i][1]-nodos[j][1])**2)**0.5
    
        gTPS = (r*self.b)**(2.0*self.q)*np.log(self.b*r)
        gTPS =np.where(r==0,0,gTPS)
        derivadagTPS =(self.b)**(2.0*self.q)*r**(2.0*self.q-1.0)*(2.0*self.q*np.log(self.b*r)+self.b)
        derivadagTPS =np.where(r==0,0,derivadagTPS)
        segundaDerivadagTPS = (self.b**(2.0*self.q))*r**(2.0*(self.q-1))*(4*self.q**2.0*np.log(self.b*r)+4.0*self.q*self.b-2*self.q*np.log(self.b*r)-self.b)
        segundaDerivadagTPS=np.where(r==0,0,segundaDerivadagTPS)
        return gTPS,derivadagTPS,segundaDerivadagTPS,r
    def numericalApproximations(self,gMQ,derivadagMQ,segundaDerivadagMQ,pesosW):
        funcionNumerica=np.matmul(gMQ,pesosW)
        derivadaNumerica=np.matmul(derivadagMQ,pesosW)
        segundaDerivadaNumerica=np.matmul(segundaDerivadagMQ,pesosW)
        return funcionNumerica,derivadaNumerica,segundaDerivadaNumerica
    def postProceso(self,x1,y1,z1,text):
        seed(1234)
        x = x1
        y =y1
        z =z1
        # define grid.
        xi = np.linspace(0,self.lx,1000)
        yi = np.linspace(0,self.ly,1000)
        # grid the data.
        zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
        plt.contourf(xi,yi,zi,50,vmin=0,vmax=100,cmap=plt.cm.jet)
        m = plt.cm.ScalarMappable(cmap=plt.cm.jet)
        m.set_array(zi)
        m.set_clim(0., 100)
        plt.colorbar(m, boundaries=np.linspace(0, 100, 10))
        plt.xlim(0,self.lx)
        plt.ylim(0,self.ly)
        plt.title(text)
        return plt.show()
    def analyticalTemperature(self,e,t):
        a=self.lx
        b=self.ly
        temperature =[]
        sumatoria=0
        
        for l in range(0,self.m*self.m):
            for i in range(1,114):
                contador=2*i-1
                sumatoria = sumatoria+(np.sin((np.pi*contador*e[l][0])/a)*np.sinh((np.pi*contador*(b-e[l][1]))/a))/(contador*(np.sinh(np.pi*b*contador/a)))
             
            #     print(i)
            # print(sumatoria, "sumatoria")
            # print("")
            temperature.append(sumatoria*(4*t)/np.pi)
            sumatoria=0
        arrayTemperatura=np.asarray(temperature)        
        return arrayTemperatura
            
