# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 18:56:48 2022

@author: josem
"""
import numpy as np
import scipy.optimize as opti
class NoLineal() :
    def __init__(self,matrix,temperature,tipo):
        self.matrix= matrix
        self.temperature=temperature
        self.matrixArmada = self.montarMatrix()
        self.tipo=tipo
    def montarMatrix(self):
      return np.insert(self.matrix,self.matrix.shape[1],self.temperature*(-1),1)
    def CalcularSolucionNoLineal(self,x):
        if self.tipo==1 :
            ll=[]
            f=open('kk.txt','w')
            for j in range (0,len(self.matrixArmada[:,0])):
                    a=""
                    for i in range(0,len(self.matrixArmada[0,:])):
                        if ( i < len(self.matrixArmada[0,:])-1 ):
                            a = a + str(self.matrixArmada[j,i])+"*"+"x"+str([i])+"+"
                        else:
                            a=a+str(self.matrixArmada[j,i])
                    f.write(a+'\n')
                    ll.append(eval(a))      
            return ll
        else:
            lineas=[]
            with open('kk.txt', 'r') as f:
                for linea in f:
                    lineas.append(eval(linea) )
            return lineas