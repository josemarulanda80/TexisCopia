# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 19:47:52 2022

@author: josem
"""

from plateTemperature import FunctionsG2d
import numpy as np
import sys

import matplotlib.pyplot as plt
class MatrixNavier(FunctionsG2d):
    def __init__(self,m,n,lx,ly,b,q,deltaPx,deltaPy,velocidad,typeF):
        self.m=m
        self.n=n
        self.lx=lx
        self.ly=ly
        self.b=b
        self.q=q
        self.deltaPx=deltaPx
        self.deltaPy=deltaPy
        self.velocidad=velocidad
        self.typeF=typeF  
    def generarPuntos(self):
        x=xc=np.linspace(0,self.lx/self.lx,self.m)
        y=yc=np.linspace(0,self.ly,self.n)
        inf_nodos=[]
        for i in range(0,self.m):
            for j in range(0,self.n):
                if x[i]==0 and y[j]!=0 and y[j]!=self.ly:
                #if x[i]==0:
                    inf_nodos.append([x[i],y[j],1,self.deltaPx]) # 1: frontera izquierda                
                elif x[i]==self.lx/self.lx and y[j]!=0 and y[j]!=self.ly:
                #elif x[i]==l_x:
                    inf_nodos.append([x[i],y[j],2,self.deltaPy])  # 2: frontera derecha               
                elif y[j]==0:
                #elif y[j]==0 and x[i]!=0 and x[i]!=l_x:
                    inf_nodos.append([x[i],y[j],3,self.velocidad['vyInferior']]) # 3: frontera inferior                
                elif y[j]==self.ly:
                #elif y[j]==l_y and x[i]!=0 and x[i]!=l_x:
                    inf_nodos.append([x[i],y[j],4,self.velocidad['vxSuperior']]) # 4: frontera superior                
                else:
                    inf_nodos.append([x[i],y[j],5,0]) # 5: interiorr 
                nodos_mx=np.asarray(inf_nodos)
        return x,xc,y,yc,nodos_mx
    
    def generarPuntosBlock(self):
        x=xc=np.linspace(0,self.lx/self.lx,self.m)
        y=yc=np.linspace(0,self.ly,self.n)
        aux=abs(x[0]-x[1])
        suma=0
        suma1=0
        while suma <= 3.00/self.lx:
             suma=suma+aux
        suma = suma -aux
        while suma1 <=6/self.lx:
            suma1=suma1+aux
        print(suma)
        yy= np.linspace(0.5,self.ly,self.n)
        
        inf_nodos=[]
        for i in range(0,self.m):
            for j in range(0,self.n):
                if x[i]==0 and y[j]!=0 and y[j]!=self.ly:
                    inf_nodos.append([x[i],y[j],1,0]) # 1: frontera izquierda                
                elif x[i]==self.lx/self.lx and y[j]!=0 and y[j]!=self.ly:
                    inf_nodos.append([x[i],y[j],2,0])  # 2: frontera derecha  
                elif y[j]==self.ly:
                    inf_nodos.append([x[i],y[j],4,10]) # 4: frontera superior 
                elif y[j]==0 and x[i]<(3/self.lx)  :
                    inf_nodos.append([x[i],y[j],3,0]) # 3: frontera inferior   
                elif y[j]==0 and x[i]>(6/self.lx)  :
                    inf_nodos.append([x[i],y[j],3,0]) # 3: frontera inferior   
                elif x[i]>=(3/self.lx) and x[i]<=(6/self.lx) :
                    if y[j]==0.5:
                      inf_nodos.append([x[i],y[j],6,0]) # Frontera inferior
                    elif y[j]>0.5:
                      inf_nodos.append([x[i],y[j],8,0]) # Frontera inferior
                elif round(x[i],2)== round(suma,2) or round(x[i],2)== round(suma1,2):
                    if y[j]<=0.5: 
                        inf_nodos.append([x[i],y[j],7,0]) 
                    else:
                            inf_nodos.append([x[i],y[j],8,0]) 
                else:
                      inf_nodos.append([x[i],y[j],5,0]) # 5: interiorr 
                      
                     
                    
                nodos_mx=np.asarray(inf_nodos)
        plt.plot(nodos_mx[:,0],nodos_mx[:,1], "o")
        return x,xc,y,yc,nodos_mx
    def matrix(self):
        x,xc,y,yc,nodos_mx=self.generarPuntos()
        print(len(nodos_mx[:,0]))
        m = self.m*self.n
        n = self.m*self.n
        r=np.zeros((m,n))
        coord1=0 # coordenada x
        coord2=1 # coordenada y                   
        with open("Eq_mx_aumentada.txt", "w") as f:
            for i in range (0,m):
                aux_vx = ""
                aux_vy = ""
                aux_interior = ""
                ecu_continuidad = ""
                aux_vx_aumen = ""            
                aux_vy_aumen = ""           
                aux_interior_aumen = ""            
                ecu_continuidad_aumen = ""
                deltax = 3.33333-0
                if nodos_mx[i][2]==4:
                    for j in range (0,n):                   
                        r = ((nodos_mx[i][coord1] - nodos_mx[j][coord1])**2 + (nodos_mx[i][coord2] - nodos_mx[j][coord2])**2)**0.5
                        ffrbf = self.functionRadial(r) 
                        aux_vx += "x"+str([j])+"*"+"("+str(ffrbf)+")"+"+"
                        aux_vy += "x"+str([j+n])+"*"+"("+str(ffrbf)+")"+"+"                        
                        # aux_vx_aumen = aux_vx
                        # aux_vy_aumen = aux_vy
                        aux_vx_aumen = aux_vx +  str(-self.velocidad['vxSuperior'])
                        aux_vy_aumen = aux_vy +  str(-self.velocidad['vySuperior']) 
                    # f.write(aux_vx_aumen+str(-self.velocidad['vxSuperior'])+"\n")
                    # f.write(aux_vy_aumen+str(-self.velocidad['vySuperior'])+"\n")
                    f.write(aux_vx_aumen+"\n")
                    f.write(aux_vy_aumen+"\n")
                elif nodos_mx[i][2]==3:
                    for j in range (0,n):
                        r = ((nodos_mx[i][coord1] - nodos_mx[j][coord1])**2 + (nodos_mx[i][coord2] - nodos_mx[j][coord2])**2)**0.5
                        ffrbf = self.functionRadial(r) 
                        aux_vx += "x"+str([j])+"*"+"("+str(ffrbf)+")"+"+"
                        aux_vy += "x"+str([j+n])+"*"+"("+str(ffrbf)+")"+"+"  
                        # aux_vx_aumen = aux_vx 
                        # aux_vy_aumen = aux_vy 
                        aux_vx_aumen = aux_vx +  str(-self.velocidad['vxInferior'])
                        aux_vy_aumen = aux_vy +  str(-self.velocidad['vyInferior'])                        
                    # f.write(aux_vx_aumen+str(-self.velocidad['vxInferior'])+"\n")                
                    # f.write(aux_vy_aumen+str(-self.velocidad['vyInferior'])+"\n")                        
                    f.write(aux_vx_aumen+"\n")                
                    f.write(aux_vy_aumen+"\n")
                else:
                    for j in range (0,n):
                        b=self.lx/self.ly
                        mu=100
                        rho=1
                    #     r=((nodos_mx[i][0]-nodos_mx[j][0])**2 + (nodos_mx[i][1]-nodos_mx[j][1])**2)**0.5
                    #     ffrbf = self.functionRadial(r)
                    #     fdrbf_dx =self.derivativeRadial(r)*self.firstCoordinate(nodos_mx[i][0],nodos_mx[j][0], r)
                    #     fdrbf_dy =self.derivativeRadial(r)*self.firstCoordinate(nodos_mx[i][1],nodos_mx[j][1], r)
                    #     fd2rbf_dx2 = self.secondDerivativeRadial(r)*(self.firstCoordinate(nodos_mx[i][0],nodos_mx[j][0],r)**2)+self.derivativeRadial(r)*self.secondCoordinate(nodos_mx[i][0],nodos_mx[j][0],r) 
                    #     fd2rbf_dy2 = self.secondDerivativeRadial(r)*(self.firstCoordinate(nodos_mx[i][1],nodos_mx[j][1],r)**2)+self.derivativeRadial(r)*self.secondCoordinate(nodos_mx[i][1],nodos_mx[j][1],r)
                    #         # self.secondDerivativeRadial(r)* (self.firstCoordinate(nodos_mx[i][0],nodos_mx[j][0],r)**2)+self.derivativeRadial(r)*self.secondCoordinate(nodos_mx[i][0],nodos_mx[j][0],r) 
                    #     aux_interior += str(1/((self.ly*100*rho)/mu))+"*"+"(" +str((b**2)*(fd2rbf_dx2 + fd2rbf_dy2)) +")"+"*"+"(" + "x"+str([j])+"+" +"x"+str([j+n])  +")" + "-"+"(" +str(b*(ffrbf * fdrbf_dx))+")"+"*"+"(" + "x"+str([j])+")"+"**"+str(2.0)+"-"+"(" +str(b*(ffrbf * fdrbf_dy))+")"+"*"+"(" + "x"+str([j+n])+")"+"**"+str(2.0)+"-"+"(" +str(ffrbf *(fdrbf_dx + fdrbf_dy)*b)+")"+"*"+"(" + "x"+str([j])+"*"+ "x"+str([j+n])+")" +"+"

                    #     ecu_continuidad += "x"+str([j])+"*"+"("+str(fdrbf_dx)+")"+"+"+"x"+str([j+n])+"*"+"("+str(fdrbf_dy)+")" +"+"
                             
                    #     aux_interior_aumen = aux_interior + str((1/1)*((b*self.deltaPx - b*self.deltaPy)/self.lx))
                    #     #aux_interior_aumen = aux_interior + str((1/1)*((b*self.deltaPx - b*self.deltaPy)/self.lx))
                    #     ecu_continuidad_aumen = ecu_continuidad + str(0.0)  
                    # f.write(aux_interior_aumen+"\n")
                    # f.write(ecu_continuidad_aumen+"\n")
                        r = ((nodos_mx[i][coord1] - nodos_mx[j][coord1])**2 + (nodos_mx[i][coord2] - nodos_mx[j][coord2])**2)**0.5
                        ffrbf = self.functionRadial(r)
                        fdrbf_dx =self.derivativeRadial(r)*self.firstCoordinate(nodos_mx[i][0],nodos_mx[j][0], r)
                        fdrbf_dy = self.derivativeRadial(r)*self.firstCoordinate(nodos_mx[i][1],nodos_mx[j][1], r)
                        fd2rbf_dx2 = self.secondDerivativeRadial(r)*(self.firstCoordinate(nodos_mx[i][0],nodos_mx[j][0],r)**2)+self.derivativeRadial(r)*self.secondCoordinate(nodos_mx[i][0],nodos_mx[j][0],r) 
                        fd2rbf_dy2 = self.secondDerivativeRadial(r)*(self.firstCoordinate(nodos_mx[i][1],nodos_mx[j][1],r)**2)+self.derivativeRadial(r)*self.secondCoordinate(nodos_mx[i][1],nodos_mx[j][1],r)
                        # aux_interior += str(1/((self.ly*10*rho)/(mu*b)))+"*"+"(" +str((b)*(fd2rbf_dx2 + fd2rbf_dy2)*(b)) +")"+"*"+"(" + "x"+str([j])+"+" +"x"+str([j+n])  +")" +\
                        # "-"+"(" +str((b)*(ffrbf * fdrbf_dx))+")"+"*"+"(" + "x"+str([j])+")"+"**"+str(2.0)+\
                        # "-"+"(" +str((b)*(ffrbf * fdrbf_dy))+")"+"*"+"(" + "x"+str([j+n])+")"+"**"+str(2.0)+\
                        # "-"+"(" +str((b)*ffrbf *(fdrbf_dx + fdrbf_dy))+")"+"*"+"(" + "x"+str([j])+"*"+ "x"+str([j+n])+")" +"+"
                        aux_interior += str(1/100)+"*"+"(" +str((fd2rbf_dx2 + (b**2)*fd2rbf_dy2)) +")"+"*"+"(" + "x"+str([j])+"+" +"x"+str([j+n])  +")" +"+"
                        #aux_interior=str(0)
                        ecu_continuidad += "x"+str([j])+"*"+"("+str(fdrbf_dx)+")"+"+"+"x"+str([j+n])+"*"+"("+str(b*fdrbf_dy)+")" +"+"
                        aux_interior_aumen = aux_interior+str((b)*((self.deltaPx - self.deltaPy)/self.lx))
                        ecu_continuidad_aumen = ecu_continuidad + str(0.0) 
                    f.write(aux_interior_aumen+"\n")
                    f.write(ecu_continuidad_aumen+"\n")
    def matrixTwo(self):
            x,xc,y,yc,nodos_mx=self.generarPuntos()
            m = len(nodos_mx[:,0])
            n = len(nodos_mx[:,0])
            r=np.zeros((m,n))
            coord1=0 # coordenada x
            coord2=1 # coordenada y                   
            with open("Eq_mx_aumentada.txt", "w") as f:
                for i in range (0,m):
                    aux_vx = ""
                    aux_vy = ""
                    aux_interior = ""
                    ecu_continuidad = ""
                    aux_vx_aumen = ""            
                    aux_vy_aumen = ""           
                    aux_interior_aumen = ""            
                    ecu_continuidad_aumen = ""
                    deltax = 3.33333-0
                    if nodos_mx[i][2]==4:
                        for j in range (0,n):                   
                            r = ((nodos_mx[i][coord1] - nodos_mx[j][coord1])**2 + (nodos_mx[i][coord2] - nodos_mx[j][coord2])**2)**0.5
                            ffrbf = self.functionRadial(r) 
                            aux_vx += "x"+str([j])+"*"+"("+str(ffrbf)+")"+"+"
                            aux_vy += "x"+str([j+n])+"*"+"("+str(ffrbf)+")"+"+"                        
                            aux_vx_aumen = aux_vx +  str(-self.velocidad['vxSuperior'])
                            aux_vy_aumen = aux_vy +  str(-self.velocidad['vySuperior']) 
                        f.write(aux_vx_aumen+"\n")
                        f.write(aux_vy_aumen+"\n")
                    elif nodos_mx[i][2]==3 or  nodos_mx[i][2]==7 or nodos_mx[i][2]==6:
                        for j in range (0,n):
                            r = ((nodos_mx[i][coord1] - nodos_mx[j][coord1])**2 + (nodos_mx[i][coord2] - nodos_mx[j][coord2])**2)**0.5
                            ffrbf = self.functionRadial(r) 
                            aux_vx += "x"+str([j])+"*"+"("+str(ffrbf)+")"+"+"
                            aux_vy += "x"+str([j+n])+"*"+"("+str(ffrbf)+")"+"+"  
                            aux_vx_aumen = aux_vx +  str(-self.velocidad['vxInferior'])
                            aux_vy_aumen = aux_vy +  str(-self.velocidad['vyInferior'])                                               
                        f.write(aux_vx_aumen+"\n")                
                        f.write(aux_vy_aumen+"\n")
                    elif nodos_mx[i][2]==8:
                         b=self.lx/self.ly
                         q=6.36612
                         mu=1
                         for j in range (0,n):
                            r = ((nodos_mx[i][coord1] - nodos_mx[j][coord1])**2 + (nodos_mx[i][coord2] - nodos_mx[j][coord2])**2)**0.5
                            ffrbf = self.functionRadial(r)
                            fdrbf_dx =self.derivativeRadial(r)*self.firstCoordinate(nodos_mx[i][0],nodos_mx[j][0], r)
                            fdrbf_dy = self.derivativeRadial(r)*self.firstCoordinate(nodos_mx[i][1],nodos_mx[j][1], r)
                            fd2rbf_dx2 = self.secondDerivativeRadial(r)*(self.firstCoordinate(nodos_mx[i][0],nodos_mx[j][0],r)**2)+self.derivativeRadial(r)*self.secondCoordinate(nodos_mx[i][0],nodos_mx[j][0],r) 
                            fd2rbf_dy2 = self.secondDerivativeRadial(r)*(self.firstCoordinate(nodos_mx[i][1],nodos_mx[j][1],r)**2)+self.derivativeRadial(r)*self.secondCoordinate(nodos_mx[i][1],nodos_mx[j][1],r)
                            aux_interior += str(1/100)+"*"+"(" +str((fd2rbf_dx2 + (b**2)*fd2rbf_dy2)) +")"+"*"+"(" + "x"+str([j])+"+" +"x"+str([j+n])  +")" +"+"
                            ecu_continuidad += "x"+str([j])+"*"+"("+str(fdrbf_dx)+")"+"+"+"x"+str([j+n])+"*"+"("+str(b*fdrbf_dy)+")" +"+"
                            aux_interior_aumen = aux_interior+str((b)*(((8*mu)/(np.pi()))*((q)/(0.5**4))))
                            ecu_continuidad_aumen = ecu_continuidad + str(0.0) 
                         f.write(aux_interior_aumen+"\n")
                         f.write(ecu_continuidad_aumen+"\n")
                                              
                    else:
                        for j in range (0,n):
                            b=self.lx/self.ly
                            r = ((nodos_mx[i][coord1] - nodos_mx[j][coord1])**2 + (nodos_mx[i][coord2] - nodos_mx[j][coord2])**2)**0.5
                            ffrbf = self.functionRadial(r)
                            fdrbf_dx =self.derivativeRadial(r)*self.firstCoordinate(nodos_mx[i][0],nodos_mx[j][0], r)
                            fdrbf_dy = self.derivativeRadial(r)*self.firstCoordinate(nodos_mx[i][1],nodos_mx[j][1], r)
                            fd2rbf_dx2 = self.secondDerivativeRadial(r)*(self.firstCoordinate(nodos_mx[i][0],nodos_mx[j][0],r)**2)+self.derivativeRadial(r)*self.secondCoordinate(nodos_mx[i][0],nodos_mx[j][0],r) 
                            fd2rbf_dy2 = self.secondDerivativeRadial(r)*(self.firstCoordinate(nodos_mx[i][1],nodos_mx[j][1],r)**2)+self.derivativeRadial(r)*self.secondCoordinate(nodos_mx[i][1],nodos_mx[j][1],r)
                            aux_interior += str(1/100)+"*"+"(" +str((fd2rbf_dx2 + (b**2)*fd2rbf_dy2)) +")"+"*"+"(" + "x"+str([j])+"+" +"x"+str([j+n])  +")" +"+"
                            ecu_continuidad += "x"+str([j])+"*"+"("+str(fdrbf_dx)+")"+"+"+"x"+str([j+n])+"*"+"("+str(b*fdrbf_dy)+")" +"+"
                            aux_interior_aumen = aux_interior+str((b)*((self.deltaPx - self.deltaPy)/self.lx))
                            ecu_continuidad_aumen = ecu_continuidad + str(0.0) 
                        f.write(aux_interior_aumen+"\n")
                        f.write(ecu_continuidad_aumen+"\n")
    def noLinearSistem(self,x):
        if(self.typeF==0):
            self.matrix()
            sys.exit()
        else:
            lineas=[]
            with open('Eq_mx_aumentada.txt', 'r') as f:
                for linea in f:
                    lineas.append(eval(linea) )
            return lineas
    def noLinearSistemTwo(self,x):
        if(self.typeF==0):
            self.matrixTwo()
            sys.exit()
        else:
            lineas=[]
            with open('Eq_mx_aumentada.txt', 'r') as f:
                for linea in f:
                    lineas.append(eval(linea) )
            return lineas
    def Graf_speed_2D(self,nx,ny,l_x,l_y,nodos_mx,velx_num,vely_num):
        from scipy.interpolate import griddata as gd
        import matplotlib.pyplot as plt
        print("Velx_num_max =", max(velx_num))
        ############## Creating arrow #################################
        x = self.lx*nodos_mx[:,0]
        y = nodos_mx[:,1]
        xi = np.linspace(0,l_x,nx)
        yi = np.linspace(0,l_y,ny)
        X,Y = np.meshgrid(xi,yi)
        v_x = gd((x,y), velx_num, (X, Y), method='cubic')
        v_y = gd((x,y), vely_num, (X, Y), method='cubic')
        fig, ax = plt.subplots(figsize =(12, 4))
        #  Escala para las flechas
        # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.quiver.html#matplotlib.axes.Axes.quiver
        h1=ax.quiver(X, Y, v_x, v_y, color='b')
        # , angles='xy', scale_units='xy', scale=1,
        #scale=1, units='xy'
        #ax.xaxis.set_ticks([])
        #ax.yaxis.set_ticks([])
        # ax.axis([0,self.lx+2,0,1.3])
        font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 16,
        }
        # ax.quiverkey(h1,                      # Mango de carcaj entrante
        #      X=X, Y=Y,       #Determina la ubicación de la etiqueta, todo limitado a [0,1]
        #      U=v_x,                    # La longitud de la flecha de referencia significa que la velocidad del viento es de 5 m / s.
        #      angle = 0,            # Ángulo de ubicación de la flecha de referencia. El valor predeterminado es 0, lo que significa ubicación horizontal
        #      label='v:5m/s',        # Suplemento de Arrow: contenido de la etiqueta + 
        #      labelpos='S',          #label es en qué dirección de la flecha de referencia; S significa sur
        #      color = 'b',labelcolor = 'b', # Color de flecha + color de etiqueta
        #      fontproperties = font,        # configuración de fuente de la etiqueta: tamaño, estilo, peso
        #      )      

        # # Debido a que el viento tiene dos direcciones de U \ V, es mejor establecer flechas de referencia en dos direcciones + etiqueta
        # ax.quiverkey(h1, X=0.07, Y = 0.071,   
        #               U = 5, 
        #               angle = 90,           # Ángulo de ubicación de la flecha de referencia, 90 significa ubicación vertical    
        #               label = 'w:5cm/s',    #labelContent
        #               labelpos='N',         #label está al norte de la flecha de referencia
        #               color = 'r',          #Flecha color
        #               labelcolor = 'r',     #labelcolor
        #               fontproperties = font)

        
        #etiquetas ejes
        plt.xlabel(r"$L$")
        plt.ylabel("$D$")
        ax.set_title('Perfil de Velocidad')  
        plt.scatter(X, Y, color='b', s=5)
        plt.savefig('PerfilVelocidad sin delta de presión_Prueba.png', dpi = 600)
        plt.show()        






