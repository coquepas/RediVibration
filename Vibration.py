#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 14:10:23 2016
@author: gwendal
"""
from scipy.integrate import odeint
from numpy import *

pho=7.85*10**3     #kg/m3
E=2.1*10**11       #Pa
L=459*10**-3       #m
a=58*10**-3        #m largeur
b=35*10**-3        #m hauteur
n=int(input("Combien d'éléments voulez-vous?"))
nu=0.3

lp=L/n              #Longueur d'une poutre m
Iq=a*b**3/12        #m4
mp=pho*a*b*lp       #kg
J=mp*(lp**2+b**2)/12
#Initialisation des matrices Masse et Raideur
M=zeros(((n+1)*2,(n+1)*2))
K=zeros((2*(n+1),2*(n+1)))

Choix="1-Modèle masses concentrées\n2-Modèle Bernoulli\n3-Modèle Timoshenko"
print(Choix)
while 1:
    entree=int(input('Quel modèle choisissez-vous?'))
    if entree not in [1,2,3]:
        print('\nMauvaise entrée!\n')
        print(Choix)
    else:
        break
if entree==1:
    Ks=array([[12,6*lp,-12,6*lp],[6*lp,4*lp**2,-6*lp,2*lp**2],
              [-12,-6*lp,12,-6*lp],[6*lp,2*lp**2,-6*lp,4*lp**2]])*E*Iq/lp**3
    Ms=array([[mp,0,0,0],[0,J,0,0],[0,0,0,0],[0,0,0,0]])
    M[-2,-2]=mp
    M[-1,-1]=J
elif entree==2:
        Ks=array([[12,6*lp,-12,6*lp],[6*lp,4*lp**2,-6*lp,2*lp**2],
              [-12,-6*lp,12,-6*lp],[6*lp,2*lp**2,-6*lp,4*lp**2]])*E*Iq/lp**3
        Ms=array([[156,22*lp,54,-13*lp],[22*lp,4*lp**2,13*lp,-3*lp**2],
                [54,13*lp,156,-22*lp],[-13*lp,-3*lp**2,-22*lp,4*lp**2]])*mp/420
else:
    #Corrections Timoshenko
    eta=(12+11*nu)/(10*(1+nu))
    phi=24*eta*Iq*(1+nu)/(a*b*L**2)
    
    #Matrice masse d'une poutre
    m1=312+588*phi+280*phi**2
    m2=(44+77*phi+35*phi**2)*lp
    m3=108+252*phi+175*phi**2
    m4=(26+63*phi+35*phi**2)*lp
    m5=(8+14*phi+7*phi**2)*lp**2
    m6=(6+14*phi+7*phi**2)*lp**2
    
    #Matrice masse d'un élément
    Ms=array([[m1,m2,m3,-m4],[m2,m5,m4,-m6],[m3,m4,m1,-m2],[-m4,-m6,-m2,m5]])*mp/840
    
    
    
    #Matrice raideur d'une poutre
    Ks=array([[12,6*lp,-12,6*lp],[6*lp,(4+phi)*lp**2,-6*lp,(2-phi)*lp**2],
                 [-12,-6*lp,12,-6*lp]
                 ,[6*lp,(2-phi)*lp**2,-6*lp,(4+phi)*lp**2]])*E*Iq/(lp**3*(1+phi))



#Définition des matrices masse et raideur
for k in range(0,2*n,2):
    for i in range(4):
        for j in range(4):
            ligne=k+i
            colonne=k+j
            K[ligne,colonne]+=Ks[i,j]
            M[ligne,colonne]+=Ms[i,j]

invM=linalg.inv(M)
A=dot(invM,K)
W,u=linalg.eig(A)

L=[]

#Recheches des modes rigides
for i,freq in enumerate(W):
    if freq<=1:
        L+=[i]
L.sort()
u=delete(u,L,1)  #Suppresion des veceurs propres associés aux modes rigides
W=delete(W,L)    #suppresion des fréquences propres "nulles"
W=[W[i].real for i in range(len(W))]  #On transforme en réel (il y a des cas
#particulier ou on a des complexes avec valeurs imaginaires
#nulle ce qui posent problème pour l'utilisation de certaines fonction après)

#Récupération des fréquences propres
f=list(1/(2*pi)*sqrt(W))
indices=[]
ftri=sorted(f)
for iwtri,elt in enumerate(ftri):
    indices+=[f.index(elt)]       #Ne fonctionne pas si 2 fréquences identiques
utri=[u[:,ind].real for ind in indices]
#Vecteur en ligne cette fois (plus simple)
#utri[0] correspond au premier vecteur propre utri[1] au deuxieme etc..
#for i in range(len(wtri)):
#    print("Fréquence f",i+1,"=",f[i]," Hz",sep="")
#    print("Vecteur propre :\n",utri[i])
if len(f)>=3:
    print(ftri[0],'Hz ',ftri[1],'Hz ',ftri[2],'Hz')
