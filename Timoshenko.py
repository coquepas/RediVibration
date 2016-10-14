#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 09:22:46 2016
Modèle Timoshenko
@author: 2014-1527
"""

from numpy import *

pho=7.85*10**3     #kg/m3
E=2.1*10**11       #Pa
L=459*10**-3       #m
a=58*1**-3         #m largeur
b=35*10**-3        #m hauteur
n=int(input("Combien d'éléments poutre voulez-vous?"))
nu=0.3

lp=L/n              #Longueur d'une poutre m
Ip=a*b**3/12        #m4
mp=pho*a*b*lp       #kg

#Corrections Timoshenko
eta=(12+11*nu)/(10*(1+nu))
phi=24*eta*Ip*(1+nu)/(a*b*L**2)

#Matrice masse d'une poutre
m1=312+588*phi+280*phi**2
m2=(44+77*phi+35*phi**2)*lp
m3=108+252*phi+175*phi**2
m4=(26+63*phi+35*phi**2)*lp
m5=(8+14*phi+7*phi**2)*lp**2
m6=(6+14*phi+7*phi**2)*lp**2

#Matrice masse d'un élément
Mp=array([[m1,m2,m3,-m4],[m2,m5,m4,-m6],[m3,m4,m1,-m2],[-m4,-m6,-m2,m5]])*mp/840



#Matrice raideur d'une poutre
Kp=array([[12,6*lp,-12,6*lp],[6*lp,(4+phi)*lp**2,-6*lp,(2-phi)*lp**2],
             [-12,-6*lp,12,-6*lp]
             ,[6*lp,(2-phi)*lp**2,-6*lp,(4+phi)*lp**2]])*E*Ip/(lp**3*(1+phi))

#Initialisation des matrices Masse et Raideur
M=zeros(((n+1)*2,(n+1)*2))
K=zeros((2*(n+1),2*(n+1)))

#Définition des matrices masse et raideur
for k in range(0,2*n,2):
    for i in range(4):
        for j in range(4):
            ligne=k+i
            colonne=k+j
            K[ligne,colonne]+=Kp[i,j]
            M[ligne,colonne]+=Mp[i,j]

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
W2=delete(W,L)    #suppresion des fréquences propres "nulles"
W2=[W2[i].real for i in range(len(W2))]  #On transforme en réel (il y a des cas
#particulier ou on a des complexes avec valeurs imaginaires
#nulle ce qui posent problème pour l'utilisation de certaines fonction après)

#Récupération des fréquences propres
w=sqrt(W2)
indices=[]
wtri=sorted(w)
f=[wtri[i]/(2*pi) for i in range(len(wtri))]
w=list(w)
for iwtri,elt in enumerate(wtri):
    indices+=[w.index(elt)]       #Ne fonctionne pas si 2 fréquences identiques
utri=[u[:,ind].real for ind in indices]
#Vecteur en ligne cette fois (plus simple)
#utri[0] correspond au premier vecteur propre utri[1] au deuxieme etc..
#for i in range(len(wtri)):
#    print("Fréquence f",i+1,"=",f[i]," Hz",sep="")
#    print("Vecteur propre :\n",utri[i])
print(f[0],'Hz ',f[1],'Hz ',f[2],'Hz')
