#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 14:53:22 2016
Modèle de Bernouilli
@author: 2014-1527
"""

from numpy import *

pho=7.85*10**3      #kg/m3
E=2.1*10**11          #Pa
L=459*10**-3        #m
a=58*10**-3         #m largeur
b=35*10**-3         #m hauteur
n=int(input("Combien d'éléments poutre voulez-vous?"))

lp=L/n              #Longueur d'une poutre mm
Ip=a*b**3/12        #m4
mp=pho*a*b*lp       #kg

#Matrice masse d'une poutre
Mp=array([[156,22*lp,54,-13*lp],[22*lp,4*lp**2,13*lp,-3*lp**2],
             [54,13*lp,156,-22*lp],[-13*lp,-3*lp**2,-22*lp,4*lp**2]])*mp/420

#Matrice raideur d'une poutre
Kp=array([[12,6*lp,-12,6*lp],[6*lp,4*lp**2,-6*lp,2*lp**2],
             [-12,-6*lp,12,-6*lp],[6*lp,2*lp**2,-6*lp,4*lp**2]])*E*Ip/(lp**3)
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
W,u=linalg.eig(A)   #Essai avec 2 matrices M et K
L=[]

#Recheches des modes rigides
for i,freq in enumerate(W):
    if freq<=1:
        L+=[i]
L.sort()
u=delete(u,L,1)  #Suppresion des veceurs propres associés aux modes rigides
W=delete(W,L)    #suppresion des fréquences propres "nulles"
W=[W[i].real for i in range(len(W))]    #On transforme en réel (il y a des cas
#particulier ou on a des complexes avec valeurs imaginaires
#nulle ce qui posent problème pour l'utilisation de certaines fonction après)

#Récupération des fréquences propres
w=sqrt(W)
indices=[]
wtri=sorted(w)
f=[wtri[i]/(2*pi) for i in range(len(wtri))]
w=list(w)
for iwtri,elt in enumerate(wtri):
    indices+=[w.index(elt)]       #Ne fonctionne pas si 2 fréquences identiques
utri=[u[:,ind].real for ind in indices]  #Vecteur en ligne cette fois (plus simple)
#utri[0] correspond au premier vecteur propre utri[1] au deuxieme etc..
#for i in range(len(wtri)):
#    print("Fréquence f",i+1,"=",f[i]," Hz",sep="")
#    print("Vecteur propre :\n",utri[i])
print(f[0],f[1],f[2])
