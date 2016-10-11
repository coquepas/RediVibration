#!/usr/bin/python3.4
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 08:49:08 2016

@author: 2014-1527
Modele masses concentrées
"""

import numpy as np

pho=7.85*10**-6     #kg/mm3
E=2*10**5           #MPa
L=459               #mm
a=58                #mm largeur
b=35                #mm hauteur
n=int(input('Combien de masses voulez-vous?'))
G=[i*L/(n-1) for i in range(n)]     #Centres de gravité
m=pho*a*b*L/n
J=m/12*((L/n)**2+b**2)
M=np.zeros((2*n,2*n))               #Initialisaiton de la matrice masse
#Définition de la matrice masse, diagonale dans ce cas
for i in range(2*n):
    if i%2==0:
        M[i][i]=m
    else:
        M[i][i]=J
K=np.zeros((2*n,2*n))               #Initialisation  de la matrice raideur
lp=L/(n-1)                          #mm
Iqua=a*b**3/12                      #mm4
Ks=np.array([[12,6*lp,-12,6*lp],[6*lp,4*lp**2,-6*lp,2*lp**2],
             [-12,-6*lp,12,-6*lp]
             ,[6*lp,2*lp**2,-6*lp,4*lp**2]])*E*Iqua/lp**3*10**-3
nel=n-1
for k in range(0,2*nel,2):
    for i in range(0,4):
        for j in range(0,4):
            ligne=k+i
            colonne=k+j
            K[ligne,colonne]=K[ligne,colonne]+Ks[i,j]

invM=np.linalg.inv(M)
A=np.dot(invM,K)
W,u=np.linalg.eig(A)
L=[]

#Recheches des modes rigides
for i,freq in enumerate(W):
    if freq<=1:
        L+=[i]
L.sort()
u=np.delete(u,L,1)  #Suppresion des veceurs propres associés aux modes rigides
W=np.delete(W,L)    #suppresion des fréquences propres "nulles"
W=[W[i].real for i in range(len(W))]    #On transforme en réel (il y a des cas
#particulier ou on a des complexes avec valeurs imaginaires
#nulle ce qui posent problème pour l'utilisation de certaines fonctions après)

#Récupération des fréquences propres
w=np.sqrt(W)
indices=[]
wtri=sorted(w)
w=list(w)
for iwtri,elt in enumerate(wtri):
    indices+=[w.index(elt)]       #Ne fonctionne pas si 2 fréquences identiques
utri=[u[:,ind].real for ind in indices]#Vecteur en ligne cette fois (plus simple)
#utri[0] correspond au premier vecteur propre utri[1] au deuxieme etc..
for i in range(len(wtri)):
    print("Fréquence f",i+1,"=",wtri[i]," Hz",sep="")
    print("Vecteur propre :\n",utri[i])
