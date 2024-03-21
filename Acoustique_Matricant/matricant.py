# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 14:16:14 2024

@author: Cyril
"""

"""
Format des shapes :
    Nombre_couches    : nombre total de couche
    Nombre_pulsations : nombre total de pulsation
    Nombre_angles     : nombre d'angle
    dimension         : dimensions sptiales
    (Nombre_couches,Nombre_pulsations,Nombre_angles,dimension,dimension)
"""

import numpy as np
import scipy as sp


class Matricant():
    
    
    def __init__(self, C, epaisseur, rho, lenteur, omega, 
                 n = np.array([1,0,0], dtype=float),
                 m = np.array([0,1,0], dtype=float)):
        
        self.C = C
        self.epaisseur = epaisseur
        self.rho = rho
        self.lenteur = lenteur
        self.omega = omega
        self.n = n
        self.m = m

        
    def nm_matrix(self,n,m):
        """
        Calcul du produit nm = a_i*c_ijkl*b_l pour toute les couches.
        Calcul de A issue de n et de B issue de m
        
        INPUT
        ----------
        n  : vecteur (3x1) ou (2x1)
        m  : vecteur (3x1) ou (2x1)
        C  : tenseur (Ncx6x6)
        Nc : nombre de couche
        
        OUTPUT
        -------
        nm : (Ncx3x3)
        """
        A = np.zeros([3,6])
        B = np.zeros([6,3])
        
        #A
        A[0,0]=n[0] ; A[0,4]=n[2] ; A[0,5]=n[1]
        A[1,1]=n[1] ; A[1,3]=n[2] ; A[1,5]=n[0]
        A[2,2]=n[2] ; A[2,3]=n[1] ; A[2,4]=n[0]
        #B Transposé
        B[0,0]=m[0] ; B[4,0]=m[2] ; B[5,0]=m[1]
        B[1,1]=m[1] ; B[3,1]=m[2] ; B[5,1]=m[0]
        B[2,2]=m[2] ; B[3,2]=m[1] ; B[4,2]=m[0]

        nm = np.matmul(A,np.matmul(self.C,B))
        return nm
    
    
    def stroh_4blocs(self):
        """
        N matrix
        
        INPUT
        ----------
        n  : vecteur (3x1) ou (2x1)
        m  : vecteur (3x1) ou (2x1)
        C  : tenseur (Ncx6x6) elasticity tensor

        OUTPUT
        -------
        N1 , N2 , N3 , N4 : les quatres blocs de N (Ncx3x3) ou (Ncx2x2) sous forme de liste
        """
        nn = self.nm_matrix(self.n,self.n)
        nm = self.nm_matrix(self.n,self.m)
        mm = self.nm_matrix(self.m,self.m)
        mn = self.nm_matrix(self.m,self.n)
        nn_inv = np.linalg.inv(nn)
        
        N1 = np.matmul(nn_inv,nm)
        N2 = -nn_inv
        # N3 = mm-np.einsum('...ij,...jk,...kl->...il',mn,nn_inv,nm)
        N3 = mm - np.matmul(mn,np.matmul(nn_inv,nm))
        N4 = np.swapaxes(N1, -1, -2)
        self.N_list = [N1 , N2 , N3 , N4]
        return [N1 , N2 , N3 , N4]

    
    def reshape(self):
        """
        Format des shapes :
            Nombre_couches    : nombre total de couche
            Nombre_pulsations : nombre total de pulsation
            Nombre_angles     : nombre d'angle
            dimension         : dimensions sptiales
            (Nombre_couches,Nombre_pulsations,Nombre_angles,dimension,dimension)
        """
        nombre_couches = self.C.shape[0]
        nombre_pulsations = self.omega.shape[0]
        nombre_angles = self.lenteur.shape[0]
        dimension = self.C.shape[-1]//2
        
        self.rho.shape         = (nombre_couches,1,1,1,1)
        self.epaisseur.shape   = (nombre_couches,1,1,1,1)
        self.omega.shape       = (1,nombre_pulsations,1,1,1)
        self.lenteur.shape     = (1,1,nombre_angles,1,1)
        
        self.N_list[0].shape = (nombre_couches,1,1,dimension,dimension)
        self.N_list[1].shape = (nombre_couches,1,1,dimension,dimension)
        self.N_list[2].shape = (nombre_couches,1,1,dimension,dimension)
        self.N_list[3].shape = (nombre_couches,1,1,dimension,dimension)
    

    def matrix_Q(self):
        """
        Calcul de Q pour chaque couche, il faut faire attention à la shape des éléments
        pour que le bradcasting se fasse correctement

        INPUT
        ----------
        omega : pulsation   (1xN_omega)
        N     : Liste des blocs de N (Ncx3x3)
        rho   : mass density(Ncx1)
        s2    : lenteur
        
        OUTPUT
        -------
        Q : matrice du pb d'EDO pour chaque couche (NcxNomegax6x6) ou (NcxNomegax4x4)
        """
        Q1 = self.lenteur*self.N_list[0]
        ones = np.ones(self.lenteur.shape)
        Q2 = -self.N_list[1]*ones
        Q3 = -self.lenteur**2*self.N_list[2] + self.rho*np.identity(Q1.shape[-1])
        Q4 = self.lenteur*np.swapaxes(self.N_list[0], -1, -2)
        Q = np.block([[Q1,Q2],[Q3,Q4]])
        return 1j*self.omega*Q
    
    def matrix_Q_ext(self):
        """
        Calcul de "Q" pour chaque couche, il n'y a pas la multiplication par j*omega
        Utile si l'on veut calculer les valeurs propres et vecteurs propres de Q

        INPUT
        ----------
        omega : pulsation   (1xN_omega)
        N     : Liste des blocs de N (Ncx3x3) ou (Ncx2x2)
        rho   : mass density(Ncx1)
        s2    : lenteur
        
        OUTPUT
        -------
        Q : matrice du pb d'EDO pour chaque couche (NcxNomegax6x6) ou (NcxNomegax4x4)
        """
        Q1 = self.lenteur*self.N_list[0]
        ones = np.ones(self.lenteur.shape)
        Q2 = -self.N_list[1]*ones
        Q3 = -self.lenteur**2*self.N_list[2] + self.rho*np.identity(Q1.shape[-1])
        Q4 = self.lenteur*np.swapaxes(self.N_list[0], -1, -2)
        Q = np.block([[Q1,Q2],[Q3,Q4]])
        return Q
    

    def matricant(self,Q):
        """
        Calcul de M pour chaque couche

        INPUT
        ----------
        Q : matrice de chaque couche (NcxNomegax6x6) ou (NcxNomegax4x4)

        OUTPUT
        -------
        M : Matricant (NcxNomegax6x6) ou (NcxNomegax4x4)
        """
        M = sp.linalg.expm(self.epaisseur*Q)
        return M
    
    
    def matricant_product(self,M,dim_couche=0):
        """
        Multilayer transfer matrix

        INPUT
        ----------
        M          : Multidimensionnal array containing the transfer matrix of N layers
        dim_couche : axis of the multidimensionnal associated with the N layers (ici au vu du sizing c'est 0)

        OUTPUT
        -------
        M_tot      : Transfer matrix of the multilayer (1xNomegax6x6) ou (1xNomegax4x4)
        """
        M_tot = np.take(M, [0], axis=dim_couche)
        for i in range(1,M.shape[dim_couche]):
            matrix_layer = np.take(M, [i], axis=dim_couche)   
            # print(matrix_layer)
            # M_tot1 = np.einsum('...ij,...jk->...ik',matrix_layer,M_tot1) #pourqupoi pas utiliser matmul
            M_tot = np.matmul(matrix_layer,M_tot)
        return M_tot
    