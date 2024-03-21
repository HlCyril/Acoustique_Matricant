# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 19:15:18 2024

@author: Cyril
"""

import numpy as np

class RT_solide():
    """
    Classe permettant de calculer les coefficients de réflexion/transmission
    pour un multicouche au sein de des solides semi-infinis
    """
    
    def __init__(self,Q_ext,M_tot):
        self.Q_ext = Q_ext
        self.M_tot = M_tot
        lambda_ext1 , xi_ext1 = self.eig_vec_solide(self.Q_ext[0,0,0,:,:])
        lambda_ext2 , xi_ext2 = self.eig_vec_solide(self.Q_ext[1,0,0,:,:])
        xi_ext1.shape = (1,1,1,6,6)
        xi_ext2.shape = (1,1,1,6,6)
        self.B = self.B_matrix(xi_ext1,xi_ext2)
        self.T = self.transmission_matrix()
        self.R = self.reflexion_matrix()
        
        
    def eig_vec_solide(self,Q_ext):
        """
        Parameters
        ----------
        Q_ext : Matrice Q
        Utilisé pour le milieu extérieur (calcul de B)

        Returns
        -------
        lambda_ord : Valeur propres de Q
        v_ord : vecteur propores de Q
        """
        Lambda , V = np.linalg.eig(Q_ext)
        
        indices_imag_p = np.where(Lambda.imag>0)
        indices_r_p = np.where((Lambda.imag==0) & (Lambda.real>0))
        indices_imag_m = np.where(Lambda.imag<0)
        indices_r_m = np.where((Lambda.imag==0) & (Lambda.real<0))
        
        i1 = np.shape(indices_imag_p)[1]
        i2 = i1 + np.shape(indices_r_p)[1]
        i3 = i2 + np.shape(indices_imag_m)[1]
        i4 = i3 + np.shape(indices_r_m)[1]
        
        lambda_ord = np.zeros(Lambda.shape,dtype=complex)
        lambda_ord[:i1]=Lambda[indices_imag_p]
        lambda_ord[i1:i2]=Lambda[indices_r_p]
        lambda_ord[i2:i3]=Lambda[indices_imag_m]
        lambda_ord[i3:i4]=Lambda[indices_r_m]
        
        v_ord=np.empty(V.shape,dtype=complex)
        v_ord[:,:i1]=V[:,indices_imag_p[0]]
        v_ord[:,i1:i2]=V[:,indices_r_p[0]]
        v_ord[:,i2:i3]=V[:,indices_imag_m[0]]
        v_ord[:,i3:i4]=V[:,indices_r_m[0]]
        return lambda_ord , v_ord


    def B_matrix(self,xi1,xi2):
        """
        Calcul de la matrice B, utile pour les coefficients de
        reflexion et de transmission
    
        INPUT
        ----------
        xi1 : matrice des vecteurs propres de Q_ext1 (milieu extérieur 1)
        M   : Matricant total
        xi2 : matrice des vecteurs propres de Q_ext2 (milieu extérieur 1)
        
        OUTPUT
        -------
        B
        """
        inv_xi1=np.linalg.inv(xi1)
        inv_M=np.linalg.inv(self.M_tot)
        B=np.einsum('...ij,...jk,...kl->...il',inv_xi1,inv_M,xi2)
        return B
    
    
    def transmission_matrix(self):
        """
        Calcul de la matrice de transmission
    
        INPUT
        ----------
        B : La matrice B
        
        OUTPUT
        -------
        T : La matrice de transmission
        """
        dim = self.B.shape[-1]
        inv_B1 = np.linalg.inv(self.B[...,:dim//2:,:dim//2])
        return inv_B1 
    
    
    def reflexion_matrix(self):
        """
        Calcul de la matrice de réflexion
    
        INPUT
        ----------
        B : La matrice B
        
        OUTPUT
        -------
        R : La matrice de réflexion
        """
        dim = self.B.shape[-1]
        inv_B1 = np.linalg.inv(self.B[...,:dim//2,:dim//2])
        B3 = self.B[...,dim//2:,:dim//2]
        R = np.einsum('...ij,...jk->...ik', B3,inv_B1)
        # R = np.matmul(B3,inv_B1)
        return R
