#Calcul du coefficient de transmission pour 
import numpy as np

class RT_fluide():
    """
    Classe permettant de calculer les coefficients de réflexion/transmission
    pour un multicouche immergé dans un fluide
    """
    
    def __init__(self,M_tot,omega,lenteur,e,rho_fluide,celerite_fluide,
                 n = np.array([1,0,0], dtype=float)):
        self.M_tot = M_tot
        self.omega = omega
        self.list_Y = self.matrice_admitance()
        self.omega = omega
        self.lenteur = lenteur
        self.e = e
        self.rho_fluide = rho_fluide
        self.celerite_fluide = celerite_fluide
        self.n = n
        self.T = self.coefficient_transmission(self.phi(np.sum(np.squeeze(self.e))),
                                               self.y_f())
        
    def matrice_admitance(self):
        """
        Calcul de la matrice Y
    
        INPUT
        ----------
        M   : Matricant total
        
        OUTPUT
        -------
        Y : Matrice d'admittance
        """
        dim = self.M_tot.shape[-1]
        M3_inv = np.linalg.inv(self.M_tot[...,dim//2:,:dim//2])
        
        Y1 = np.matmul(M3_inv,self.M_tot[...,dim//2:,dim//2:])*(-1j)
        Y2 = M3_inv*(-1j)
        Y3 = ( (np.matmul(self.M_tot[...,:dim//2,:dim//2],
                        np.matmul(M3_inv,self.M_tot[...,dim//2:,dim//2:])) 
                        - self.M_tot[...,:dim//2,dim//2:])*(-1j) )
        
        Y4 = np.matmul(self.M_tot[...,:dim//2,:dim//2],M3_inv)*(-1j)
        return np.array([Y1 , Y2 , Y3 , Y4])
    
    
    def phi(self,e_tot):
        """
        Calcul de phi
    
        INPUT
        ----------
        e_tot : Épaisseur totale du multicouche
        
        OUTPUT
        -------
        phi : Étape du calcul des coef de transmission/réflexions
        """
        return e_tot*self.omega*np.sqrt(1-(self.lenteur*self.celerite_fluide)**2)/self.celerite_fluide
    
    
    def y_f(self):
        """
        Calcul de y_f
        
        OUTPUT
        -------
        y_f : Étape du calcul des coef de transmission/réflexions
        """
        return -1j*np.sqrt(1-(self.lenteur*self.celerite_fluide)**2)/(self.rho_fluide*self.celerite_fluide)
    
    
    def coefficient_transmission(self,phi,y_f):
        '''
        Calcul du coefficient de transmission pour
        
        INPUT
        ----------
        phi    : voir document de Sasha
        y_f    : voir document de Sasha
    
        Returns
        -------
        T      : coefficient de transmission
        '''
        Y1_n = np.einsum('i,...i',self.n,np.matmul(self.list_Y[0],self.n))
        Y2_n = np.einsum('i,...i',self.n,np.matmul(self.list_Y[1],self.n))
        Y4_n = np.einsum('i,...i',self.n,np.matmul(self.list_Y[3],self.n))
        
        N_omega = Y1_n.shape[1]
        N_angles = Y1_n.shape[2]
        N_angles = 1
        Y1_n.shape = (1,N_omega,N_angles,1,1)
        Y2_n.shape = (1,N_omega,N_angles,1,1)
        Y4_n.shape = (1,N_omega,N_angles,1,1)
        
        T = -(2*np.conj(Y2_n)*y_f)*np.exp(-1j*phi)/((Y1_n-y_f)*(Y4_n-y_f)-abs(Y2_n)**2)
        return T