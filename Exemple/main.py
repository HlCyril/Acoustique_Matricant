import Acoustique_Matricant
import yaml
import numpy as np


def RT_coefficient_periodic_analytic_duralumin_epoxy(omega,Ncell):
    """
    Analytical amplitude transmisson and reflexion coefficients (Sasha Shuvalov)
    Configuration :
   
    Water | Duralumin | Epoxy | Duralumin | Epoxy .... | Duralumin | Epoxy | Water

    INPUT
    ----------
   
    OUTPUT
    -------
    """

    COS = np.cos(0.0785*omega)*np.cos(0.2*omega)
    SIN = np.sin(0.0785*omega)*np.sin(0.2*omega)
    ACOS = np.arccos(COS-(1/2)*6.64*SIN)
    SIN_COS = np.sin(0.0785*omega)*np.cos(0.2*omega)
    COS_SIN = np.cos(0.0785*omega)*np.sin(0.2*omega)

    DR1 = 4*np.sin(ACOS)**2
    DR2 = (40.1*SIN**2+(11.936*SIN_COS+1.3135*COS_SIN)**2)*np.sin(Ncell*np.arccos(COS-(1/2)*6.64*SIN))**2
    R_mod_carre = 1/(1+DR1/(DR2))
    T_mod_carre = 1/(1+DR2/(DR1))
   
    return R_mod_carre,T_mod_carre


if __name__ == "__main__":
    
    test = 2
    
    if test == 1:
    
        N_cell = 6
        
        with open('deck.yaml','r') as fichier:
            data_materiaux = yaml.safe_load(fichier)
        
        constantes_physique = Acoustique_Matricant.Proprietes_Physiques(data_materiaux,'isotrope')
        dictionnaire = constantes_physique.dic
        
        rho , cij , e = constantes_physique.periodic_multilayer(N_cell)
        
        with open('deck_ext.yaml','r') as fichier:
            data_exterieur = yaml.safe_load(fichier)
        
        rho_ext = data_exterieur['1']['rho']
        celerite = data_exterieur['1']['c']
    
        freq       = np.linspace(0.,5,num=1001,dtype=float)
        omega      = freq[1:]*(2*np.pi)
        N_omega = omega.shape[0]
        
        theta = np.array([0.])
        lenteur = np.sin(theta)/celerite
        
        calculs_matricant = Acoustique_Matricant.Matricant(cij, e, rho, lenteur, omega)
       
        N = calculs_matricant.stroh_4blocs()
        
        calculs_matricant.reshape()
        Q = calculs_matricant.matrix_Q()
        M = calculs_matricant.matricant(Q)
        M_tot=calculs_matricant.matricant_product(M)     
        
        rt = Acoustique_Matricant.RT_fluide(M_tot,omega,lenteur,e,rho_ext,celerite)
        T = np.abs(np.squeeze(rt.T))**2
        
        R_sq , T_sq = RT_coefficient_periodic_analytic_duralumin_epoxy(np.squeeze(omega),N_cell)
    
        afficher_figure = Acoustique_Matricant.Ploter()
        xy = [[[freq[1:],T],[freq[1:],T_sq]]]
        fontsize = 20
        figsize=(9,9)
        title = 'Transmission & Réflexion théorique (formule Sasha) vs calculs'
        xy_labels = [['Fréquence (MHz)','|T|'],['Fréquence (MHz)','|T|']]
        
        labels = [[['Code de calcul'],["Modèle Théorique"]]]
        extent = [0,1,0,freq[-1]]
        cmap = 0
        
        afficher_figure.multi_plot(xy, 1,1,figsize,title,labels,extent,cmap,xy_labels,fontsize)
        
    elif test == 2:
        N_cell = 6
        
        with open('deck.yaml','r') as fichier:
            data_materiaux = yaml.safe_load(fichier)
        
        constantes_physique = Acoustique_Matricant.Proprietes_Physiques(data_materiaux,'isotrope')
        dictionnaire = constantes_physique.dic
        
        rho , cij , e = constantes_physique.periodic_multilayer(N_cell)
        
        with open('deck_ext.yaml','r') as fichier:
            data_exterieur = yaml.safe_load(fichier)
        
        rho_ext = data_exterieur['1']['rho']
        celerite = data_exterieur['1']['c']
    
        freq       = np.linspace(0.,10.,num=201,dtype=float)
        omega      = freq[1:]*(2*np.pi)
        N_omega = omega.shape[0]
        
        theta = np.linspace(0.,0.35,num=101,dtype=float)
        lenteur = np.sin(theta)/celerite
        
        calculs_matricant = Acoustique_Matricant.Matricant(cij, e, rho, lenteur, omega)
       
        N = calculs_matricant.stroh_4blocs()
        
        calculs_matricant.reshape()
        Q = calculs_matricant.matrix_Q()
        M = calculs_matricant.matricant(Q)
        M_tot=calculs_matricant.matricant_product(M)     
        
        rt = Acoustique_Matricant.RT_fluide(M_tot,omega,lenteur,e,rho_ext,celerite)
        T = np.abs(np.squeeze(rt.T))**2
                
        max_f = np.max(freq)
        min_f = np.min(freq)
        max_theta = np.max(theta*180/np.pi)
        min_theta = np.min(theta*180/np.pi)
        
        
        afficher_figure = Acoustique_Matricant.Ploter()
        T = np.swapaxes(T, 0, 1)
        xy=[[[T]]]
        fontsize = 20
        figsize=(9,9)
        title = f"Transmission fonction des angles d'incidence Eau_[Aluminium][Plexi]_Eau pour {N_cell} cellule(s)"
        xy_labels = [['Fréquence (MHz)','|T|'],['Fréquence (MHz)','|T|']]
        
        labels = [['T']]
        extent = [min_f,max_f,min_theta,max_theta]
        cmap = 'jet'
        afficher_figure.multi_plot(xy, 1,1,figsize,title,labels,extent,cmap,xy_labels,fontsize)
