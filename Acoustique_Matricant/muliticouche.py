#Classe pour définir les propriétés du/des matériaux

import numpy as np

class Proprietes_Physiques():
    
    def __init__(self,dic,isotropic):
        self.dic = dic
        if isotropic == 'isotrope':
            self.dic = self.dic_isotrope()
    
    #Il faudrait ajouter les axes cristallographiques de chaque couche pour faire les rotations
    def unit_cell(self,interieur=True):
        """
        Arrays containing the physical properties of the unit cell 

        INPUT
        ----------
        dic : {'Unit cell' : ('1','2', ...) tuple of keys that represents the material layering,           
                          '1': dict of the material properties of material n°1 {'Cij':...,'rho':...,'thickness':...}
                          '2': dict of the material properties of material n°2 {'Cij':...,'rho':...,'thickness':...}
                          ...}
        
        OUTPUT
        -------
        rho : Array (N_layerx1)  of the density of the layers of the unit cell
        C   : Array (N_layerx6x6)of the elasticity matric/tensor (i,j) of the layers of the unit cell
        e   : Array (N_layerx1)  of the thickness of the layers of the unit cell
        """
        C = [self.dic[mat]['Cij'] for mat in self.dic['Unit cell']]
        rho = [self.dic[mat]['rho'] for mat in self.dic['Unit cell']]
        if interieur:
            e = [self.dic[mat]['thickness'] for mat in self.dic['Unit cell']]
            return rho , C , e
        else :
            return rho , C
        
        
    def dic_isotrope(self):
        """
        dictionnaire avec les matrices de rigidités à la place de Nu        

        Parameters
        ----------
        dic : dictionnaire
            

        Returns
        -------
        """
        E = np.asanyarray([self.dic[mat]['E'] for mat in self.dic['Unit cell']])
        Nu = np.asarray([self.dic[mat]["Nu"] for mat in self.dic['Unit cell']])
        nbr_couches = E.shape[0]

        c11 = (1.-Nu)*E/((1.+Nu)*(1.-2.*Nu))
        c12 = Nu*E/((1.+Nu)*(1.-2.*Nu))

        cij = np.zeros((nbr_couches,6,6))
        c11.shape = (nbr_couches,1,1)
        c12.shape = (nbr_couches,1,1)
        cij[:,:3,:3] = c12
        np.squeeze(c11)
        np.squeeze(c12)
        
        for i in range(nbr_couches):
            np.fill_diagonal(cij[i,:3,:3], c11[i])
            np.fill_diagonal(cij[i,3:,3:], (c11[i]-c12[i])/2)
            del self.dic[str(i+1)]['E'],self.dic[str(i+1)]['Nu']
            self.dic[str(i+1)]['Cij'] = cij[i].tolist()
        return self.dic
    
    def milieu_exterieur_fluide(self):
        c = [self.dic[mat]['c'] for mat in self.dic['Unit cell']]
        rho = [self.dic[mat]['rho'] for mat in self.dic['Unit cell']]
        return c,rho
    
    
    def periodic_multilayer(self,N_cell):
        """
        Arrays containing the physical properties of the periodic multilayer

        INPUT
        ----------
        dic_unit_cell : {'Unit cell' : ('1','2', ...) tuple of keys that represents the material layering,           
                          '1': dict of the material properties of material n°1 {'Cij':...,'rho':...,'thickness':...}
                          '2': dict of the material properties of material n°2 {'Cij':...,'rho':...,'thickness':...}
                          ...}
               }
        N_cell : number of unit cells 
        
        OUTPUT
        -------
        rho : list of the density of the layers of the multilayer
        C   : list of the elasticity matric/tensor (i,j) of the layers of the unit cell
        e   : list of the thickness of the layers of the multilayer
        
        """
        rho_unit_cell,C_unit_cell,e_unit_cell = self.unit_cell()
        rho = np.array(N_cell*rho_unit_cell)
        C = np.array(N_cell*C_unit_cell)
        e = np.array(N_cell*e_unit_cell)
        return rho,C,e
    