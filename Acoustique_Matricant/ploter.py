# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 15:41:05 2024

@author: Cyril
"""

import numpy as np
import matplotlib.pyplot as plt

class Ploter():

    
    def multi_plot(self,xydata,ncols,nrows,figsize,title,labels,extent,cmap,xy_labels,fontsize):
        #title : chaine de caractères qui est le titre de tous les graphes
        #lables : titre de chaque graphe. Syntaxe : [['T_th','T_exp'],['R_th','R_exp']] --> 2 graphes avec 2 courbes
        #Syntaxe xydata : xy=[[[freq[1:],T_th],[freq[1:],T_f]],[[freq[1:],R_th],[freq[1:],R_f]]] --> 2 graphes avec 2 courbes
        fig,axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=figsize)
        try :
            iter(axes)
        except :
            axes=np.asarray(axes)
        for axe , xydata_in_axe , label , xy_label in zip(axes.flat,xydata,labels,xy_labels):
            for xy in xydata_in_axe:
                if np.shape(xy[0][1])==():
                    #on trace des courbes
                    axe.plot(xy[0],xy[1])
                    axe.set_xlabel(xy_label[0],fontsize=fontsize)
                    axe.set_ylabel(xy_label[1],fontsize=fontsize)
                    # print(xy[0])
                    axe.tick_params(axis='both', which='major', labelsize=25)
                    axe.legend(labels=label,fontsize=fontsize*2/3,loc='upper right') 
                else :
                    #on trace la matrice
                    im = axe.imshow(xy[0], extent=extent, aspect='auto', origin='lower',  cmap=cmap)
                    cbar = fig.colorbar(im, ax=axe)
                    axe.set_xlabel(xy_label[0],fontsize=fontsize)
                    axe.set_ylabel(xy_label[1],fontsize=fontsize)
                    axe.tick_params(labelsize=25)
                    cbar.ax.tick_params(labelsize=fontsize)
        plt.suptitle(title,fontsize=fontsize)
        fig.tight_layout()
        plt.show()
        pass
