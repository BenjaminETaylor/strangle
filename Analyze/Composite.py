# -*- coding: utf-8 -*-
"""
===============================================================================
    StrEngL.Analyze.Composite
    Copyright (C) 2016  Benjamin E. Taylor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
===============================================================================
"""

#from __future__ import division
import numpy as np
import pylab
from StrEngL.Composite.Laminate import _T_sig, _T_eps

def _repeatIndices(List):
    """creates a new List from given List with each term repeated except for
    the first and last indices"""
    newList = [List[0]]
    for item in List[1:-1]:
        newList.append(item)
        newList.append(item)
    newList.append(List[-1])
    return newList


class Analysis():
    """defines an analysis object"""
    def __init__(self, laminate, loads):
        """laminate = a laminate object"""
        self.lam = laminate
        self.loads = np.array(loads)

#    def lamLoads(self, strains):
#        """returns laminate running loads and moments resulting from given
#        applied strains ([[ex], [ey], [exy], [kx], [ky], [kxy]])"""
#        return self.lam.ABD.dot(strains)

    def lamMidStrains(self):
        """returns laminate mid-plane strains and curvatures resulting from
        given applied loads ([[Nx], [Ny], [Nxy], [Mx], [My], [Mxy]])"""
        return self.lam.abd.dot(self.loads)

    def globalStrain(self, z):
        """returns the laminate strains resulting from given applied loads
        at location z (thru the thickness)"""
        strains = self.lamMidStrains()
        return strains[:3] + strains[3:] * z

    def globalPlyStress(self, plyNumber, strain):
        """returns ply stress for plyNumber under given strain"""
        Qbar = self.lam.plies[plyNumber].Qbar()
        return Qbar.dot(strain)
        
    def globalStrainList(self):
        """returns a List of ply boundary strains from bottom to top of
        laminate (useful for graphing laminate strains)"""
        strains = []
        for z in self.lam.z:
             strains.append(self.globalStrain(z))
        return strains

    def globalStressList(self, component=0):
        """returns a List of global ply stresses from bottom to top of
        laminate (useful for graphing laminate stresses)"""
        stresses = []
        for k in range(len(self.lam.plies)):
            stresses.append(self.globalPlyBotStress(k))
            stresses.append(self.globalPlyTopStress(k))
        return [ply[component] for ply in stresses]

    def globalPlyMidstrain(self, plyNumber):
        """returns the global ply strains at the mid ply for plyNumber
        (starting from 0)"""
        z = (self.lam.z[plyNumber] + self.lam.z[plyNumber+1]) / 2
        return self.globalStrain(z)

    def globalPlyMidStress(self, plyNumber):
        """returns the global ply stresses at mid ply for plyNumber
        (starting from 0)"""
        strain = self.globalPlyMidstrain(plyNumber)
        return self.globalPlyStress(plyNumber, strain)

    def globalPlyBotStress(self, plyNumber):
        """returns the global ply stresses at bottom of plyNumber
        (starting from 0)"""
        z = self.lam.z[plyNumber]
        strain = self.globalStrain(z)
        return self.globalPlyStress(plyNumber, strain)

    def globalPlyTopStress(self, plyNumber):
        """returns the global ply stresses at top of plyNumber
        (starting from 0)"""
        z = self.lam.z[plyNumber+1]
        strain = self.globalStrain(z)
        return self.globalPlyStress(plyNumber, strain)

    def localPlyMidStresses(self):
        """returns a List of local ply stresses, beginning with first ply,
        from middle of each ply (useful for mid-ply stress analysis)"""
        localStressList = []
        for k,ply in enumerate(self.lam.plies):
            midStress = self.globalPlyMidStress(k)
            locStress = _T_sig(ply.angle).dot(midStress)
            localStressList.append([locStress])
        return localStressList

    def localPlyMidStrains(self):
        """returns a List of local ply strains, beginning with first ply,
        from middle of each ply (useful for mid-ply strain analysis)"""
        localStrainList = []
        for k,ply in enumerate(self.lam.plies):
            midStrain = self.globalPlyMidstrain(k)
            localStrain = _T_eps(ply.angle).dot(midStrain)
            localStrainList.append([localStrain])
        return localStrainList
        
    def localPlyTopBotStresses(self):
        """returns a List of local ply stresses, beginning with the bottom
        of the first ply, from the bottom and top of each ply (useful for
        top-bottom stress analysis)"""
        localTopBotStresses = []
        for k,ply in enumerate(self.lam.plies):
            botStress = self.globalPlyBotStress(k)
            localBotStress = _T_sig(ply.angle).dot(botStress)
            topStress = self.globalPlyTopStress(k)
            localTopStress = _T_sig(ply.angle).dot(topStress)
            localTopBotStresses.append([localBotStress, localTopStress])
        return localTopBotStresses
        
    def localPlyTopBotStrains(self):
        """returns a List of local ply strains, beginning with the bottom
        of the first ply, from the bottom and top of each ply (useful for
        top-bottom strain analysis)"""
        localTopBotStrains = []
        strains = self.globalStrainList()
        for k,ply in enumerate(self.lam.plies):
            localBotStrain = _T_eps(ply.angle).dot(strains[k])
            localTopStrain = _T_eps(ply.angle).dot(strains[k+1])
            localTopBotStrains.append([localBotStrain, localTopStrain])
        return localTopBotStrains

    def plotGlobalStress(self, component=0):
        """plots global laminate stress component"""
        pylab.plot(self.globalStressList(component),
                   _repeatIndices(self.lam.z))
        pylab.xlabel('stress')
        pylab.ylabel('thru thickness')
        pylab.title('Laminate Stresses')
        pylab.show()

    def plotGlobalStrain(self, component=0):
        """plots global laminate strain component"""
        strains = [strain[component] for strain in self.globalStrainList()]
        pylab.plot(strains,self.lam.z)
        pylab.xlabel('strain')
        pylab.ylabel('thru thickness')
        pylab.title('Laminate Strains')
        pylab.show()


class Thermal(Analysis):
    """defines a thermal analysis"""
    def __init__(self, lam, loads, deltaT=0):
        Analysis.__init__(self, lam, loads)
        self.deltaT = deltaT
        self.thermLoads = np.array(np.bmat([[self.lam.Nt(self.deltaT)],
                                            [self.lam.Mt(self.deltaT)]]))

    def lamMidStrains(self):
        """returns laminate mid-plane strains and curvatures resulting from
        given applied loads (includes thermal loads)"""
        return self.lam.abd().dot(self.loads + self.thermLoads)

    def _plyStress(self, plyNumber, strain):
        """returns ply stress for plyNumber under given strain (considering
        thermal strains)"""
        ply = self.lam.plies[plyNumber]
        thermStrain = ply.alphaBar() * self.deltaT
        return ply.Qbar().dot(strain - thermStrain)

class Strength(Analysis):
    """defines a strength analysis"""
    def minMargin(self, FIs):
        """
        FIs == array of arrays full of failure indices for each ply
        """
        return 1/np.max(FIs)-1

class MaxStrain(Strength):
    """defines a maximum strain analysis"""
    def analyze(self, midply=False):
        """analyzes each ply in the laminate for strain failure"""
        FI = []
        
        if midply:
            plyStrains = self.localPlyMidStrains()
        else:
            plyStrains = self.localPlyTopBotStrains()
        
        for k,plies in enumerate(self.lam.plies):
            for strain in plyStrains[k]:
                if strain[0] > 0:
                    FI1 = strain[0] / self.lam.plies[k].mat.e1t
                else:
                    FI1 = strain[0] / self.lam.plies[k].mat.e1c
                if strain[1] > 0:
                    FI2 = strain[1] / self.lam.plies[k].mat.e2t
                else:
                    FI2 = strain[1] / self.lam.plies[k].mat.e2c
                FI12 = abs(strain[2]) / self.lam.plies[k].mat.e12
                FI.append([FI1,FI2,FI12])
                
        return np.array(FI)


class MaxStress(Strength):
    """defines a maximum stress analysis"""
    def analyze(self, midply=False):
        """analyzes each ply in the laminate for stress failure at mid-ply"""
        FI = []
        
        if midply:
            plyStresses = self.localPlyMidStresses()
        else:
            plyStresses = self.localPlyTopBotStresses()
                    
        for k,plies in enumerate(self.lam.plies):
            for stress in plyStresses[k]:
                if stress[0] > 0:
                    FI1 = stress[0] / self.lam.plies[k].mat.Xt
                else:
                    FI1 = stress[0] / self.lam.plies[k].mat.Xc
                if stress[1] > 0:
                    FI2 = stress[1] / self.lam.plies[k].mat.Yt
                else:
                    FI2 = stress[1] / self.lam.plies[k].mat.Yc
                FI12 = abs(stress[2]) / self.lam.plies[k].mat.S12
                FI.append([FI1,FI2,FI12])
                
        return np.array(FI)

class TsaiWu(Strength):
    """defines a Tsai-Wu analysis"""
    def analyze(self, midply=False):
        """analyzes each ply in the laminate for failure using
        stress based Tsai-Wu failure criterion"""
        FI = []
        
        if midply:
            plyStresses = self.localPlyMidStresses()
        else:
            plyStresses = self.localPlyTopBotStresses()
        
        for k,ply in enumerate(self.lam.plies):
            F1 = 1/self.lam.plies[k].mat.Xt + 1/self.lam.plies[k].mat.Xc
            F2 = 1/self.lam.plies[k].mat.Yt + 1/self.lam.plies[k].mat.Yc
            F11 = -1/(self.lam.plies[k].mat.Xt * self.lam.plies[k].mat.Xc)
            F22 = -1/(self.lam.plies[k].mat.Yt * self.lam.plies[k].mat.Yc)
            F66 = 1/self.lam.plies[k].mat.S12
            F12 = 1/(2 * self.lam.plies[k].mat.Xt * self.lam.plies[k].mat.Xc)
            
            for stress in plyStresses[k]:
                FIp = F1*stress[0] + F2*stress[1] + F11*stress[0]**2 + \
                      F22*stress[1]**2 + 2*F12*stress[0]*stress[1] + \
                      F66*stress[2]**2
                FI.append([FIp])
            
        return np.array(FI)