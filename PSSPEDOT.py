# define constants
import math
k_B = 8.617333 * 10**(-5) # ev/K
pi = math.pi
hbar = 6.582119569 * 1e-16 # eV*s
'''
change directory and frame number
before use
'''
###### import packages ####################################
import os
import xml.etree.ElementTree as ET
import numpy as np

# plot the data with matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import pylab

import scipy as sp

########################################################
## reorganization energy by orca
Ess_nC = -25419.93474 # eV
Ess_nN = -25420.27921 # eV
Ess_cN = -25411.21219 # eV
Ess_cC = -25411.21265381581 # eV
Re_ss_ss = Ess_nC - Ess_nN + Ess_cN - Ess_cC

Ess_noh_nC = -25405.30295 # eV
Ess_noh_nN = -25406.016087 # eV
Ess_noh_cN = -25401.31412 # eV
Ess_noh_cC = -25401.90196 # eV
Re_ss_ss_noh = Ess_nC - Ess_nN + Ess_noh_cN - Ess_noh_cC
Re_ss_noh_ss = Ess_noh_nC - Ess_noh_nN + Ess_cN - Ess_cC

Re_f = np.zeros(15)
for i in range(15):
    if (i+1)%4==0:
        Re_f[i] = Re_ss_noh_ss
    elif (i+1)%4==3:
        Re_f[i] = Re_ss_ss_noh
    else:
        #ii = ((i+1)//4)*2 + (i+1)%4 -1
        Re_f[i] = Re_ss_ss

Re_b = np.zeros(15)
for i in range(15):
    if (i+1)%4==0:
        Re_b[i] = Re_ss_ss_noh
    elif (i+1)%4==3:
        Re_b[i] = Re_ss_noh_ss
    else:
        #ii = ((i+1)//4)*2 + (i+1)%4 -1
        Re_b[i] = Re_ss_ss

Re_f_NN = np.zeros(14)
for i in range(14):
    if (i+1)%2==1:
        Re_f_NN[i] = Re_ss_ss
    if (i+1)%4==2:
        Re_f_NN[i] = Re_ss_ss_noh
    if (i+1)%4==0:
        Re_f_NN[i] = Re_ss_noh_ss

Re_b_NN = np.zeros(14)
for i in range(14):
    if (i+1)%2==1:
        Re_b_NN[i] = Re_ss_ss
    if (i+1)%4==2:
        Re_b_NN[i] = Re_ss_noh_ss
    if (i+1)%4==0:
        Re_b_NN[i] = Re_ss_ss_noh

Re_f_2NN = np.zeros(13)
for i in range(13):
    if (i+1)%4==1:
        Re_f_2NN[i] = Re_ss_ss_noh
    if (i+1)%4==2 or (i+1)%4==3:
        Re_f_2NN[i] = Re_ss_ss
    if (i+1)%4==0:
        Re_f_2NN[i] = Re_ss_noh_ss

Re_b_2NN = np.zeros(13)
for i in range(13):
    if (i+1)%4==1:
        Re_b_2NN[i] = Re_ss_noh_ss
    if (i+1)%4==2 or (i+1)%4==3:
        Re_b_2NN[i] = Re_ss_ss
    if (i+1)%4==0:
        Re_b_2NN[i] = Re_ss_ss_noh

####### define vaccum values, in eV ##########################
IP_ss_vac = 0.33316 * 27.2116
IP_ssnoh_vac = 0.15117 * 27.2116 # the segment with 19 atoms
EA_ss_vac = 0.00390 * 27.2116
EA_ssnoh_vac = -0.15156 * 27.2116
print('IP of ss_vac: ', 0.33316 * 27.2116)

####### define HOMO/LUMO vaccum values, in eV#################
HOMO_ss_vac = -IP_ss_vac 
HOMO_ss_noh_vac =  -IP_ssnoh_vac 
HOMO_pss_vac = np.zeros(16)
for j in range(16):
    if ((j+1)%4 == 0):
        HOMO_pss_vac[j] = HOMO_ss_noh_vac 
    else:
        HOMO_pss_vac[j] = HOMO_ss_vac
        
LUMO_ss_vac = -EA_ss_vac 
LUMO_ss_noh_vac = -EA_ssnoh_vac
LUMO_pss_vac = np.zeros(16)
for i in range(16):
    if ((i+1)%4 == 0):
        LUMO_pss_vac[i] = -EA_ssnoh_vac 
    else:
        LUMO_pss_vac[i] = -EA_ss_vac 
        
####### define vaccum values for PEDOT, in eV###########
IP_edoth_vac = 0.28616 * 27.2116 # in eV
EA_edoth_vac = -0.04875 * 27.2116
# define HOMO/LUMO vaccum values, in eV
HOMO_edoth_vac = np.zeros(16)
for i in range(16):
    HOMO_edoth_vac[i] = -IP_edoth_vac
LUMO_edoth_vac = np.zeros(16)
for i in range(16):
    LUMO_edoth_vac[i] = -EA_edoth_vac

### class definiation #################################
class PSSPEDOT:
    def __init__(self, filename1, filename2):
        self.filename1 = filename1
        self.filename2 = filename2
    @staticmethod
    def get_E(self):
        tree = ET.parse('%s/ewald_jobs.xml' %self.filename1)
        root = tree.getroot()
        rootLen = len(root) 
        JobNum = int(rootLen/3)
        En = np.zeros(JobNum)
        Ee = np.zeros(JobNum)
        Eh = np.zeros(JobNum)
        self.xcoor = np.zeros(JobNum)
        self.ycoor = np.zeros(JobNum)
        self.zcoor = np.zeros(JobNum)
        for i in range(JobNum):
            En[i] = float(root[3*i][6][0][2].text)
            Ee[i] = float(root[3*i + 1][6][0][2].text)
            Eh[i] = float(root[3*i + 2][6][0][2].text)
            self.xcoor[i] = float(root[3*i][6][0][1].text.split(" ")[0])
            self.ycoor[i] = float(root[3*i][6][0][1].text.split(" ")[1])
            self.zcoor[i] = float(root[3*i][6][0][1].text.split(" ")[2])
            
        # HOMO = HOMO_vac + correction
        self.HOMO_pss = np.zeros(16)
        for j in range(16):
            if ((j+1)%4 == 0):
                self.HOMO_pss[j] = HOMO_pss_vac[j] + (En[j] - Eh[j])
            else:
                self.HOMO_pss[j] =  HOMO_pss_vac[j] + (En[j] - Eh[j])
        self.LUMO_pss = np.zeros(16)
        for j in range(16):
            if ((j+1)%4 == 0):
                self.LUMO_pss[j] = LUMO_pss_vac[j] + (Ee[j] - En[j]) 
            else:
                self.LUMO_pss[j] = LUMO_pss_vac[j] + (Ee[j] - En[j])
                
        self.HOMO_PEDOT = np.zeros(16)
        self.LUMO_PEDOT = np.zeros(16)
        for j in range(16):
            self.HOMO_PEDOT[j] = HOMO_edoth_vac[j] + (En[j+16] - Eh[j+16])
            self.LUMO_PEDOT[j] = LUMO_edoth_vac[j] + (Ee[j+16] - En[j+16])
        
        self.HOMO_pss_del = self.HOMO_pss[1:]-self.HOMO_pss[0:-1]
        self.HOMO_pss_del = np.append(self.HOMO_pss_del, -self.HOMO_pss_del, axis=None)
        self.LUMO_pss_del = self.LUMO_pss[1:]-self.LUMO_pss[0:-1]
        self.LUMO_pss_del = np.append(self.LUMO_pss_del, -self.LUMO_pss_del, axis=None)
        
        self.HOMO_pss_del_NN = self.HOMO_pss[2:]-self.HOMO_pss[0:-2]
        self.HOMO_pss_del_NN =  np.append(self.HOMO_pss_del_NN, -self.HOMO_pss_del_NN, axis=None)
        self.LUMO_pss_del_NN = self.LUMO_pss[2:]-self.LUMO_pss[0:-2]
        self.LUMO_pss_del_NN =  np.append(self.LUMO_pss_del_NN, -self.LUMO_pss_del_NN, axis=None)
        
        self.HOMO_pss_del_2NN = self.HOMO_pss[3:]-self.HOMO_pss[0:-3]
        self.HOMO_pss_del_2NN =  np.append(self.HOMO_pss_del_2NN, -self.HOMO_pss_del_2NN, axis=None)
        self.LUMO_pss_del_2NN = self.LUMO_pss[3:]-self.LUMO_pss[0:-3]
        self.LUMO_pss_del_2NN =  np.append(self.LUMO_pss_del_2NN, -self.LUMO_pss_del_2NN, axis=None)

    @staticmethod    
    def get_J(self):
        tree = ET.parse('%s/iqm_ss_ss.jobs' % self.filename2)
        root = tree.getroot()
        self.J_f = np.zeros(15)
        for i in range(15):
            if (i+1)%4==0 or (i+1)%4==3:
                continue
            else:
                ii = ((i+1)//4)*2 + (i+1)%4 -1
                self.J_f[i] = float(root[ii][6][0][0][0].attrib['j'])
                
        tree = ET.parse('%s/iqm_ss_ss_noh.jobs' % self.filename2)
        root = tree.getroot()
        for i in range(15):
            if (i+1)%4==1 or (i+1)%4==2:
                continue
            else:
                ii = ((i+1)//4)*2 + (4-(i+1)%4)%4 -1  # a mapping of the seg id to job id
                self.J_f[i] = float(root[ii][6][0][0][0].attrib['j'])
        self.J_b = self.J_f
        
        # read overlap integral for NN
        self.J_f_NN = np.zeros(14)
        tree = ET.parse('%s/iqm_ss_ss_NN.jobs' % self.filename2)
        root = tree.getroot()
        for i in range(14):
            if (i+1)%2==0:
                continue
            else:
                ii = (i)/2  # a mapping of the seg id to job id
                ii = int(ii)
                self.J_f_NN[i] = float(root[ii][6][0][0][0].attrib['j'])
        tree = ET.parse('%s/iqm_ss_ss_noh_NN.jobs' % self.filename2)
        root=tree.getroot()
        for i in range(14):
            if (i+1)%2==1:
                continue
            else:
                ii = (i+1)/2-1  # a mapping of the seg id to job id
                ii = int(ii)
                self.J_f_NN[i] = float(root[ii][6][0][0][0].attrib['j'])
        self.J_b_NN =self.J_f_NN
        '''
        # read overlap integral for 2NN
        self.J_f_2NN = np.zeros(13)
        tree = ET.parse('%s/iqm_ss_ss_2NN.jobs' % self.filename2)
        root = tree.getroot()
        for i in range(13):
            if ((i+1)%4==1) or ((i+1)%4==0):
                continue
            else:
                ii = ((i+1)//4) *2 + (i+1)%4-2  # have to understand in a graph
                self.J_f_2NN[i] = float(root[ii][6][0][0][0].attrib['j'])
        tree = ET.parse('%s/iqm_ss_ss_noh_2NN.jobs' % self.filename2)
        root=tree.getroot()
        for i in range(14):
            if ((i+1)%4==2) or ((i+1)%4==3):
                continue
            else:
                ii=(i+1)//4 *2 + ((i+1)%4)%2-1  # have to understand in a graph, +1)%2 is used to tell 0 or 1
                ii = int(ii)
                self.J_f_2NN[i] = float(root[ii][6][0][0][0].attrib['j'])
        self.J_b_2NN =self.J_f_2NN
    '''
    @staticmethod
    def get_Rate_Matrix(self):
        self.rate_f = np.zeros(15)
        self.rate_f = 2*pi/hbar*self.J_f**2/np.sqrt(4*pi*Re_f*k_B*300)\
                *np.exp(-(self.HOMO_pss_del[0:15]-Re_f)**2/(4*Re_f*k_B*300))
        self.rate_b = np.zeros(15)
        self.rate_b = 2*pi/hbar*self.J_b**2/np.sqrt(4*pi*Re_b*k_B*300)\
                *np.exp(-(self.HOMO_pss_del[15:30]-Re_f)**2/(4*Re_b*k_B*300))
            
        self.rate_f_NN = np.zeros(14)
        self.rate_f_NN = 2*pi/hbar*self.J_f_NN**2/np.sqrt(4*pi*Re_f_NN*k_B*300)\
                *np.exp(-(self.HOMO_pss_del_NN[0:14]-Re_f_NN)**2/(4*Re_f_NN*k_B*300))

        self.rate_b_NN = np.zeros(14)
        self.rate_b_NN = 2*pi/hbar*self.J_b_NN**2/np.sqrt(4*pi*Re_b_NN*k_B*300)\
                *np.exp(-(self.HOMO_pss_del_NN[14:28]-Re_f_NN)**2/(4*Re_b_NN*k_B*300))
            
                        
        '''self.rate_f_2NN = np.zeros(13)
        self.rate_f_2NN = 2*pi/hbar*self.J_f_2NN**2/np.sqrt(4*pi*Re_f_2N*k_B*300)\
                *np.exp(-(self.HOMO_pss_del_2NN[0:14]-Re_f_2NN)**2/(4*Re_f_2NN*k_B*300))

        self.rate_b_2NN = np.zeros(15)
        self.rate_b_2NN = 2*pi/hbar*self.J_b_2NN**2/np.sqrt(4*pi*Re_b_2NN*k_B*300)\
                *np.exp(-(self.HOMO_pss_del_2NN[14:28]-Re_f_2NN)**2/(4*Re_b_2NN*k_B*300))
                '''

        self.W_HMBE = np.zeros((17,17))
            # forward rate, lower matrix element
        for i in range(17):
            if i < 15:
                self.W_HMBE[i+1,i] = self.rate_f[i]
            elif i == 15:
                self.W_HMBE[i+1,i] = 1e30
        for i in range(17):
            if i < 14:
                self.W_HMBE[i+2,i] = self.rate_f_NN[i]  
        '''
        for i in range(17):
            if i <13:
                self.W_HMBE[i+3,i] = self.rate_f_2NN[i]
                '''
           
            # backward rate, upper matrix element
        for i in range(17):
            if i < 15:
                self.W_HMBE[i,i+1] = self.rate_b[i]
        for i in range(17):
            if i < 14:
                self.W_HMBE[i,i+2] = self.rate_b_NN[i]
        '''for i in range(17):
            if i <13:
                self.W_HMBE[i,i+3] = self.rate_b_2NN[i]
                '''
        
            # diagonal element
        for i in range(17):
            self.W_HMBE[i,i] = -np.sum(self.W_HMBE[:,i])
