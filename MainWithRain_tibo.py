import numpy as np
import random
import matplotlib.pyplot as plt
from copy import copy
import pandas as pd
import itertools
import seaborn as sns
from scipy.ndimage import convolve, generate_binary_structure
from itertools import product
from matplotlib import colors
import matplotlib.patches as mpatches
from matplotlib import cm
from CellObject import Cell
from PIL import Image
import sys


# TODO, Qds (surface water detension in ht and Qe0 calculation)
# I (infiltration rate)
# k (decay rate)
# Rain
# Tune parameters (hv, temp)

dfElevation = pd.read_csv('Data/MNT_ascii_Oued_BoukhalefCSV.csv',header=0,sep=';|:|,',engine='python')
dfElevation = pd.DataFrame(dfElevation)

dfLandUse = pd.read_csv('Data/Landuse_Oued_BoukhalefCSV.csv',header=0,sep=';|:|,',engine='python')
dfLandUse = pd.DataFrame(dfLandUse)

# cmap = colors.ListedColormap(['white','red',"magenta","green","cyan","orange","brown","yellow","black","grey",'blue'])
# bounds=[0,1,2,3,4,5,6,7,8,9,10,11]
# norm = colors.BoundaryNorm(bounds, cmap.N)
# plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))

# tell imshow about color map so that only set colors are used
# img1 = plt.imshow(dfLandUse,cmap=cmap, norm=norm)

#print(dfElevation)
#dfElevation = dfElevation.replace(-9999,-1)
#print(dfElevation)
#img2 = plt.imshow(dfElevation, vmin=-100, vmax=400)
#plt.show()


npCells = dfElevation.to_numpy()  #numpy array filled with cell objects
npCells = np.empty_like(npCells,dtype=Cell)
nSquares = 0
#print(npCells.shape)
for i in range(npCells.shape[0]):
    for j in range(npCells.shape[1]):
        if (dfElevation.iloc[i,j] != -9999) and (dfLandUse.iloc[i,j] != -9999):  #only create cell if it is in the region of data
            npCells[i,j] = Cell(i,j,dfElevation.iloc[i,j], dfLandUse.iloc[i,j])
            nSquares += 1

#prints to see if there are -9999's inside region
#for i in range(195,205,1):
#    for j in range(50,100,1):
#        print('i,j,elevation:', i, j, npCells[i,j].elevation)
#DOES NOT SHOW ANY ISSUES

#print('Number of cells:', nSquares)
#print(npCells)
#print(npCells[150,150].x,npCells[150,150].y,npCells[150,150].elevation)
#print(npCells[150,150].cLeft)


nRows = npCells.shape[0]
nCols = npCells.shape[1]
for i in range(nRows):
    for j in range(nCols):
        if npCells[i,j] != None:  #only set neighbors if cell exist there
            cList = []
            for n in [[-1,0],[0,1],[1,0],[0,-1]]:
                if (i+n[0] < 0) or (j+n[1] < 0) or (i+n[0] > nRows - 1) or (j+n[1] > nCols - 1):
                    c = None
                else:
                    c = npCells[i+n[0],j+n[1]]
                cList.append(c)
            if i == 0 and j == 0:
                print(cList)
            npCells[i,j].setNeighbors(cList[0],cList[1],cList[2],cList[3])


# Erik test script--------------------------------------------

# Initialization
T = int(144*4)# No. time steps (144 = 24h)
p = 0.1*(24+25+27+29+31+32+31+30+28+26+24+12)/(12*100) # Mean amount of sunlight in Morocco in a year ref paper tibo
temp = 20
Qr = 0
Qe = 0
Qi = 0
Qps = 0
Qpi = 0
Qvi = 0
Qvs = 0
dt = 10*60  #delta t is 10 minutes as per paper

# II = "forest"
II = "settlement"
#II = "culture"
#II = "Real_Land_Use"



waterTot = 54631*Qps

# Plot data as a function of t
Qr_tot = np.zeros((T+1))
Qpi_tot = Qr_tot.copy()
Qe_tot = Qr_tot.copy()
Qi_tot = Qr_tot.copy()
Qps_tot = Qr_tot.copy()
Qvs_tot = Qr_tot.copy()
Qvi_tot = Qr_tot.copy()


rainRegion_iMin = 0
rainRegion_iMax = nRows
rainRegion_jMin = 0
rainRegion_jMax = nCols
rainIntensity = 8  #mm/h  (Slight: < 0.5, Moderate: 0.5 ~ 4, Heavy : 4 ~ 8, Very heavy: > 8) - see also definitions of showers in link:
#https://water.usgs.gov/edu/activity-howmuchrain-metric.html#:~:text=Moderate%20rain%3A%20Greater%20than%200.5,than%202%20mm%20per%20hour.
rainTimeStart = 0
rainTimeStop = 0.5*T  #let it rain for first half of simulation period


#to visualize water distribution on surface
npQpsValues = np.empty_like(npCells) #Qps set to NaN if no cell exist
npQpiValues = np.empty_like(npCells)  
npStates = np.empty_like(npCells)

for i in range(nRows):
    for j in range(nCols):
        if npCells[i,j] != None:
            npCells[i,j].setQps(Qps)
            npCells[i,j].setQr(Qr)
            npCells[i,j].setQe(Qe)
            npCells[i,j].setQpi(Qpi)
            npCells[i,j].setQi(Qi)
            npCells[i,j].setQvi(Qvi)
            npCells[i,j].setQvs(Qvs)
            
            npQpsValues[i,j] = Qps
            npQpiValues[i,j] = Qpi
            #npCells[i,j].set_Kc_r(npCells[i,j].landUse)
            npCells[i,j].set_Kc_r_I_Qs_settlement(npCells[i,j].landUse)
            npCells[i,j].setQds()
            #npCells[i,j].setQs()

#print('Qps values before updates:\n', npQpsValues)
npQpsValues = npQpsValues.astype(float) / (25*25) * 1000  #convert from Qps to mm on surface
npQpiValues = npQpiValues.astype(float) / (25*25) * 1000  #convert from Qpi to mm on surface
cmap = colors.LinearSegmentedColormap.from_list('custom blue', ['#b2bac2','#1f6dc4'], N=256)
plt.matshow(npQpsValues,cmap=cmap)
plt.title('t = 0')
plt.colorbar()
plt.clim(0,120)
# plt.savefig('plots/testInfEva/rainfall/t=0.pdf',format='pdf')
# plt.savefig('plotsQps/t=' + str((0)*dt) + '.jpeg',format='jpeg')
#plt.show()
plt.clf()
plt.close()


# Simulation
for t in range(T):

    # Qe_temp for t+1
    for i in range(nRows):
        for j in range(nCols):
            if npCells[i,j] != None:
                npCells[i,j].set_ht()

    for i in range(nRows):
        for j in range(nCols):
            if npCells[i,j] != None:
                npCells[i,j].set_all_p()
    # Qe
    for i in range(nRows):
        for j in range(nCols):
            if npCells[i,j] != None:
                Qe0 = npCells[i,j].Qps + npCells[i,j].Qr - npCells[i,j].Qi - npCells[i,j].Qvs # - npCells[i,j].Qds
                h = Qe0*npCells[i,j].hv
                npCells[i,j].set_all_v(h) # We understood v was calculated only for points s belonging to E (lower down)
                Qe1 = (npCells[i,j].vUp + npCells[i,j].vRight +  npCells[i,j].vDown +  npCells[i,j].vLeft) * h * npCells[i,j].hv

                if Qe0 > 0:
                    d = min((npCells[i,j].r / dt) * Qe0,Qe1)
                else:
                    d = 0

                npCells[i,j].setQe_temp(d*dt)
                
    # Qr_temp for t+1 (WITH OPTIONAL RAIN as function of area and time)
    for i in range(nRows):
        for j in range(nCols):
            if npCells[i,j] != None:
                npCells[i,j].set_all_lambda()
                npCells[i,j].setQr_s()

                #adding rain
                if (i >= rainRegion_iMin and
                    i <= rainRegion_iMax and
                    j >= rainRegion_jMin and
                    j <= rainRegion_jMax and
                    t >= rainTimeStart and
                    t <= rainTimeStop):
                    npCells[i,j].setQr_temp(rainIntensity / 60 * 10 * (25*25) / 1000)  #mm/h to m^3 / 10min as per Qps units

    for i in range(nRows):
        for j in range(nCols):
            if npCells[i,j] != None:
                # Boundary conditions for getting Qr from neighbors,
                #changed to match new structure of no cells outside region
                try:
                    npCells[i,j].setQr_temp(npCells[i-1,j].QrDown)
                except:
                    pass
                try:
                    npCells[i,j].setQr_temp(npCells[i+1,j].QrUp)
                except:
                    pass
                try:
                    npCells[i,j].setQr_temp(npCells[i,j-1].QrRight)
                except:
                    pass
                try:
                    npCells[i,j].setQr_temp(npCells[i,j+1].QrLeft)
                except:
                    pass

    # Evaporated Water
    for i in range(nRows):
        for j in range(nCols):
            if npCells[i,j] != None:
                Ep = npCells[i,j].Kc*p*(0.46*temp+8.13)
                #Ep = 0  #JUST TO TEST
                Qv0 = npCells[i,j].Qps + npCells[i,j].Qr
                Qv1 = np.minimum(npCells[i,j].Qpi, (Ep-Qv0))
                if Qv0 >= Ep:
                    npCells[i,j].setQvi(0)
                    npCells[i,j].setQvs(Ep)
                else:
                    npCells[i,j].setQvi(Qv1)
                    npCells[i,j].setQvs(Qv0)
                #if i == 150:
                 #   print('\nEp:', Ep, '\nQps:',  npCells[i,j].Qps, '\nQr:',  npCells[i,j].Qr, '\nQvs:', npCells[i,j].Qvs_temp,  '\nQr:',  npCells[i,j].Qr, '\nQv0:', Qv0)
                #if t < 2 and i == 150:
                 #   print(Ep, Qv0, Qv1, npCells[i,j].Qpi, npCells[i,j].Qr, npCells[i,j].Qps)

    # Infiltrated Water
    for i in range(nRows):
        for j in range(nCols):
            if npCells[i,j] != None:
                Qi0 = npCells[i,j].Qps + npCells[i,j].Qr - npCells[i,j].Qvs
                Qi1 = min(Qi0,npCells[i,j].I*dt)
                Qi2 = npCells[i,j].Qs - npCells[i,j].Qpi + npCells[i,j].Qvi
                npCells[i,j].setQi_temp(min(Qi1,Qi2))
                #if t < 3 and i == 150:
                 #   print('\nEp:', Ep, '\nQi0:', Qi0, '\nQi1:', Qi1, '\nQi2:', Qi2, '\nQps:',  npCells[i,j].Qps, '\nQr:',  npCells[i,j].Qr, '\nQvs:',  npCells[i,j].Qvs)
                    #print('\nQi0:', Qi0, '\nQi1:', Qi1, '\nQi2:', Qi2, '\nQi_temp:', npCells[i,j].Qi_temp)
    
    # Updating dynamics: Qe, Qr, Qi, Qvi, Qvs, Qpi, Qps for t+1
    for i in range(nRows):
        for j in range(nCols):
            if npCells[i,j] != None:
                npCells[i,j].setQr(npCells[i,j].Qr_temp)
                npCells[i,j].setQe(npCells[i,j].Qe_temp)
                npCells[i,j].setQi(npCells[i,j].Qi_temp)

                npCells[i,j].Qpitemp =npCells[i,j].Qpi
                npCells[i,j].Qpstemp =npCells[i,j].Qps
                
                npCells[i,j].setQpi(npCells[i,j].Qpi + npCells[i,j].Qi - npCells[i,j].Qvi)
                npCells[i,j].setQps(npCells[i,j].Qps + npCells[i,j].Qr - npCells[i,j].Qe - npCells[i,j].Qi - npCells[i,j].Qvs)
                #test prints
                #if t % 10 == 0 and i == 150:
                 #   print('Ep:', Ep, '\nQpi:', npCells[i,j].Qpi,  '\nQr:', npCells[i,j].Qr,  '\nQps:', npCells[i,j].Qps,  '\nQi:', npCells[i,j].Qi, '\nQvs:', npCells[i,j].Qvs,  '\nQvi:', npCells[i,j].Qvi, '\nQe:', npCells[i,j].Qe, '\nr:', npCells[i,j].r)
                #if t > 1:
                 #   print(npCells[i,j].Qvi_temp, npCells[i,j].Qvs_temp)
                #if t > 1 and npCells[i,j].Qv != None:
                #    print(npCells[i,j].Qi, npCells[i,j].Qv)
                npCells[i,j].resetQr_temp()
                
                # Plot water quantities
                Qr_tot[t+1] += npCells[i,j].Qr
                Qpi_tot[t+1] += npCells[i,j].Qpitemp
                Qe_tot[t+1] += npCells[i,j].Qe
                Qi_tot[t+1] += npCells[i,j].Qi
                Qps_tot[t+1] += npCells[i,j].Qpstemp
                Qvs_tot[t+1] += npCells[i,j].Qvs
                Qvi_tot[t+1] += npCells[i,j].Qvi

    # To see if we ever reach state s=1 : saturation with no water on surface.
    ls0 = []
    ls1 = []
    ls2 = []
    ls3 = []
    
    for i in range(nRows):
        for j in range(nCols):
            if npCells[i,j] != None:
                npQpsValues[i,j] = npCells[i,j].Qps
                npQpiValues[i,j] = npCells[i,j].Qpi
                npCells[i,j].setState()
                npStates[i,j] = npCells[i,j].state
                if npStates[i,j] == 0 :
                    ls0.append(npStates[i,j])
                elif npStates[i,j] == 1 :
                    ls1.append(npStates[i,j])
                elif npStates[i,j] == 2 :
                    ls2.append(npStates[i,j])
                elif npStates[i,j] == 3 :
                    ls3.append(npStates[i,j])
                if i > 0 and j == 240 and t % 10 == 0:
                    #print('Qps:', npCells[i,j].Qps, '\tElevation:', npCells[i,j].elevation)  #print to see if water flow in the right directions
                    print('\nQr:', npCells[i,j].Qr, '\nQe:', npCells[i,j].Qe, '\nQps:', npCells[i,j].Qps, '\nQi:', npCells[i,j].Qi, '\nQvi:', npCells[i,j].Qvi, '\nQvs:', npCells[i,j].Qvs,'\nQpi:', npCells[i,j].Qpi)
                    #print('selfElevation:', npCells[i,j].elevation, '\tupNeighborElevation:', npCells[i,j].cUp.elevation, '\tself_ht:', npCells[i,j].ht, '\tupNeighbor_ht:', npCells[i,j].cUp.ht)
                    #print('Lambdas:',npCells[i,j].lambdaUp,npCells[i,j].lambdaDown,npCells[i,j].lambdaRight,npCells[i,j].lambdaLeft)

    if (t+1) % 1 == 0:
        print('t = ' + str(t+1))
        # # print('Qps values after updates:\n', npQpsValues)
        # plt.figure(2)
        npQpsValues = npQpsValues.astype(float) / (25*25) * 1000
        # plt.matshow(npQpsValues,cmap=cmap)
        # plt.title('Water on the surface ('+II+'), t = ' + str((t+1)*dt))
        # plt.colorbar()
        # plt.clim(0,120)
        # plt.savefig('plotsQps_folder/plotsQps_'+ II +'/'+II+'_t=' + str((t+1)*dt) + '.jpeg',format='jpeg')
        # # plt.show()
        # plt.close()
        # plt.clf()
        
        # plt.figure(3)
        npQpiValues = npQpiValues.astype(float) / (25*25) * 1000
        # plt.matshow(npQpiValues,cmap=cmap)
        # plt.title('Water in the soil ('+II+'), t = ' + str((t+1)*dt))
        # plt.colorbar()
        # plt.clim(0,120)
        # plt.savefig('plotsQpi_folder/plotsQpi_' + II +'/'+II+'_t=' + str((t+1)*dt) + '.jpeg',format='jpeg')
        # # plt.show()
        # plt.close()
        # plt.clf()
        
        #Code for plotting the states of map.
        plt.figure(4)
        bounds2=[0,1,2,3]
        cmap2 = colors.ListedColormap(['red',"magenta","green","orange"])
        norm2 = colors.BoundaryNorm(bounds2, cmap2.N)
        plt.colorbar(cm.ScalarMappable(norm=norm2, cmap=cmap2))
        # Creating legend with color box
        red_patch = mpatches.Patch(color='red', label='No saturation, no water on surface')
        magenta_patch = mpatches.Patch(color='magenta', label='Saturation, no water on surface')
        green_patch = mpatches.Patch(color='green', label='No saturation, water on surface')
        orange_patch = mpatches.Patch(color='orange', label='Saturation and water on surface')
        npStates = npStates.astype(float)
        plt.matshow(npStates,cmap=cmap2)
        plt.legend(handles=[red_patch,magenta_patch,green_patch,orange_patch],loc='upper center')#, bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=2)
        # plt.show()
        plt.title('Map of the cells state ('+II+'),  t = ' + str((t+1)*dt))
        plt.savefig('Plots/States/tibo/'+II+'_t='+ str((t+1)*dt) + '.jpeg',format="jpeg")#,bbox_to_anchor=(0.5, -0.05), bbox_inches='tight')
        plt.close()
        plt.clf()
        


    #for i in range(195,205,1):
     #   for j in range(50,100,1):
      #      if npCells[i,j].Qps < 0.9:
       #         print('i,j,elevation,Qps:', i, j, npCells[i,j].elevation, npCells[i,j].Qps)

print('Qr: ', Qr_tot)
print('Qe: ', Qe_tot)
print('Qps: ', Qps_tot)
print('Qi: ', Qi_tot)
print('Qpi: ', Qpi_tot)
print('Qvi: ', Qvi_tot)
print('Qvs: ', Qvs_tot)

print("(Qps+Qpi+Qvs)/waterTot",(Qps_tot+Qpi_tot+Qvs_tot+Qvi_tot)/waterTot)


fig2,ax2 = plt.subplots()
t_lin = np.linspace(0,T,T+1)
ax2.plot(t_lin, Qr_tot.astype(float), label='Qr')
ax2.plot(t_lin, Qe_tot, label='Qe')
ax2.plot(t_lin, Qi_tot, label='Qi')
ax2.plot(t_lin, Qpi_tot, label='Qpi')
ax2.plot(t_lin, Qps_tot, label='Qps')
ax2.plot(t_lin, Qvi_tot, label='Qvi')
ax2.plot(t_lin, Qvs_tot, label='Qvs')
ax2.plot(t_lin, Qvs_tot+Qvi_tot, label='Qv')
ax2.legend()
plt.show()
plt.xlabel("time")
plt.title("Total Surface, Infiltrated, Evaporated water, "+II)
plt.savefig('Plots/States/tibo/Tot_'+II+'/'+II+'_t=' + str((T)*dt) + '.jpeg',format="jpeg")
plt.close()

#--------------------------------------------------------------------------