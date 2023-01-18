# Cell Object Class:
# Creates a cell as an object and assigns properties to a cell from functions.
# Details about the mathematics behind the functions are refrenced to in the report.

import numpy as np

class Cell:
    def __init__(self,xPos,yPos,elevation,landUse):
        #intrinsic terrain values
        self.x = xPos
        self.y = yPos
        self.cUp = None  #neighboring cells
        self.cRight = None
        self.cDown = None
        self.cLeft = None
        self.distUp = 25  #|c_ij s(i,j)|: pixels in data are 25x25 m
        self.distRight = 25
        self.distDown = 25
        self.distLeft = 25
        self.elevation = elevation
        self.soilType = None
        self.landUse = landUse

        #calculated intrinsic terrain values
        self.pStarUp = None  #topographic slopes
        self.pStarRight = None
        self.pStarDown = None
        self.pStarLeft = None
        self.I = 0.001  #infiltration intensity
        self.k = 0.001
        self.Qs = 0.1 #self.I / self.k #soil saturation in water, CHECK eq 16 replace 1 with K
        self.Qds = 0  #surface water detention (0 for this project)
        self.r = None  #water flow resistance on the soil
        self.Kc = None

        #dynamic flow values
        self.Qpi = 0  #water quantity present in the soil
        self.Qps = 1  #water quantity present on the surface
        self.ht = None  #water level on the surface
        self.hv = 4*4 / (3*25**2)  #water surface level constant: see eq (20)
        self.pUp = None  #slopes with water level
        self.pRight = None
        self.pDown = None
        self.pLeft = None
        self.vUp = None  #potential flow velocities toward downstream neighboring cells
        self.vRight = None
        self.vDown = None
        self.vLeft = None
        self.lambdaUp = None  #proportionality constants of the slopes to neighbouring cells
        self.lambdaRight = None
        self.lambdaDown = None
        self.lambdaLeft = None
        self.QrUp = None  #receiving water quanitity of neighbouring cells
        self.QrRight = None
        self.QrDown = None
        self.QrLeft = None



        #calculated dynamic flow values
        self.d = None  #rate of runoff
        self.Qv = None  #quantity of evaporated water
        self.Qi = None  #quantity of infiltrated water
        self.Qr = None  #quantity of received water
        self.Qe = None  #quantity of drained water
        self.Qr_temp = 0 #quantity of received water, temporary for synchronous update
        self.Qe_temp = 0 #quantity of drained water, temporary for synchronous update
        self.Qi_temp = 0 #quantity of drained water, temporary for synchronous update
        self.Qvs_temp = 0 #quantity of evaporated water, temporary for synchronous update
        self.Qvi_temp = 0 #quantity of evaporated water, temporary for synchronous update

    def setNeighbors(self,cUp,cRight,cDown,cLeft):
        self.cUp = cUp
        self.cRight = cRight
        self.cDown = cDown
        self.cLeft = cLeft

    def getElevation(self):
        return self.elevation

    def set_ht(self):
        self.ht = max(0,self.hv*(self.Qps - self.Qds))

    def get_ht(self):
        return self.ht
        
    def setQr(self,value):
        self.Qr = value

    def setQr_s(self): # Consider positive lambda values? Positive slopes -> No Qr_s!
        self.QrUp = self.lambdaUp * self.Qe
        self.QrRight = self.lambdaRight * self.Qe
        self.QrDown = self.lambdaDown * self.Qe
        self.QrLeft = self.lambdaLeft * self.Qe
        #if self.QrDown > 1000:
            #print('self.lambdaDown:', self.lambdaDown, 'self.Qe:', self.Qe)

    def setQr_temp(self,value):
        #if value > 1000:
            #print('Value > 1000:', value)
        self.Qr_temp += value

    def resetQr_temp(self):
        self.Qr_temp = 0

    def setQe_temp(self,value):
        self.Qe_temp = value
        
    def setQe(self,value):
        self.Qe = value

    def setQi(self,value):
        self.Qi = value
        
    def setQi_temp(self,value):
        self.Qi_temp = value
    
    def setQpi(self,value):
        self.Qpi = value

    def setQvi_temp(self,value):
        self.Qvi_temp = value

    def setQvi(self,value):
        self.Qvi = value

    def setQvs_temp(self,value):
        self.Qvs_temp = value

    def setQvs(self,value):
        self.Qvs = value

    def setQps(self,value):
        self.Qps = value if value > 0 else 0 

    def set_all_p(self):
        if self.cUp is not None:
            self.pUp = np.arctan((self.cUp.getElevation() + self.cUp.get_ht() - self.elevation - self.ht) / self.distUp)
        if self.cRight is not None:
            self.pRight = np.arctan((self.cRight.getElevation() + self.cRight.get_ht() - self.elevation - self.ht) / self.distRight)
        if self.cDown is not None:
            self.pDown = np.arctan((self.cDown.getElevation() + self.cDown.get_ht() - self.elevation - self.ht) / self.distDown)
        if self.cLeft is not None:
            self.pLeft = np.arctan((self.cLeft.getElevation() + self.cLeft.get_ht() - self.elevation - self.ht) / self.distLeft)

    def set_all_lambda(self):
        pUpAbs = np.abs(self.pUp) if self.pUp is not None else 0
        pRightAbs = np.abs(self.pRight) if self.pRight is not None else 0
        pDownAbs = np.abs(self.pDown) if self.pDown is not None else 0
        pLeftAbs = np.abs(self.pLeft) if self.pLeft is not None else 0

        pAbsSum = 0
        try:  #try / excepts are to avoid error where no neighbors exist. If so, the pass exception just continues
            if self.pUp < 0:
                pAbsSum += pUpAbs
        except:
            pass
        try:
            if self.pRight < 0:
                pAbsSum += pRightAbs
        except:
            pass
        try:
            if self.pDown < 0:
                pAbsSum += pDownAbs
        except:
            pass
        try:
            if self.pLeft < 0:
                pAbsSum += pLeftAbs
        except:
            pass

        try:  #update lambdas
            if self.pUp < 0:
                self.lambdaUp = pUpAbs / pAbsSum
            else:
                self.lambdaUp = 0
        except:
            self.lambdaUp = 0

        try:
            if self.pRight < 0:
                self.lambdaRight = pRightAbs / pAbsSum
            else:
                self.lambdaRight = 0
        except:
            self.lambdaRight = 0

        try:
            if self.pDown < 0:
                self.lambdaDown = pDownAbs / pAbsSum
            else:
                self.lambdaDown = 0
        except:
            self.lambdaDown = 0

        try:
            if self.pLeft < 0:
                self.lambdaLeft = pLeftAbs / pAbsSum
            else:
                self.lambdaLeft = 0
        except:
            self.lambdaLeft = 0


    def set_all_v(self,h):  #set the potential flow velocity for neighbors s in E
        try:  #try / excepts are to avoid error where no neighbors exist. If so, the pass exception just continues
            if self.pUp < 0:
                # h = Qe0*self.hv
                self.vUp = self.r * np.cbrt(h**2) * np.sqrt(np.abs(self.pUp))
            else:
                self.vUp = 0
        except:
            #pass
            self.vUp = 0
        try:
            if self.pRight < 0:
                # h = Qe0*self.hv
                self.vRight = self.r * np.cbrt(h**2) * np.sqrt(np.abs(self.pRight))
            else:
                self.vRight = 0
        except:
            #pass
            self.vRight = 0
        try:
            if self.pDown < 0:
                # h = Qe0*self.hv
                self.vDown = self.r * np.cbrt(h**2) * np.sqrt(np.abs(self.pDown))
            else:
                self.vDown = 0
        except:
            #pass
            self.vDown = 0
        try:
            if self.pLeft < 0:
                # h = Qe0*self.hv
                self.vLeft = self.r * np.cbrt(h**2) * np.sqrt(np.abs(self.pLeft))
            else:
                self.vLeft = 0
        except:
            #pass
            self.vLeft = 0

    def set_Kc_r_I_Qs(self,value):
        
        I_forest = 0.1
        I_culture = 0.01
        I_settlement = 0.001

        if value == 1 : #Fake forest
            self.Kc = 0.5
            self.r = 0.8
            self.I = I_forest
        elif value == 2 : #Infrastructure
            self.Kc = 0.3
            self.r = 0.99
            self.I = I_settlement
        elif value == 3 : #Forest
            self.Kc = 0.5
            self.r = 0.8
            self.I = I_forest
        elif value == 4 : #Vegetation
            self.Kc = 0.75
            self.r = 0.65
            self.I = I_forest
        elif value == 5 : #Culture 1
            self.Kc = 0.65
            self.r = 0.65
            self.I = I_culture
        elif value == 6 : #Culture 2
            self.Kc = 0.65
            self.r = 0.65
            self.I = I_culture
        elif value == 7 : # Soil
            self.Kc = 0.25
            self.r = 0.99
            self.I = I_culture
        elif value == 8 : # Dense forest
            self.Kc = 0.55
            self.r = 0.7
            self.I = I_forest
        elif value == 9 : # Dense vegetation
            self.Kc = 0.75
            self.r = 0.5
            self.I = I_forest
        elif value == 10 : # Building
            self.Kc = 0.3
            self.r = 0.99
            self.I = I_settlement  

        self.Qs = self.I / self.k

    def set_Kc_r_I_Qs_forest(self,value):
        I_forest = 0.1

        self.Kc = 0.5
        self.r = 0.8
        self.I = I_forest

    def set_Kc_r_I_Qs_agriculture(self,value):
        I_culture = 0.01

        self.Kc = 0.65
        self.r = 0.65
        self.I = I_culture

    def set_Kc_r_I_Qs_settlement(self,value):
        I_settlement = 0.001

        self.Kc = 0.3
        self.r = 0.99
        self.I = I_settlement