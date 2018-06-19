## SCC0251 1st Semester 2018
## Final Project - BPCS Steganography
##   Guilherme dos Santos Marcon     9293564     ICMC-USP

import numpy as np
import imageio
import os
import time

##Global var
## all blocks and planes are 8x8
dAlpha = 0.3
## 8x8 matrix with [0, 0] = 0 and alternating
Wc = np.zeros((8, 8), dtype=np.uint8)
Wc[1::2, 0::2] = 1
Wc[0::2, 1::2] = 1

## max changes of values with neighbors, e.g.:
## [0 1 0; 1 0 1; 0 1 0]- (0, 0) has 2 variations, (0, 1) has 3 and (1, 1) has 4
maxChanges = 2*4 + 3*(8-2)*4 + np.power(8-2, 2)*4
## maxChanges = max changes  in corners + max changes in non-corner borders
##                  + max changes in the middle

bitWiseArray = np.array([128, 64, 32 ,16, 8, 4, 2, 1], dtype=np.uint8)
bitWiseOrder = np.array([7,   6,  5,  4,  3, 2, 1, 0], dtype=np.uint8)

## using classes to make the code more understandable
## VesselIterator is a class to iterate through and modify the vessel image
## It contais methods like: nextComplexPlane(), getPlane() and insertPlane(plane)
class VesselIterator:
    def __init__(self, vesselImage):
        ## Preparing the vessel, separating the image bit plane wise
        
        self.vessel = np.empty((vesselImage.shape[0],
                                  vesselImage.shape[1],
                                  vesselImage.shape[2],
                                  8), dtype=np.uint8)
        
        for i in range(0, 8):
            np.bitwise_and(vesselImage, bitWiseArray[i], out=self.vessel[:,:,:,i])
            self.vessel[:,:,:,i] = self.vessel[:,:,:,i] >> bitWiseOrder[i]

        self.finished = False
        self.__x = 0
        self.__y = 0
        self.__z = 0
        self.__B = 0

        if not isComplex(self.getPlane()):
            self.nextComplexPlane()

    def getVessel(self, oldVessel=None):
        if oldVessel is None:
            oldVessel = np.empty([self.vessel.shape[0],
                                  self.vessel.shape[1],
                                  self.vessel.shape[2]])
        oldVessel = np.multiply(self.vessel[:,:,:,0], bitWiseArray[0])
        for i in range(1, 8):
            oldVessel += np.multiply(self.vessel[:,:,:,i], bitWiseArray[i])
        return oldVessel

    ## Iterate to the next vessel plane
    def nextPlane(self):
        self.__x += 8
        if self.__x+8 > self.vessel.shape[0]:
            self.__x = 0
            self.__y += 8
            if self.__y+8 > self.vessel.shape[1]:
                self.__y = 0
                self.__z += 1
                if self.__z >= self.vessel.shape[2]:
                    self.__z = 0
                    self.__B += 1
                    if self.__B >= 8:
                        self.__B = 0
                        self.finished = True

    ## Iterate through vessel image to find the next complex plane
    def nextComplexPlane(self):
        while not self.finished:
            self.nextPlane()
            if isComplex(self.getPlane()): break

    ## Returns the 8x8 current vessel block
    def getBlock(self):
        return self.vessel[self.__x:self.__x+8, self.__y:self.__y+8, self.__z,:]

    ## Returns the 8x8 current vessel plane
    def getPlane(self):
        return self.vessel[self.__x:self.__x+8, self.__y:self.__y+8, self.__z, bitWiseOrder[self.__B]]

    ## Returns a string of the current coordinate, debug purpose
    def getCoord(self):
        return "("+str(self.__x)+", "+str(self.__y)+", "+str(self.__z)+", "+str(bitWiseOrder[self.__B])+")"

    ## Insert the 1 bit target plane at the Bth bit plane of the vessel
    ## Returning the boolean if the target plane was conjugated or not
    def insertPlane(self, plane):
        if not isComplex(plane):
            plane = conjugate(plane)
            wasConjugated = 1
        else:
            wasConjugated = 0
        self.insertBits(plane)
        return wasConjugated

    ## Properly insert the 1 bit target plane at the Bth bit plane of vessel
    def insertBits(self, plane):
        self.vessel[self.__x:self.__x+8, self.__y:self.__y+8, self.__z,bitWiseOrder[self.__B]] = plane

## TargetIterator is a class to iterate through the target image
## It contains methods like: nextPlane()
class TargetIterator:
    def __init__(self, targetImage):
        self.__target = targetImage
        self.__iter = np.nditer(self.__target, flags=['c_index'], op_flags=['readwrite'])
        self.finished = False

    ## Return a 8x8 binary bit plane of target's 'columns' next numbers
    ## if columns < 8, the first '8-columns' columns will be random values
    def nextPlane(self, columns=8):
        plane = (np.random.rand(8, 8)*2).astype(np.uint8)

        for i in range(8-columns, 8):
            bin_number = np.binary_repr(self.__iter[0], width=8)
            ## plane[:,i] = ord(bin_number[:])-ord('0') doesn't work
            for B in range(0, 8):
                plane[B, i] = ord(bin_number[B])-ord('0')
            self.next()
            if self.finished: break

        return plane

    ## Unhide methods
    def set(self, pixel):
        self.__iter[0] = pixel

    def next(self):
        self.__iter.iternext()
        if self.__iter.finished:
            self.finished = True

## ConjugationMap is the class that marks if determinated plane is conjugated or not
## It contains methods like: set(bit) and next()
class ConjugationMap:
    def __init__(self, targetShape):
        self.maxTargetBlocks = int(np.ceil(targetShape[0]*
                                       targetShape[1]*
                                       targetShape[2]/8))
        ## each bit of the conjugation map represents to a target image block
        self.maxConjBlocks = int(np.ceil(self.maxTargetBlocks/512))
        self.conjMap = np.zeros((self.maxConjBlocks, 8, 8), dtype=np.uint8)
        self.reset()

    def reset(self):
        self.__iter = np.nditer(self.conjMap, flags=['c_index'], op_flags=['readwrite'])
        self.finished = False
        self.__count = 0
        self.__B = 0

    def set(self, bit):
        self.__iter[0] += (bit*(np.power(2, self.__B))).astype(np.uint8)

    def get(self):
        return (np.bitwise_and(self.__iter[0], np.power(2, self.__B)) >> self.__B)

    def next(self):
        self.__iter.iternext()
        self.__count += 1
        if self.__iter.finished:
            self.__iter.reset()
            self.__B += 1
            if self.__B >= 8:
                self.__B = 0
                self.finished = True
        if self.__count >= self.maxTargetBlocks:
            self.finished = True        

## Hides T inside of V, by the BPCS method
## V: vessel image
## T: target image
def BPCS_hide(V, T):
    ## transforming from Pure Binary Code to Canonical Gray Code
    ## this make the insertion of the BPCS planes less intrusive
    print("Transforming vessel to CGC")
    V = PBCtoCGC(V)

    print("Preparing vessel image")
    Viter = VesselIterator(V)

    print("Inserting target shape")
    ## creating and inserting the initPlanes
    initPlanes = createInitPlanes(T)
    for i in range(0, 2):
        if Viter.finished: return None
        Viter.insertBits(initPlanes[i,:,:])
##        finitPlanes.write("I: "+str(i)+" - "+Viter.getCoord()+"\n")
        Viter.nextComplexPlane()

    if Viter.finished: return None

    print("Inserting image")
    ## creating and inserting the target image planes
    Titer = TargetIterator(T)
    CMiter = ConjugationMap(T.shape)
    while not CMiter.finished and not Viter.finished:
        CMiter.set(Viter.insertPlane(Titer.nextPlane()))
        
        CMiter.next()
        Viter.nextComplexPlane()

    if Viter.finished: return None
    
    print("Inserting conjugation map")
    CMnewIter = TargetIterator(CMiter.conjMap)
    while not CMnewIter.finished and not Viter.finished:
        plane = CMnewIter.nextPlane(columns=7)
        plane[0,0] = 0 ## small complexity change error
        if not isComplex(plane):
            plane = conjugate(plane)
            plane[0,0] = 1
        else:
            plane[0,0] = 0
        Viter.insertBits(plane)
        Viter.nextComplexPlane()

    if Viter.finished and not CMnewIter.finished: return None

    ## transforming back to Pure Binary Code
    print("Transforming vessel back to PBC")
    return CGCtoPBC(Viter.getVessel(oldVessel=V))


class HiddenIterator:
    def __init__(self, Tshape):
        self.__Tshape = Tshape
        self.maxPlanes = int(np.ceil(Tshape[0]*
                                     Tshape[1]*
                                     Tshape[2]/8))
        self.__planes = np.zeros((self.maxPlanes, 8, 8), dtype=np.uint8)
        self.reset()

    def reset(self):
        self.__i = 0
        self.finished = False

    def nextPlane(self):
        self.__i += 1
        if self.__i >= self.maxPlanes:
            self.__i = 0
            self.finished = True

    def getPlane(self):
        return self.__planes[self.__i, :, :]

    def setPlane(self, plane):
        self.__planes[self.__i, :, :] = plane

    def getImage(self, Ishape, columns=8):
        if not self.finished: return None
        
        T = np.zeros(Ishape, dtype=np.uint8)
        Titer = TargetIterator(T)

        while not Titer.finished:
            for i in range(8-columns, 8):
##                print("Setting: ", np.multiply(self.__planes[self.__i,:,i], auxVec).sum())
                Titer.set(np.multiply(self.__planes[self.__i,:,i], bitWiseArray).sum())
                Titer.next()
                if Titer.finished: break
            self.nextPlane()
        
        return T

def BPCS_unhide(V):
    print("Transforming vessel to CGC")
    V = PBCtoCGC(V)

    print("Preparing vessel image")
    Viter = VesselIterator(V)

    print("Recovering target shape")
    Tshape = recoverTargetShape(Viter)

    print("Checking if it's possible for vessel to contain hidden image")
    Tmult = 1
    for i in range(0, len(Tshape)):
        if Tshape[i] <= 0: return None
        Tmult *= Tshape[i]
    Vmult = 1
    for i in range(0, len(V.shape)):
        Vmult *= V.shape[i]

    if Tmult/Vmult >= 0.6: return None

    print("Recovering image")
    Hiter = HiddenIterator(Tshape)
    while not Hiter.finished and not Viter.finished:
        Hiter.setPlane(Viter.getPlane())
        Hiter.nextPlane()
        Viter.nextComplexPlane()

    if Viter.finished: return None

    print("Recovering conjugation map")
    maxTargetBlocks = int(np.ceil(Tshape[0]*Tshape[1]*Tshape[2]/8))
    maxConjBlocks = int(np.ceil(maxTargetBlocks/512))
    ## 1 original conj block makes (9+1/7) vessel planes substitution
    CMiter = HiddenIterator((int(np.ceil(maxConjBlocks*(9+1/7))), 8, 1))
    while not CMiter.finished and not Viter.finished:
        plane = Viter.getPlane()
        if plane[0,0] == 1:
            plane = conjugate(plane)
        CMiter.setPlane(plane)
        CMiter.nextPlane()
        Viter.nextComplexPlane()

    if not CMiter.finished and Viter.finished: return None

    print("Conjugating back the image")
    CMnewIter = ConjugationMap(Tshape)
    CMnewIter.conjMap = CMiter.getImage((maxConjBlocks, 8, 8), columns=7)
    CMnewIter.reset()
    Hiter.reset()
    while not Hiter.finished:
        if CMnewIter.get() == 1:
            Hiter.setPlane(conjugate(Hiter.getPlane()))
        CMnewIter.next()
        Hiter.nextPlane()
    
    print("Transforming back image")
    return Hiter.getImage(Tshape)
    

## Creates the initPlanes, which are 2 planes that contains T.shape
## The first half of the first plane contains which block is conjugated
def createInitPlanes(T):
    baseVec = (np.random.rand(16, 1, 1)*256).astype(np.uint8)

    for i in range(0, 3):       ## each T.shape[i]
        for j in range(0, 4):   ## each int8 of int32
            ## separating the 3*32 bit int into 3*4*8 bit int
            baseVec[4*(i+1)+j, 0, 0] = (np.bitwise_and(T.shape[i], 255 << j*8) >> j*8)

    initPlanes = np.zeros((2, 8, 8), dtype=np.uint8)
    baseVecIter = TargetIterator(baseVec)

    for i in range(0, 2):
        initPlanes[i,:,:] = baseVecIter.nextPlane()

        if not isComplex(initPlanes[i,:,:]):
            initPlanes[i,:,:] = conjugate(initPlanes[i,:,:])
            initPlanes[0,0,i] = 1
        else:
            initPlanes[0,0,i] = 0

    return initPlanes

## Unhide part of the create init planes part, will get the init planes and
## return the shape of target image
def recoverTargetShape(Viter):
    initPlanes = np.zeros((2, 8, 8), dtype=np.uint8)
    for i in range(0, 2):
        if Viter.finished: return None
        initPlanes[i,:,:] = Viter.getPlane()
        Viter.nextComplexPlane()
    
    if(initPlanes[0,0,1] == 1): initPlanes[1,:,:] = conjugate(initPlanes[1,:,:])
    if(initPlanes[0,0,0] == 1): initPlanes[0,:,:] = conjugate(initPlanes[0,:,:])

    shape = np.zeros((3), dtype=int)
    auxVec = np.array([128, 64, 32 ,16, 8, 4, 2, 1])
    for i in range(0, 3):
        fColumn = 4 if i%2 == 0 else 0
        for j in range(0, 4):
            shape[i] += (np.multiply(initPlanes[int(np.ceil(i/2)),:,fColumn+j], auxVec).sum() << j*8)

    return shape
    
## Returns True or False, if the binary plane is complex or not
def isComplex(binPlane):
    binPlane = binPlane.astype(np.int8) ## int 8
    L = binPlane[:, 0:-1]   ## all columns except the last
    R = binPlane[:, 1:]     ## all columns except the first
    T = binPlane[0:-1, :]   ## all lines except the last
    B = binPlane[1:, :]     ## all lines except the first

    count = 0.0
    count += np.absolute(L-R).sum()
    count += np.absolute(R-L).sum()
    count += np.absolute(T-B).sum()
    count += np.absolute(B-T).sum()

    return count/maxChanges > dAlpha

## Conjugates the binary plane, means bitwise_xor with Wc
def conjugate(binPlane):
    return np.bitwise_xor(binPlane, Wc)

## Converts the array/matrix from Pure Binary Code to Canonical Gray Code
def PBCtoCGC(PBC):
    CGC = np.zeros(PBC.shape, dtype=PBC.dtype)
    ## first column is equal
    CGC[:,0,:] = PBC[:,0,:]
    ## the others are equal PBC(column i-1) xor PBC(column i)
    for i in range(1, PBC.shape[1]):
        CGC[:,i,:] = np.bitwise_xor(PBC[:,i-1,:], PBC[:,i,:])
    return CGC

## Converts the array/matrix from Canonical Gray Code to Pure Binary Code
def CGCtoPBC(CGC):
    PBC = np.zeros(CGC.shape, dtype=CGC.dtype)
    ## first column is equal
    PBC[:,0,:] = CGC[:,0,:]
    ## the others are equal CGC(column i) xor PBC(column i-1)
    for i in range(1, CGC.shape[1]):
        PBC[:,i,:] = np.bitwise_xor(CGC[:,i,:], PBC[:,i-1,:])
    return PBC

def compareRMSE(H, HL):
    # HL = ^H
    return np.sqrt((np.power(H-HL, 2).sum())/(H.shape[0]*H.shape[1]))
    
## Main
print("BPCS Steganography\nInsert operation:\nHide -> 1\n"+
      "Unhide -> 2\nCompare RMSE -> 3\nDEMO -> 4\n")
operation = int(input())

if operation == 1:
    print("Vessel name:")
    vessel_name = str(input()).rstrip()
    V = imageio.imread(vessel_name).astype(np.uint8)
    Vmult = V.shape[0]*V.shape[1]*V.shape[2]
    
    print("Target name: ")
    target_name = str(input()).rstrip()
    T = imageio.imread(target_name).astype(np.uint8)
    Tmult = T.shape[0]*T.shape[1]*T.shape[2]
    
    print("Final image name:")
    final_name = str(input()).rstrip()

    t = time.time()
    F = BPCS_hide(V, T)
    t = time.time()-t

    f = open(final_name.split('.')[0]+".txt", "w+")
    f.write("Vessel: "+vessel_name+"\n\tSize: "+str(os.path.getsize(vessel_name))+" bytes\n")
    f.write("Target: "+target_name+"\n\tSize: "+str(os.path.getsize(target_name))+" bytes\n")
    f.write("Percentage:\t"+str.format("%.5f" % (100*Tmult/Vmult))+"\n")
    
    if F is None :
        f.write("Insufficient complex blocks to image be inserted.")
        print("Insufficient complex blocks to image be inserted.")
    else:
        imageio.imwrite(final_name, F)
        f.write("Time:\t\t"+str.format("%.5f" % t)+" seconds\n")
        f.write("RMSE:\t\t"+str.format("%.5f" % compareRMSE(V, F))+"\n")

    f.close()
    
elif operation == 2:
    print("Vessel name:")
    vessel_name = str(input()).rstrip()
    V = imageio.imread(vessel_name).astype(np.uint8)
    
    print("Final image name:")
    final_name = str(input()).rstrip()
    
    t = time.time()
    F = BPCS_unhide(V)
    t = time.time()-t
    
    if F is None: print("Error... Vessel image contains no message.")
    else:
        f = open(final_name.split('.')[0]+".txt", "w+")
        f.write("Time: "+str.format("%.5f" % t)+" seconds\n")
        imageio.imwrite(final_name, F)
        f.close()

elif operation == 3:
    print("First image:")
    first_name = str(input()).rstrip()
    print("Second image:")
    second_name = str(input()).rstrip()

    print("RMSE: ", compareRMSE(imageio.imread(first_name).astype(np.uint8),
                                imageio.imread(second_name).astype(np.uint8)))

elif operation == 4:
    print("Executing DEMO")
    print("Hiding image")
    print("Vessel image: images/1-mountain.png")
    print("Target image: images/1-colorful-smoke.png")
    V = imageio.imread("images/1-mountain.png").astype(np.uint8)
    T = imageio.imread("images/1-colorful-smoke.png").astype(np.uint8)

    H = BPCS_hide(V, T)
    
    print("Saving vessel with hidden image: result_hidden/DEMO.png")
    imageio.imwrite("result_hidden/DEMO.png", H)

    print()
    print("Unhiding image")
    U = BPCS_unhide(H)

    print("Saving unhidden image")
    imageio.imwrite("result_unhidden/DEMO.png", U)

    print("\nDEMO finished")
