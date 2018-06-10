## SCC0251 1st Semester 2018
## Final Project - BPCS Steganography
##   Guilherme dos Santos Marcon     9293564     ICMC-USP

import numpy as np
import imageio

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

## using classes to make the code more understandable
## VesselIterator is a class to iterate through the vessel image
## It contais methods like: nextComplexPlane(), getPlane() and insertPlane(plane)
class VesselIterator:
    def __init__(self, vesselImage):
        self.__vessel = vesselImage
        self.finished = False
        self.__x = 0
        self.__y = 0
        self.__z = 0
        self.__B = 0

        if not isComplex(self.getPlane()):
            self.nextComplexPlane()

    ## Iterate to the next vessel plane
    def nextPlane(self):
        self.__x += 8
        if self.__x+8 > self.__vessel.shape[0]:
            self.__x = 0
            self.__y += 8
            if self.__y+8 > self.__vessel.shape[1]:
                self.__y = 0
                self.__z += 1
                if self.__z >= self.__vessel.shape[2]:
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
        return self.__vessel[self.__x:self.__x+8, self.__y:self.__y+8, self.__z]

    ## Returns the 8x8 current vessel plane
    def getPlane(self):
        return (np.bitwise_and(self.getBlock(), np.power(2, self.__B)) >> self.__B)

    ## Returns a string of the current coordinate, debug purpose
    def getCoord(self):
        return "("+str(self.__x)+", "+str(self.__y)+", "+str(self.__z)+", "+str(self.__B)+")"

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
        self.__vessel[self.__x:self.__x+8,
                      self.__y:self.__y+8,
                      self.__z] = (np.bitwise_and(
                          self.getBlock(),
                          255-np.power(2, self.__B))) + (plane << self.__B)

## TargetIterator is a class to iterate through the target image
## It contains methods like: nextPlane()
class TargetIterator:
    def __init__(self, targetImage):
        self.__target = targetImage
        self.__iter = np.nditer(self.__target, flags=['c_index'], op_flags=['readonly'])
        self.finished = False

    ## Return the binary bit plane of target's 8 next numbers
    def nextPlane(self):
        return self.nextCustomPlane(8)

    ## Return a 8x8 binary bit plane of target's 'columns' next numbers
    ## if columns < 8, the first '8-columns' columns will be random values
    def nextCustomPlane(self, columns):
        plane = (np.random.rand(8, 8)*2).astype(np.uint8)

        for i in range(8-columns, 8):
            bin_number = np.binary_repr(self.__iter[0], width=8)
            for B in range(0, 8):
                plane[B, i] = ord(bin_number[B])-ord('0')
            self.__iter.iternext()
            if self.__iter.finished:
                self.finished = True
                break

        return plane

## ConjugationMap is the class that marks if determinated plane is conjugated or not
## It contains methods like: set(bit) and next()
class ConjugationMap:
    def __init__(self, targetImage):
        self.maxPlanes = int(np.ceil(targetImage.shape[0]*
                                       targetImage.shape[1]*
                                       targetImage.shape[2]/8))
        ## each bit of the conjugation map represents to a target image block
        self.maxConjBlocks = int(np.ceil(self.maxPlanes/512))
        self.conjMap = np.zeros((self.maxConjBlocks, 8, 8), dtype=np.uint8)
        self.__iter = np.nditer(self.conjMap, flags=['c_index'], op_flags=['readwrite'])
        self.finished = False
        self.__count = 0;
        self.__B = 0;

    def set(self, bit):
        self.__iter[0] += (bit*(np.power(2, self.__B))).astype(np.uint8)

    def next(self):
        self.__iter.iternext()
        self.__count += 1
        if self.__iter.finished:
            self.__iter.reset()
            self.__B += 1
            if self.__B >= 8:
                self.__B = 0
                self.finished = True
        if self.__count >= self.maxPlanes:
            self.finished = True        

## Hides T inside of V, by the BPCS method
## V: vessel image
## T: target image
def BPCS_hide(V, T):
    ## transforming from Pure Binary Code to Canonical Gray Code
    ## this make the insertion of the BPCS planes less intrusive
    V = PBCtoCGC(V)

    Viter = VesselIterator(V)

    finitPlanes = open("./initPlanes.txt", "w+")
    finitPlanes.write("T.shape: "+str(T.shape[0])+" - "+str(T.shape[1])+" - "+str(T.shape[2])+"\n")

    print("Inserting target shape")
    ## creating and inserting the initPlanes
    initPlanes = createInitPlanes(T)
    for i in range(0, 2):
        Viter.nextComplexPlane()
        if Viter.finished: return None
        Viter.insertBits(initPlanes[i,:,:])
        finitPlanes.write("I: "+str(i)+" - "+Viter.getCoord()+"\n")

    print("Inserting image")
    ## creating and inserting the target image planes
    Titer = TargetIterator(T)
    CMiter = ConjugationMap(T)
    while not CMiter.finished:
        Viter.nextComplexPlane()
        if Viter.finished: return None
        CMiter.set(Viter.insertPlane(Titer.nextPlane()))
        CMiter.next()

    print("Inserting conjugation map")
    CMnewIter = TargetIterator(CMiter.conjMap)
    while not CMnewIter.finished:
        Viter.nextComplexPlane()
        if Viter.finished: return None
        plane = CMnewIter.nextPlane()
        if not isComplex(plane):
            plane = conjugate(plane)
            plane[0,0] = 1
        else:
            plane[0,0] = 0
        Viter.insertBits(plane)

    print("Finished")
    
    finitPlanes.write("maxPlanes: "+str(CMiter.maxPlanes)+" maxConjBlocks: "+str(CMiter.maxConjBlocks)+"\n")
    finitPlanes.write(np.array2string(initPlanes))
    finitPlanes.close()
    
    ## transforming back to Pure Binary Code
    return CGCtoPBC(V)


class HiddenIterator:
    def __init__(self, Tshape):
        self.__Tshape = Tshape
        self.maxPlanes = int(np.ceil(Tshape[0]*
                                     Tshape[1]*
                                     Tshape[2]/8))
        self.__planes = np.zeros((maxPlanes, 8, 8), dtype=np.uint8)
        self.__i = 0
        self.finished = 0

    def nextPlane(self):
        self.__i += 1
        if self.__i >= self.maxPlanes:
            self.__i = 0
            self.finished = True

    def setPlane(self, plane):
        self.__planes[self.__i, :, :] = plane

    def getImage(self, conjMap):
        if not self.finished: return None
        T = np.zeros(self.__Tshape, dtype=np.uint8)
        
        return T


def BPCS_unhide(V):
    V = PBCtoCGC(V)

    Viter = VesselIterator(V)

    print("Recovering target shape")
    shape = recoverTargetShape(Viter)
    print(shape)

    print("Recovering image")
    Titer = HiddenIterator(shape)
    while not Titer.finished:
        Viter.nextComplexPlane()
        if Viter.finished: return None
        Titer.setPlane(Viter.getPlane())
        Titer.nextPlane()

    print("Recovering conjugation map")
    maxPlanes = int(np.ceil(targetImage.shape[0]*
                                       targetImage.shape[1]*
                                       targetImage.shape[2]/8))
    maxConjBlocks = int(np.ceil(self.maxPlanes/512))
    CMiter = HiddenIterator((maxConjBlocks, 8, 8))
    while not CMiter.finished:
        Viter.nextComplexPlane()
        if Viter.finished: return None
        CMiter.setPlane(Viter.getPlane())
        CMiter.nextPlane()
    
    return None

    return CGCtoPBC(V)
    

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
        Viter.nextComplexPlane()
        if Viter.finished: return None
        print(Viter.getPlane())
        initPlanes[i,:,:] = Viter.getPlane()
    
    if(initPlanes[0,0,1] == 1): initPlanes[1,:,:] = conjugate(initPlanes[1,:,:])
    if(initPlanes[0,0,0] == 1): initPlanes[0,:,:] = conjugate(initPlanes[0,:,:])

    print(initPlanes)

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
print("BPCS Steganography\nInsert operation:\nHide -> 1\nUnhide -> 2")
operation = int(input())

if operation == 1:
    print("Vessel name:")
    vessel_name = str(input()).rstrip()
    print("Target name: ")
    target_name = str(input()).rstrip()
    print("Final image name:")
    final_name = str(input()).rstrip()

    V = imageio.imread(vessel_name).astype(np.uint8)
    T = imageio.imread(target_name).astype(np.uint8)
    F = BPCS_hide(V, T)
    if F is None: print("Error...")

    imageio.imwrite(final_name, F)
    print(compareRMSE(V, F))
    
elif operation == 2:
    print("Vessel name:")
    vessel_name = str(input()).rstrip()
    print("Final image name:")
    final_name = str(input()).rstrip()

    V = imageio.imread(vessel_name).astype(np.uint8)
    F = BPCS_unhide(V)
    if F is None: print("Error...")

    imageio.imwrite(final_name, F)

elif operation == 3:
    V = imageio.imread('reduced-road.png').astype(np.uint8)
    print(V.shape)
    Viter = VesselIterator(V)
    Viter.B = 7
    while not Viter.finished:
        print("["+str(Viter.x)+", "+str(Viter.y)+", "+str(Viter.z)+", "+str(Viter.B)+"]")
        Viter.nextPlane()
    print(V.shape)
    
elif operation == 4:
    name = str(input()).rstrip()
    V = imageio.imread(name).astype(np.uint8)
    imageio.imwrite(name, V)
