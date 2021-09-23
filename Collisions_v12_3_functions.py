import math, random, numpy
from random import randrange

#Parameters
red = (255,0,0)

#Base atom definition
class atom:
    def __init__(self, rad, color, maxVel, boundary):
        self.ID = -1
        self.rad = rad
        self.mass = self.rad**2*math.pi
        self.color = color
        self.vel = [random.uniform(-maxVel, maxVel), random.uniform(-maxVel, maxVel)] 
        self.x = randrange(boundary[0] + 20 + self.rad + math.ceil(maxVel), boundary[1] - 20 - self.rad - math.ceil(maxVel))
        self.y = randrange(boundary[2] + 20 + self.rad + math.ceil(maxVel), boundary[3] - 20 - self.rad - math.ceil(maxVel))


#Takes a list of atoms and checks if a new atom overlaps anything, appending if no overlaps
def addAtom(atoms, sections, newAtom):
    for i in range(len(atoms)):
        dist = ((newAtom.x-atoms[i].x)**2 + (newAtom.y-atoms[i].y)**2)**0.5
        if dist <= newAtom.rad + atoms[i].rad:
            return False

    colorAtoms(newAtom,red)
    newAtom.ID = len(atoms) + 1
    atoms.append(newAtom)

    return True    


#Color atom based on starting position, left half of the screen starts as red
def colorAtoms(atom, color):
    if atom.x < 500:
        atom.color = color


#Calculates the angle between two 2D vectors, A and B, in radians
def vectorAngle(A, B):
    alphaRad = numpy.arccos((A[0]*B[0]+A[1]*B[1])/((A[0]**2+A[1]**2)**0.5*(B[0]**2+B[1]**2)**0.5))
    return alphaRad


#Calculates the magnitude of a 2D vector
def vectorMag(v):
    magnitude = ((v[0]**2+v[1]**2)**0.5)
    return magnitude


#Elastic collisions with other atoms, dependent on mass. Gets two atoms, returns two atoms with updated parameters
def atomCollision(a1, a2):
    #Velocity vectors
    VA1 = [a1.vel[0], a1.vel[1]]
    VA2 = [a2.vel[0], a2.vel[1]]
    velDiff = [VA1[0]-VA2[0], VA1[1]-VA2[1]]
    velDiffMag = vectorMag(velDiff)

    #Calculate time t to backtrack to collision point
    distVec = [a1.x-a2.x,a1.y-a2.y]
    thetaO = vectorAngle(distVec,[velDiff[0],velDiff[1]])
    phiO = math.asin(vectorMag(distVec)*math.sin(thetaO)/(a1.rad+a2.rad))
    chiO = math.pi-thetaO-phiO
    timeOverlap = math.sin(chiO)/(math.sin(thetaO)/(a1.rad+a2.rad)*velDiffMag)
    
    #Subtract original velocities * timeOverlap
    a1.x -= timeOverlap * a1.vel[0]
    a1.y -= timeOverlap * a1.vel[1]
    a2.x -= timeOverlap * a2.vel[0]
    a2.y -= timeOverlap * a2.vel[1]

    #Calculate collision mechanics at the actual point of collision
    #Line between two atom centerpoints
    colVector = [a1.x-a2.x,a1.y-a2.y]
    colDist = ((colVector[0]**2 + colVector[1]**2)**0.5)
            
    #Angle between collision vector and velocity difference
    theta = vectorAngle(velDiff, colVector)
    tMag = velDiffMag*math.cos(theta)
    velTransfer = [colVector[0]*tMag/colDist,colVector[1]*tMag/colDist]
    
    #Coefficients for adjusting velocity transfer based on mass, 1D collision with the second object stationary relative to velocity transfer
    m1Ratio2 = 2*a2.mass/(a1.mass+a2.mass)
    m2Ratio1 = 2*a1.mass/(a1.mass+a2.mass)

    #Calculate new velocity vector
    a1.vel[0] -= m1Ratio2*velTransfer[0]
    a1.vel[1] -= m1Ratio2*velTransfer[1]
    a2.vel[0] += m2Ratio1*velTransfer[0]
    a2.vel[1] += m2Ratio1*velTransfer[1]
    
    #Add new velocities * timeOverlap
    a1.x += timeOverlap * a1.vel[0]
    a1.y += timeOverlap * a1.vel[1]
    a2.x += timeOverlap * a2.vel[0]
    a2.y += timeOverlap * a2.vel[1]


#Calculate elastic collisions off of the boundary walls, given an atom object and a rectangular boundary[xMin,xMax,yMin,yMax]
def boundaryCollision(atom, boundary):
    #Left wall
    if atom.x < boundary[0] + atom.rad :
        atom.vel[0] = -atom.vel[0]
        #Adjust posiiton for how much the atom should have bounced off the wall
        atom.x += 2*(boundary[0] + atom.rad - atom.x)
    #Right wall
    elif atom.x > boundary[1] - atom.rad:
        atom.vel[0] = -atom.vel[0]
        atom.x += 2*(boundary[1] - atom.rad - atom.x)
    #Top wall
    if atom.y < boundary[2] + atom.rad:
        atom.vel[1] = -atom.vel[1]
        atom.y += 2*(boundary[0] + atom.rad - atom.y)
    #Bottom wall
    elif atom.y > boundary[3] - atom.rad:
        atom.vel[1] = -atom.vel[1]
        atom.y += 2*(boundary[1] - atom.rad - atom.y)


#Creates a 2D array of size numX by numY
def sectionsBuild(numX,numY):
    sects = []
    for i in range(int(numX)):
        sects.append([])
        for j in range(int(numY)):
            sects[i].append([])

    return sects


#Defines location coefficients based on section definitions
def locCoefficients(atom, boundary, baseDim, overlapDim):
    coeffs = [[],[]]
    #Base atom locations
    coeffs[0].append((atom.x - boundary[0])/baseDim)
    coeffs[1].append((atom.y - boundary[2])/baseDim)
    #Positions to check to see if atoms are in neighboring overlaps
    coeffs[0].append(((atom.x-boundary[0]) - overlapDim)/baseDim)
    coeffs[0].append(((atom.x-boundary[0]) + overlapDim)/baseDim)
    coeffs[1].append(((atom.y-boundary[2]) - overlapDim)/baseDim)
    coeffs[1].append(((atom.y-boundary[2]) + overlapDim)/baseDim)

    return coeffs


#Allocates an atom to appropriate collision sections given coefficients based on section definitions
def sectionsAllocate(atom,sections,locX1,locX2,locX3,locY1,locY2,locY3):
    #Add to base section
    sections[int(locX1)][int(locY1)].append(atom)       
    #Add to left 3 overlapping sections if located within them
    if int(locX2) != int(locX1):
        sections[int(locX2)][int(locY1)].append(atom)
        if int(locY2) != int(locY1):
            sections[int(locX2)][int(locY2)].append(atom)
        elif int(locY3) != int(locY1):
            sections[int(locX2)][int(locY3)].append(atom)

    #Add to right 3 overlapping sections if located within them
    elif int(locX3) != int(locX1):
        sections[int(locX3)][int(locY1)].append(atom)
        if int(locY2) != int(locY1):
            sections[int(locX3)][int(locY2)].append(atom)
        elif int(locY3) != int(locY1):
            sections[int(locX3)][int(locY3)].append(atom)

    #Add to top and bottom middle overlap sections if located within them
    if int(locY2) != int(locY1):
        sections[int(locX1)][int(locY2)].append(atom)
    elif int(locY3) != int(locY1):
        sections[int(locX1)][int(locY3)].append(atom)
  















