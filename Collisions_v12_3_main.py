import pygame as pyg
import pygame.freetype, decimal, math, sys
import Collisions_v12_3_functions as funcs
from random import randrange

pygame.init()

#Base simulation parameters
xDim = 1300
yDim = 1000
win = pygame.display.set_mode((xDim, yDim))
pygame.display.set_caption("Collisions")

timeStep = 10
red = (255,0,0)
green = (0,255,0)
blue = (0,0,255)
black = (0,0,0)
boundary = [100,900,100,900]
boundaryWidth = 1
maxVel = 1
maxRad = 30
atomCounter = 0
areaCapacity = 0.35
capacity = 0

#Collision section parameters
baseSectionDim = maxRad*5
overlapDim = maxRad
numXSections = ((boundary[1]-boundary[0])/baseSectionDim) + 1
numYSections = ((boundary[3]-boundary[2])/baseSectionDim) + 1

#Slider control variables
pastButtonState = 0
baseSliderHeight = 110
clickX = 0
clickY = 0
objSliderXMin = 955
objSliderXMax = 1045
objSliderMult = 1
radSliderXMin = 1055
radSliderXMax = 1145
radSliderMult = 0.2
velSliderXMin = 1155
velSliderXMax = 1245
velSliderMult = 0.01

font = pygame.freetype.SysFont(pygame.font.get_default_font(), 18)

atoms = []

#Manually set initial conditions for testing purposes
"""
atoms.append(funcs.atom(50, red, 0, boundary))
atoms[0].ID = 0
atoms[0].mass = atoms[0].rad**2*math.pi
atoms[0].vel = [-2,0]
atoms[0].x = 650
atoms[0].y = 450
atomCounter += 1

atoms.append(funcs.atom(50, blue, 0, boundary))
atoms[1].ID = 1
atoms[1].mass = atoms[1].rad**2*math.pi
atoms[0].vel = [2,0]
atoms[1].x = 350
atoms[1].y = 530
atomCounter += 1
"""

#Simulation loop
run = True
while run:
    #Creates time delay of timeStep in ms
    pygame.time.delay(timeStep)
    
    #Check for ESC key pressed or pygame window closed, to quit
    for event in pygame.event.get():
        if event.type == pyg.QUIT:
            pygame.quit()
            sys.exit()
        elif event.type == pyg.KEYDOWN:
            if event.key == pyg.K_ESCAPE:
                pygame.quit()
                sys.exit()

    #Refresh screen and draw boundary
    win.fill((255, 255, 255))
    pygame.draw.rect(win, red, pygame.Rect(100,100,800,800),boundaryWidth)


    #Object number slider
    if clickX >= objSliderXMin and clickX <= objSliderXMax:
        #Removing atoms
        if clickY*objSliderMult < atomCounter:
            for i in range(atomCounter - clickY*objSliderMult):
                if len(atoms) > 0:
                    atoms.pop(0)
                    atomCounter -= 1
        #Adding atoms if they get placed without overlapping anything
        elif clickY*objSliderMult > atomCounter:
            for i in range(clickY*objSliderMult - atomCounter):
                #Don't add any more atoms if the screen is too full
                areaFilled = 0
                areaTotal = (boundary[1]-boundary[0])*(boundary[3]-boundary[2])
                for j in range(len(atoms)):
                    areaFilled += math.pi*atoms[j].rad**2
                    capacity = areaFilled/areaTotal
                if areaFilled/areaTotal >= areaCapacity:
                    break
                tempAtom = funcs.atom(maxRad, blue, maxVel, boundary)
                #See if atom overlaps anything before adding it
                if funcs.addAtom(atoms, sections, tempAtom):
                    atomCounter += 1
                
        #Update collision section parameters to work for the largest object on the screen
        radList = [1]
        for i in range(len(atoms)):
            radList.append(atoms[i].rad)
        biggestRad = max(radList)
        baseSectionDim = math.ceil(biggestRad/5)*25
        overlapDim = biggestRad
        numXSections = ((boundary[1]-boundary[0])/baseSectionDim) + 1
        numYSections = ((boundary[3]-boundary[2])/baseSectionDim) + 1

                  
    #New object radius slider
    if clickX >= radSliderXMin and clickX <= radSliderXMax:
        if clickY*radSliderMult != maxRad and clickY*radSliderMult > 0:
            maxRad = math.ceil(clickY*radSliderMult)


    #New object velocity slider
    if clickX >= velSliderXMin and clickX <= velSliderXMax:
        if clickY*velSliderMult != maxVel and clickY*velSliderMult > 0:
            maxVel = math.floor(clickY*velSliderMult)


    #Update atoms position
    for i in range(len(atoms)):
        atoms[i].x += atoms[i].vel[0]
        atoms[i].y += atoms[i].vel[1]

        
    #Make an empty 2D array of all collision sections
    sections = funcs.sectionsBuild(numXSections, numYSections)

    #Calculate boundary collisions and allocate atoms to sections based on position to calculate atom-atom collisions
    for i in range(len(atoms)):
        
        #Boundary collisions
        funcs.boundaryCollision(atoms[i],boundary)

        #Location coefficient calculation
        locCoeffs = funcs.locCoefficients(atoms[i], boundary, baseSectionDim, overlapDim)

        #Allocate atoms to appropriate sections based on location
        funcs.sectionsAllocate(atoms[i],sections,locCoeffs[0][0],locCoeffs[0][1],locCoeffs[0][2],locCoeffs[1][0],locCoeffs[1][1],locCoeffs[1][2])

 
    #Check for collisions within subsections with direct check between all atoms in subsection
    for section in sections:
        for subSection in section:
            if len(subSection) > 1:
                for a in range(len(subSection)):
                    for b in range(a+1,len(subSection)):
                        dist = ((subSection[a].x-subSection[b].x)**2 + (subSection[a].y-subSection[b].y)**2)**0.5
                        if dist <= subSection[a].rad + subSection[b].rad:
                            funcs.atomCollision(subSection[a],subSection[b])
                                  

    #Control slider interaction. Get x coord only when first clicking, get y coord continuously
    button = pygame.mouse.get_pressed()
    if button[0]:
        pos = pygame.mouse.get_pos()
        clickY = yDim - baseSliderHeight- pos[1]
        if pastButtonState == 0:
            pastButtonState = 1
            clickX = pos[0]
    #Reset clicking variables when releasing button
    elif not button[0] and pastButtonState == 1:
        pastButtonState = 0
        clickX = 0
        clickY = 0


    #Draw sliders
    pygame.draw.rect(win,black,pyg.Rect(objSliderXMin-5,yDim-baseSliderHeight-atomCounter/objSliderMult,objSliderXMax-objSliderXMin,10))
    pygame.draw.rect(win,black,pyg.Rect(radSliderXMin,yDim-baseSliderHeight-(maxRad-1)/radSliderMult,radSliderXMax-radSliderXMin,10))
    pygame.draw.rect(win,black,pyg.Rect(velSliderXMin,yDim-baseSliderHeight-(maxVel)/velSliderMult,velSliderXMax-velSliderXMin,10))

    #Write slider labels
    font.render_to(win, (objSliderXMin, yDim-50), "Atoms: " + str(atomCounter), black)
    font.render_to(win, (radSliderXMin+15, yDim-50), "Radii: " + str(maxRad), black)
    font.render_to(win, (velSliderXMin, yDim-50), "Velocity: " + str(maxVel), black)

    
    #Draw all the atoms
    for i in range(len(atoms)):
        pygame.draw.circle(win, atoms[i].color, (atoms[i].x, atoms[i].y), atoms[i].rad)

    #Refreshes the window
    pygame.display.update()

# closes the pygame window
pygame.quit()


