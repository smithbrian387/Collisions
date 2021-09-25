"""This is the main module for an elastic collision physics simulation
of circular objects.
"""

import decimal
import math
import sys
from random import randrange

import pygame as pyg
import pygame.freetype

import collisions_functions as funcs


pygame.init()

# Base simulation parameters
XDIM = 1300
YDIM = 1000
win = pygame.display.set_mode((XDIM, YDIM))
pygame.display.set_caption("Collisions")

TIMESTEP = 10
red = (255, 0, 0)
green = (0, 255, 0)
blue = (0, 0, 255)
black = (0, 0, 0)
boundary = [100, 900, 100, 900]
BOUNDARY_WIDTH = 1
max_vel = 1
max_rad = 25
atom_counter = 0
AREA_CAPACITY = 0.20
capacity = 0

# Collision section parameters
base_section_dim = max_rad * 5
overlap_dim = max_rad
num_x_sections = int(((boundary[1] - boundary[0]) / base_section_dim) + 1)
num_y_sections = int(((boundary[3] - boundary[2]) / base_section_dim) + 1)

# Slider control variables
past_button_state = 0
BASE_SLIDER_Y = 110
click_x = 0
click_y = 0
OBJ_SLIDER_X_MIN = 955
OBJ_SLIDER_X_MAX = 1045
OBJ_SLIDER_MULT = 1
RAD_SLIDER_X_MIN = 1055
RAD_SLIDER_X_MAX = 1145
RAD_SLIDER_MULT = 0.2
VEL_SLIDER_X_MIN = 1155
VEL_SLIDER_X_MAX = 1245
VEL_SLIDER_MULT = 0.01

font = pygame.freetype.SysFont(pygame.font.get_default_font(), 18)

atoms = []

# Manually set initial conditions for testing purposes
"""
atoms.append(funcs.Atom(50, red, 0, boundary))
atoms[0].ID = 0
atoms[0].mass = atoms[0].rad**2 * math.pi
atoms[0].vel = [-2, 0]
atoms[0].x = 650
atoms[0].y = 450
atom_counter += 1

atoms.append(funcs.Atom(50, blue, 0, boundary))
atoms[1].ID = 1
atoms[1].mass = atoms[1].rad**2 * math.pi
atoms[0].vel = [2, 0]
atoms[1].x = 350
atoms[1].y = 530
atom_counter += 1
"""


# Simulation loop
run = True
while run:
    # Creates time delay of TIMESTEP in ms
    pygame.time.delay(TIMESTEP)
    
    # Check for ESC key pressed or pygame window closed, to quit
    for event in pygame.event.get():
        if event.type == pyg.QUIT:
            pygame.quit()
            sys.exit()
        elif event.type == pyg.KEYDOWN:
            if event.key == pyg.K_ESCAPE:
                pygame.quit()
                sys.exit()

    # Refresh screen and draw boundary
    win.fill((255, 255, 255))
    pygame.draw.rect(win, red, pygame.Rect(100, 100, 800, 800), BOUNDARY_WIDTH)


    # Object number slider
    if click_x >= OBJ_SLIDER_X_MIN and click_x <= OBJ_SLIDER_X_MAX:
        # Removing atoms
        if click_y*OBJ_SLIDER_MULT < atom_counter:
            for i in range(atom_counter - click_y * OBJ_SLIDER_MULT):
                if len(atoms) > 0:
                    atoms.pop(0)
                    atom_counter -= 1
        # Adding atoms if they get placed without overlapping anything
        elif click_y * OBJ_SLIDER_MULT > atom_counter:
            for i in range(click_y * OBJ_SLIDER_MULT - atom_counter):
                # Don't add any more atoms if the screen is too full
                area_filled = 0
                area_total = ((boundary[1] - boundary[0])
                             * (boundary[3] - boundary[2]))
                for j in range(len(atoms)):
                    area_filled += math.pi * atoms[j].rad**2
                    capacity = area_filled / area_total
                if area_filled / area_total >= AREA_CAPACITY:
                    break
                temp_atom = funcs.Atom(max_rad, blue, max_vel, boundary)
                # See if atom overlaps anything before adding it
                if funcs.add_atom(atoms, sections, temp_atom):
                    atom_counter += 1
                
        # Update section parameters to work for largest object in atoms
        rad_list = [1]
        for i in range(len(atoms)):
            rad_list.append(atoms[i].rad)
        biggest_rad = max(rad_list)
        base_section_dim = math.ceil(biggest_rad / 5) * 25
        overlap_dim = biggest_rad
        num_x_sections = int(((boundary[1] - boundary[0])
                              / base_section_dim) + 1)
        num_y_sections = int(((boundary[3] - boundary[2])
                              / base_section_dim) + 1)

                  
    # New object radius slider
    if click_x >= RAD_SLIDER_X_MIN and click_x <= RAD_SLIDER_X_MAX:
        if (click_y * RAD_SLIDER_MULT != max_rad and
                click_y * RAD_SLIDER_MULT > 0):
            max_rad = math.ceil(click_y * RAD_SLIDER_MULT)


    # New object velocity slider
    if click_x >= VEL_SLIDER_X_MIN and click_x <= VEL_SLIDER_X_MAX:
        if (click_y * VEL_SLIDER_MULT != max_vel and
                click_y * VEL_SLIDER_MULT > 0):
            max_vel = math.floor(click_y * VEL_SLIDER_MULT)


    # Update atoms position
    for i in range(len(atoms)):
        atoms[i].x += atoms[i].vel[0]
        atoms[i].y += atoms[i].vel[1]

        
    # Make an empty 2D array of all collision sections
    sections = funcs.sections_build(num_x_sections, num_y_sections)

    # Calculate boundary collisions and allocate atoms to sections
    # based on position to calculate atom-atom collisions
    if len(sections) > 1:
        for i in range(len(atoms)):
            # Boundary collisions
            funcs.boundary_collision(atoms[i], boundary)

            # Location coefficient calculation
            loc_coeffs = funcs.loc_coefficients(atoms[i], boundary,
                                              base_section_dim, overlap_dim)

            # Allocate atoms to appropriate sections based on location
            funcs.sections_allocate(atoms[i], sections,
                                   loc_coeffs[0][0], loc_coeffs[0][1],
                                   loc_coeffs[0][2], loc_coeffs[1][0],
                                   loc_coeffs[1][1], loc_coeffs[1][2])

 
    # Check for collisions within sub_sections
    # Does a direct check between all atoms in sub_section
    for section in sections:
        for sub_section in section:
            if len(sub_section) > 1:
                for a in range(len(sub_section)):
                    for b in range(a+1, len(sub_section)):
                        dist = ((sub_section[a].x - sub_section[b].x)**2
                                +(sub_section[a].y - sub_section[b].y)**2)**0.5
                        if dist <= sub_section[a].rad + sub_section[b].rad:
                            funcs.atom_collision(sub_section[a],
                                                 sub_section[b])
                                  

    # Control slider interaction
    # Get x coord only when first clicking, get y coord continuously
    button = pygame.mouse.get_pressed()
    if button[0]:
        pos = pygame.mouse.get_pos()
        click_y = YDIM - BASE_SLIDER_Y - pos[1]
        if past_button_state == 0:
            past_button_state = 1
            click_x = pos[0]
    # Reset clicking variables when releasing button
    elif not button[0] and past_button_state == 1:
        past_button_state = 0
        click_x = 0
        click_y = 0


    # Draw sliders
    pygame.draw.rect(win,black,pyg.Rect(OBJ_SLIDER_X_MIN - 5,
                                        YDIM - BASE_SLIDER_Y - 5 - atom_counter
                                        / OBJ_SLIDER_MULT,
                                        OBJ_SLIDER_X_MAX - OBJ_SLIDER_X_MIN,
                                        10))
    pygame.draw.rect(win,black,pyg.Rect(RAD_SLIDER_X_MIN,
                                        YDIM - BASE_SLIDER_Y - 5
                                        - (max_rad - 1) / RAD_SLIDER_MULT,
                                        RAD_SLIDER_X_MAX - RAD_SLIDER_X_MIN,
                                        10))
    pygame.draw.rect(win,black,pyg.Rect(VEL_SLIDER_X_MIN,
                                        YDIM - BASE_SLIDER_Y - 5 - (max_vel)
                                        / VEL_SLIDER_MULT,
                                        VEL_SLIDER_X_MAX - VEL_SLIDER_X_MIN,
                                        10))

    # Write slider labels
    font.render_to(win, (OBJ_SLIDER_X_MIN, YDIM - 50),
                   "Atoms: " + str(atom_counter), black)
    font.render_to(win, (RAD_SLIDER_X_MIN + 15, YDIM - 50),
                   "Radii: " + str(max_rad), black)
    font.render_to(win, (VEL_SLIDER_X_MIN + 5, YDIM - 50),
                   "Max Speed: " + str(max_vel), black)

    
    # Draw all the atoms
    for i in range(len(atoms)):
        pygame.draw.circle(win, atoms[i].color, (atoms[i].x, atoms[i].y),
                           atoms[i].rad)

    # Refreshes the window
    pygame.display.update()

pygame.quit()
