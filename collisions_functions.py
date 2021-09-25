"""This module defines atom class and functions used in an elastic
collision physics simulation.
"""

import math
import random
import numpy
from random import randrange, uniform


# Parameters
red = (255,0,0)


class Atom:
    """ Defines an atom structure. """
    def __init__(self, rad, color, max_vel, boundary):
        self.ID = -1
        self.rad = rad
        self.mass = self.rad**2 * math.pi
        self.color = color
        self.vel = [uniform(-max_vel, max_vel), uniform(-max_vel, max_vel)] 
        self.x = randrange(boundary[0] + self.rad + math.ceil(max_vel),
                           boundary[1] - self.rad - math.ceil(max_vel))
        self.y = randrange(boundary[2] + self.rad + math.ceil(max_vel),
                           boundary[3] - self.rad - math.ceil(max_vel))


def add_atom(atoms, sections, new_atom):
    """Takes a list of atoms and checks if a new atom overlaps anything.
    new_atom is appended to the atoms list if no overlaps.
    """
    for i in range(len(atoms)):
        dist = ((new_atom.x-atoms[i].x)**2 + (new_atom.y-atoms[i].y)**2)**0.5
        if dist <= new_atom.rad + atoms[i].rad:
            return False

    color_atoms(new_atom, red)
    new_atom.ID = len(atoms) + 1
    atoms.append(new_atom)

    return True


def color_atoms(atom, color):
    """Color atom based on starting position."""
    if atom.x < 500:
        atom.color = color


def vector_angle(a, b):
    """
    Calculates the angle between two 2D vectors, A and B, in radians.
    """

    alpha_rad = numpy.arccos((a[0]*b[0] + a[1]*b[1])
                            /((a[0]**2 + a[1]**2)**0.5
                              * (b[0]**2 + b[1]**2)**0.5))
    return alpha_rad


def vector_mag(v):
    """Calculates the magnitude of a 2D vector"""
    magnitude = ((v[0]**2 + v[1]**2)**0.5)
    return magnitude


def atom_collision(a1, a2):
    """Elastic collisions with other atoms, dependent on mass.
    Gets two atoms, returns two atoms with updated parameters
    """
    # Velocity vectors
    v_a1 = [a1.vel[0], a1.vel[1]]
    v_a2 = [a2.vel[0], a2.vel[1]]
    vel_diff = [v_a1[0] - v_a2[0], v_a1[1] - v_a2[1]]
    vel_diff_mag = vector_mag(vel_diff)

    # Calculate time t to backtrack to collision point
    dist_vec = [a1.x - a2.x,a1.y - a2.y]
    theta_overlap = vector_angle(dist_vec, [vel_diff[0], vel_diff[1]])
    phi_overlap = math.asin(vector_mag(dist_vec) * math.sin(theta_overlap)
                            / (a1.rad + a2.rad))
    chi_overlap = math.pi - theta_overlap - phi_overlap
    time_overlap = (math.sin(chi_overlap) /
                    (math.sin(theta_overlap)*vel_diff_mag / (a1.rad + a2.rad)))
    
    # Subtract original velocities * time_overlap
    a1.x -= time_overlap * a1.vel[0]
    a1.y -= time_overlap * a1.vel[1]
    a2.x -= time_overlap * a2.vel[0]
    a2.y -= time_overlap * a2.vel[1]

    # Calculate collision mechanics at the actual point of collision
    # Line between two atom centerpoints
    col_vector = [a1.x - a2.x,a1.y - a2.y]
    col_dist = ((col_vector[0]**2 + col_vector[1]**2)**0.5)
            
    # Angle between collision vector and velocity difference
    theta = vector_angle(vel_diff, col_vector)
    t_mag = vel_diff_mag * math.cos(theta)
    vel_transfer = [col_vector[0]*t_mag/col_dist,col_vector[1]*t_mag/col_dist]
    
    # Coeffs for adjusting velocity transfer based on mass, 1D collision
    # with the second object stationary relative to velocity transfer
    m1_ratio2 = 2 * a2.mass / (a1.mass + a2.mass)
    m2_ratio1 = 2 * a1.mass / (a1.mass + a2.mass)

    # Calculate new velocity vector
    a1.vel[0] -= m1_ratio2 * vel_transfer[0]
    a1.vel[1] -= m1_ratio2 * vel_transfer[1]
    a2.vel[0] += m2_ratio1 * vel_transfer[0]
    a2.vel[1] += m2_ratio1 * vel_transfer[1]
    
    # Add new velocities * time_overlap
    a1.x += time_overlap * a1.vel[0]
    a1.y += time_overlap * a1.vel[1]
    a2.x += time_overlap * a2.vel[0]
    a2.y += time_overlap * a2.vel[1]


def boundary_collision(atom, boundary):
    """Calculate elastic collisions off of the boundary walls.
    Takes an atom object and a rectangular boundary[xMin,xMax,yMin,yMax]
    """
    # Left wall
    if atom.x < boundary[0] + atom.rad :
        atom.vel[0] = -atom.vel[0]
        # Adjust posiiton for how much the atom should have bounced
        atom.x += 2*(boundary[0] + atom.rad - atom.x)
    # Right wall
    elif atom.x > boundary[1] - atom.rad:
        atom.vel[0] = -atom.vel[0]
        atom.x += 2*(boundary[1] - atom.rad - atom.x)
    # Top wall
    if atom.y < boundary[2] + atom.rad:
        atom.vel[1] = -atom.vel[1]
        atom.y += 2*(boundary[0] + atom.rad - atom.y)
    # Bottom wall
    elif atom.y > boundary[3] - atom.rad:
        atom.vel[1] = -atom.vel[1]
        atom.y += 2*(boundary[1] - atom.rad - atom.y)


def sections_build(num_x,num_y):
    """Creates a 2D array of size num_x by num_y"""
    sects = []
    for i in range(num_x + 1):
        sects.append([])
        for j in range(num_y + 1):
            sects[i].append([])

    return sects


def loc_coefficients(atom, boundary, base_dim, overlap_dim):
    """Defines location coefficients based on section definitions"""
    coeffs = [[],[]]
    # Base atom locations
    coeffs[0].append((atom.x - boundary[0]) / base_dim)
    coeffs[1].append((atom.y - boundary[2]) / base_dim)
    # Positions to check to see if atoms are in neighboring overlaps
    coeffs[0].append(((atom.x-boundary[0]) - overlap_dim) / base_dim)
    coeffs[0].append(((atom.x-boundary[0]) + overlap_dim) / base_dim)
    coeffs[1].append(((atom.y-boundary[2]) - overlap_dim) / base_dim)
    coeffs[1].append(((atom.y-boundary[2]) + overlap_dim) / base_dim)

    return coeffs


def sections_allocate(atom, sections, loc_x1, loc_x2,
                      loc_x3, loc_y1, loc_y2, loc_y3):
    """
    Allocates atom to collision sections given location coefficients
    """
    # Add to base section
    sections[int(loc_x1)][int(loc_y1)].append(atom)       
    # Add to left 3 overlapping sections if located within them
    if int(loc_x2) != int(loc_x1):
        sections[int(loc_x2)][int(loc_y1)].append(atom)
        if int(loc_y2) != int(loc_y1):
            sections[int(loc_x2)][int(loc_y2)].append(atom)
        elif int(loc_y3) != int(loc_y1):
            sections[int(loc_x2)][int(loc_y3)].append(atom)

    # Add to right 3 overlapping sections if located within them
    elif int(loc_x3) != int(loc_x1):
        sections[int(loc_x3)][int(loc_y1)].append(atom)
        if int(loc_y2) != int(loc_y1):
            sections[int(loc_x3)][int(loc_y2)].append(atom)
        elif int(loc_y3) != int(loc_y1):
            sections[int(loc_x3)][int(loc_y3)].append(atom)

    # Add to top and bottom middle overlap sections if located in them
    if int(loc_y2) != int(loc_y1):
        sections[int(loc_x1)][int(loc_y2)].append(atom)
    elif int(loc_y3) != int(loc_y1):
        sections[int(loc_x1)][int(loc_y3)].append(atom)
