#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 17:55:59 2022

@author: james
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import Body as bd

class Circles(object):
    
    #Constructor of Circles class
    #Param: bodies | a list of type Body
    def __init__(self, bodies):
        self.bodies = bodies
        
        
    def init(self):
        return self.patches
        
    def animate(self, i): #i is the frame number
            #loops through every body
            #updates the position of each planet in accordance to the frame no.
            for x in range(len(self.bodies)):
                self.patches[x].center = (self.bodies[x].positions[i][0], 
                                          self.bodies[x].positions[i][1])

            #returns the now updated list of patches in their new position
            return self.patches
        
    def display(self):
        fig = plt.figure() #creates a figure
        ax = plt.axes() #creaates an axes
        
        #creates empty list of Patches
        #each patch represents one celestial body
        self.patches = [] 
        
        #loops through all the bodies and appends a circle containing its starting
        #position to the patches list
        for x in range(len(self.bodies)):
            self.patches.append(plt.Circle((self.bodies[x].positions[0][0], self.bodies[x].positions[0][0]), self.bodies[x].size, color = self.bodies[x].colour, animated = True))
        
        
        #loops through all patches and adds it to the axes
        for x in range(len(self.patches)):
            ax.add_patch(self.patches[x])
        
        # <--- SETS SCALE OF PLOT --->
        ax.axis('scaled')
        ax.set_xlim(-3e11,3e11)
        ax.set_ylim(-3e11,3e11)
        # <-------------------------->
    
        
        numFrames = len(self.bodies[0].positions)
        self.anim = FuncAnimation(fig, self.animate, frames = numFrames, 
                                  repeat = True, interval = 1, blit = True)
        print(type(self.anim))
        #self.animate is the name of the method causing the animation
        plt.show()
        
def main():
   print("ran")
    
    
       
if __name__ == '__main__':
    main()
   
        
        