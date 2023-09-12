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
    def __init__(self, positions, positions2, positions3, positions4, positions5):
        self.positions = positions
        self.positions2 = positions2
        self.positions3 = positions3
        self.positions4 = positions4
        self.positions5 = positions5

        
# =============================================================================
#     def runSim(self):
#         theta = np.linspace(0, 2*np.pi, 250)
#         self.xpos = np.cos(theta) #in ours, x pos is just the first element in the array
#         self.ypos = np.sin(theta) #in ours, y pos is just the 2nd elemnt in the array
# =============================================================================
        
    def init(self):
        return self.patch, self.patch2, self.patch3, self.patch4, self.patch5
        
    def animate(self, i): #i is the frame number
# =============================================================================
#         for x in range(len(self.bodies)):
#             currentBody = self.bodies[x]
#             self.patches[x].center = (currentBody.positions[x][i], currentBody.positions[x][i])
# =============================================================================
           # self.patches[0].center = (self.bodies[1].positions[i][0], self.bodies.positions[i][1])
            self.patch.center = (self.positions[i][0], self.positions[i][1])
            self.patch2.center = (self.positions2[i][0], self.positions2[i][1])
            self.patch3.center = (self.positions3[i][0], self.positions3[i][1])
            self.patch4.center = (self.positions4[i][0], self.positions4[i][1])
            self.patch5.center = (self.positions5[i][0], self.positions5[i][1])


            return self.patch, self.patch2, self.patch3, self.patch4, self.patch5
        
    def display(self):
        fig = plt.figure()
        ax = plt.axes()
        
        self.patches = []
        
        
        
        self.patch = plt.Circle((self.positions[0][0], self.positions[0][1]), 1e10, color='b', animated = True)
        self.patch2 = plt.Circle((self.positions2[0][0], self.positions2[0][1]), 0.6e10, color='#9e0909', animated = True)
        self.patch3 = plt.Circle((self.positions3[0][0], self.positions3[0][1]), 1.5e10, color='black', animated = True)
        self.patch4 = plt.Circle((self.positions4[0][0], self.positions4[0][1]), 0.4e10, color='0.6', animated = True)
        self.patch5 = plt.Circle((self.positions5[0][0], self.positions5[0][1]), 0.8e10, color='#85640c', animated = True)

        ax.add_patch(self.patch)
        ax.add_patch(self.patch2)
        ax.add_patch(self.patch3)
        ax.add_patch(self.patch4)
        ax.add_patch(self.patch5)

        
        ax.axis('scaled')
        ax.set_xlim(-3e11,3e11)
        ax.set_ylim(-3e11,3e11)
    
        numFrames = len(self.positions)
        self.anim = FuncAnimation(fig, self.animate, frames = numFrames, repeat = True, interval = 20, blit = True)
        #self.animate is the name of the method causing the animation
        plt.show()
        
def main():
   print("ran")
    
    
       
if __name__ == '__main__':
    main()
   
        
        