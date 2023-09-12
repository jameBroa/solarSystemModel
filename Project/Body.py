#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 12:07:39 2022

@author: james
"""

import numpy as np
from numpy.linalg import norm
import Animation as anim
import Simulation as sim
import math
import matplotlib.pyplot as plt

G = 6.67408*(10**-11) 
listOfBodies = []


class Body(object):
    
    #<-- Constructor for Body object -->
    #name: name of the body
    #mass: mass of the body
    #position: radius of the body relative to sun
    #colour: colour of body (for simulation)
    #size: size of body (for simulation)
    def __init__(self, name, mass, position, colour, size):
        self.name = name
        self.mass = mass
        self.position = np.array([position, 0]) 
        self.positions = [] 
        if(position != 0):
            #calculates the initial velocity of the body based on relative 
            #position to Sun
            self.currentVelocity = np.array([0, math.sqrt((G*1.98e30)/position)])
            
        else:
            #Else statement is called when input is 0
            #this is a special case for the suns instantiation
            self.currentVelocity = np.array([0,0])
        
        self.newVelocity = 0
        
        self.currentAcceleration = 0
        self.previousAcceleration = 0
        self.newAcceleration = 0
        
        
        self.orbitalPeriod = 0
        self.kineticEnergy = 0.5*mass*(self.currentVelocity**2)
        self.potentialEnergy = 0
        
        self.colour = colour
        self.size = size
        
        
        #This function calculates the current acceleration of all bodies
        #This is performed before the first timestep
    def calculateInitialAccelerations(self):
        accelerationSum = 0 #creates sum variable for total acceleration
        for x in range(len(listOfBodies)): #loops through list of bodies
            #prevents finding acceleration with itself
            if(listOfBodies[x] != self):
                #calculates vector position relative to other body
                vectorPosition = self.position - listOfBodies[x].position
                #calculates new acceleration
                acceleration = -1*G * ((listOfBodies[x].mass)/
                                       (norm(vectorPosition)**3) * vectorPosition)
                #adds acceleration between two bodies total acceleration
                accelerationSum += acceleration 
        #since this sets initial acceleration, both previousAcceleration and 
        #current acceleration are set to the Sum
        self.previousAcceleration = accelerationSum #sets a(t-dt)
        self.currentAcceleration = accelerationSum #sets a(t)
        
        
        #This function performs exactly the same as 
        #calculateInitialAccelerations() except this is used to set
        #the new acceleration or a(t+dt)
    def calculateAccelerations(self):
        accelerationSum = 0
        for x in range(len(listOfBodies)):
            if(listOfBodies[x] != self):
                vectorPosition = self.position - listOfBodies[x].position
                acceleration = -1*G * ((listOfBodies[x].mass)/(norm(vectorPosition)**3) 
                                       * vectorPosition)
                accelerationSum = accelerationSum + acceleration
        
        self.newAcceleration = accelerationSum #updates a(t+dt)
       
    
    
        #This function calculates the new velocity given current velocity
        #using the Beeman integration scheme
    def calculateNewVelocity(self, timeStep):
         newVelocity = self.currentVelocity + (1/6)*(
             2*self.newAcceleration + 
             (5*self.currentAcceleration) - self.previousAcceleration)*timeStep
         self.newVelocity = newVelocity #sets v(t+dt)
         
         #This function calculates the new position using Beeman integration
    def calculateNewPosition(self, timeStep):
        position = self.position + self.currentVelocity*timeStep + (1/6)*(
            4*self.currentAcceleration-self.previousAcceleration) * timeStep**2
        self.position = position #updates position | sets r(t+dt)
        self.positions.append(position) #appends position to position list of body
        
        #This function is resposible for adjusting previous/current accelerations
        #or velocities. This is so in the next timestep, the same calculation
        #isn't performed
    def resetVariables(self):
        self.previousAcceleration = self.currentAcceleration #sets a(t-dt) = a(t)
        self.currentAcceleration = self.newAcceleration #sets a(t) = a(t+dt)
        self.currentVelocity = self.newVelocity #sets v(t) = v(t+dt)
        
        #This function calculates the orbital period of the body
        #A prequisite is that at least one orbit is simulated after a set no.
        #of timesteps have ocurred
    def calculateOrbitalPeriod(self, timeStep):
        #variable which will hold what timestep the orbit is complete
        numIterForOnePeriod = 0 
        
        for x in range(len(self.positions)):
            #necessary if statement to prevent list out of range
            if (x < (len(self.positions))-1):
                #if statement checks to see when the y-position of the object 
                #is negative and the next one is positive. This implies a full
                #orbit is completed.
                if(self.positions[x][1] <= 0 and self.positions[x+1][1] >= 0):
                    #sets the number of timesteps until this is achieved
                    numIterForOnePeriod = x 
                    #returns the number of timesteps x timeStep which gives
                    #the time taken in SECONDS for one orbit.
                    return (numIterForOnePeriod * timeStep)
                
        #Calculates kinetic energy of the body 
    def calculateKineticEnergy(self):
        self.kineticEnergy = 0.5*self.mass*(norm(self.currentVelocity)**2)
        
    
    #calculates Potential energy of the entire system 
def calculatePotentialEnergy():
    #initializes variable for sum
    potentialEnergySum = 0 
    #loops through list of bodies
    for x in range(len(listOfBodies)):
        currentBody = listOfBodies[x]
        sum = 0
        #sets another variable for sum
        #we do this because we find the potential of one body in relation to 
        #all the others and sums it
        #we do this for all bodies
        for y in range(len(listOfBodies)):
            #checks to see its not finding potential with itself
            if(listOfBodies[y] != currentBody):  
                #calcualtes distance between two bodies
                distance = norm(listOfBodies[y].position - currentBody.position)
                #calculates potential between bodies
                calc = (G * currentBody.mass * listOfBodies[y].mass)/distance
                #adds it to sum
                sum += calc
        #adds potential of that body with all other bodies to solar system sum
        potentialEnergySum += sum 
        
    #returns half of sum due to double counting
    return -0.5*potentialEnergySum 
                
            
            
# <-- initializes empty lits for energies --> 
kineticEnergies = []
potentialEnergies = []
totalEnergies = []
# <----------------------------------------->

def performTimeStep(timeStep, iterCount): #performs one timesteo

    for x in range(len(listOfBodies)): #updates position
        currentBody = listOfBodies[x]
        currentBody.calculateNewPosition(timeStep)
    
    for x in range(len(listOfBodies)): #updates accelerations
        currentBody = listOfBodies[x]
        currentBody.calculateAccelerations()
    
    for x in range(len(listOfBodies)): #updates velocities
        currentBody = listOfBodies[x]
        currentBody.calculateNewVelocity(timeStep)
        
    kineticEnergySum = 0 #sum for total K.E. in solar system
    
    for x in range(len(listOfBodies)): #updates kinetic energy and sums it
        currentBody = listOfBodies[x]
        currentBody.calculateKineticEnergy()
        kineticEnergySum += currentBody.kineticEnergy
        
        
    #appends kinetic energy list of the kinetic energy of the system
    #at a particular timestep
    kineticEnergies.append(np.array([iterCount, kineticEnergySum]))

    #calculates potential energy and appends value to list of potential energies
    #given the timestep
    potentialEnergy = calculatePotentialEnergy()
    potentialEnergies.append(np.array([iterCount, potentialEnergy]))
    
    #calculates total energy and appends value to list of potential energies
    #given the timestep
    totalEnergy = kineticEnergySum + potentialEnergy
    totalEnergies.append(np.array([iterCount, totalEnergy]))
    
        
    for x in range(len(listOfBodies)): #update new/old/current variables
        currentBody = listOfBodies[x]
        currentBody.resetVariables()
        
        
    #setsUp all bodies accelerations before starting starting simulation
    #sets up a(t) and a(t-dt)
def setupBodies():
    for x in range(len(listOfBodies)):
        listOfBodies[x].calculateInitialAccelerations()
        
    #Prints position of all bodies at the time of function call
def printAllPositions():
    print("---")
    for x in range(len(listOfBodies)):
        print(str(listOfBodies[x].name) + " is in position: " + 
              str(listOfBodies[x].position))
    print("---")
    
    #Peforms 'numIterations' number of timeSteps at the call of one function
def performAllTimeSteps(numIterations, timeStep):
    for x in range(numIterations):
        performTimeStep(timeStep, x)
    with open("energyInfo.txt", "w") as f:
        for x in range(len(totalEnergies)):
            f.write(str(totalEnergies[x]))
            
    f.close()
        
    #Calculates orbital period of all bodies in the solar system    
def calculateOrbitalPeriods(timeStep):
    #loops through all bodies and calcualtes orbital period
    for x in range(len(listOfBodies)):
        period = listOfBodies[x].calculateOrbitalPeriod(timeStep)
        listOfBodies[x].orbitalPeriod = period
        
    #Prints orbital periods of all bodies in relation to a particular time
    #In this case this is the time taken for one Earth orbit
def printOrbitalPeriods(referenceTime):
    for x in range(len(listOfBodies)):
        time = (listOfBodies[x].orbitalPeriod) / referenceTime
        if(referenceTime == 84600):
            print("Orbital period of " + listOfBodies[x].name + " is: " +  '%.2f'%time + 
                  " Earth days.")
        elif(referenceTime == 84600*365):
            print("Orbital period of " + listOfBodies[x].name + " is: " +  '%.2f'%time + 
                  " Earth years.")
    
    #Finds the initial velocity of the satellite which gets closest to Mars
    #input parameter is the body in question
def findVelocityForClosestDistance(sat, timeStep):
    #initializes a list of velocities from the calcaulted initial velocity to 30,000
    
    velocities = np.linspace(sat.currentVelocity, np.array([0,30000]), 10)
    
    #creates variable minimum distance set to a large number
    #this is so the first calculated distance will definately be smaller
    #this variable is for the minimumDistance over all velocities
    minimumDistance = 10e12
    minimumDistanceEarth = 10e12

    #initializes idealVelocity to be 0.
    idealVelocity = 0
    #sets time
    numIterationsToReachClosestD = 0
    
    #loops through all possible velocities
    for x in range(len(velocities)):
        #sets the current velocity to a velocity in the velocities list
        sat.currentVelocity = velocities[x]
        #adds satellite to listOfBodies
        listOfBodies.append(sat)
        print("current v: " + str(velocities[x]))
        
        # <-- sets params for simulation -->
        numIterations = 5*365
        timeStep = 86400
        # <-------------------------------->
        
        #sets variable for minimum
        #this minimum is for finding the minimum distance for a given velocity
        minimum = 10e11

        
        setupBodies() #sets up initial accelerations of bodies
        
        #loops through and performs 'numIterations' number of timesteps
        #at this given initial velocity and finds the closest distance
        for y in range(numIterations):
            performTimeStep(timeStep, y)
            vectorPosition = listOfBodies[4].position - listOfBodies[5].position
            vectorPositionToEarth = listOfBodies[3].position - listOfBodies[5].position
            distance = norm(vectorPosition)
            
            #finds closest distance of satellite at this given velocity
            if(distance < minimum):
                print(distance)
                minimum = distance
                print("temp: " + str(y))
                if(minimum < minimumDistance):
                    numIterationsToReachClosestD =y
            
                
            
        #if the closest distance at this given velocity is smaller than previous
        #distances at different velocities, it beocmes the new minimum
        if(minimum < minimumDistance):
            minimumDistance = minimum
            idealVelocity = velocities[x] #sets the new ideal velocity
     
            
    
    print("This is the closest overrall distance: " + str(minimumDistance))
    print("This is the ideal velocity: " + str(idealVelocity))
    print("This is how long is took to get the closest distance: " + str((numIterationsToReachClosestD * timeStep)/86400))
    print("This is ho")
    #returns ideal velocity
    return idealVelocity
    
    #function creates graphs of different types of Energy over time for the system
    #input parameter is the type of graph it outputs
def createGraphs(graphNo):
    #--kinetic energy--
    xpos =[]
    ypos = []
    #------------------
    #--potential energy--
    xpos2 = []
    ypos2 = []
    #--------------------
    #--total energy --
    xposTotalEnergy = []
    yposTotalEnergy = []
    #-----------------
    
    #loops through all data points for energies and adds each data point to coordinate
    #lists
    #since all energies are appended in the same function, the len(kineticEnergies)
    #is the same as len(potentialEnergies) and len(totalEnergies)
    for x in range(len(totalEnergies)):
        if (x > 5):
            xpos.append(kineticEnergies[x][0])
            ypos.append(kineticEnergies[x][1])
        
            xpos2.append(potentialEnergies[x][0])
            ypos2.append(potentialEnergies[x][1])
        
            xposTotalEnergy.append(totalEnergies[x][0])
            yposTotalEnergy.append(totalEnergies[x][1])
            
            
    #only outputs the graph of the one selected on input
    if(graphNo == 1):
        plt.plot(xpos, ypos, color='r', label='Kinetic Energy')
        plt.show()
    elif(graphNo == 2):       
        plt.plot(xpos2,ypos2, color='b', label='Potential Energy')
        plt.show()
    elif(graphNo == 3):    
        plt.plot(xposTotalEnergy, yposTotalEnergy, color = 'g', label='Total Energy')
        plt.xlabel("Time (Earth Days)")
        plt.ylabel("Total Energy (J)")
        plt.show()
    
    #creates Body objects from text file
def setUpPlanets():
    f = open("planetInfo.txt", "r") #opens text file and reads by line
    for x in f: #loops through all lines
        temp = x.split(",") #splits each line by "," and is stored in a list
        #creates Body object
        celestialBody = Body(temp[0], float(temp[1]), float(temp[2]), temp[3], 
                             float(temp[4]))
        #adds Body to list of bodies
        listOfBodies.append(celestialBody)
    #closes text file
    f.close()
    
    #creates Satellite from text file
    #functions same as setUpPlanets but for one body
def addSatelliteWithNormalV():
    f = open("satelliteInfo.txt", "r")
    for x in f:
        temp = x.split(",")
    celestialBody = Body(temp[0], float(temp[1]), float(temp[2]), temp[3], 
                         float(temp[4]))
    listOfBodies.append(celestialBody)
    f.close()
    
    #Creates Satellite from text file but sets the current velocity to be
    #velocity which gets it closest to Mars
def addSatelliteWithOptimalV(timeStep):
    f = open("satelliteInfo.txt", "r")
    for x in f:
        temp = x.split(",")
    celestialBody = Body(temp[0], float(temp[1]), float(temp[2]), temp[3], 
                         float(temp[4]))
    newV = findVelocityForClosestDistance(celestialBody, timeStep) #calcualtes ideal velocity
    celestialBody.currentVelocity = newV #sets satellites velocity to the ideal velocity
    listOfBodies.append(celestialBody) 
    f.close()#closes text file
       
    #function returns whether the body gets close to Earth
def returnsToEarth(body):
    #loops thorugh all the positions of Earth
    for x in range(len(listOfBodies[3].positions)):
        #Sets Earth position at this timestep
        earthPosition = np.array([listOfBodies[3].positions[x][0], 
                                  listOfBodies[3].positions[x][1]])
        #Sets Satellite position at this timestep
        satellitePosition = np.array([body.positions[x][0], body.positions[x][1]])
        #calculates distance between bodies
        distanceToEarth = norm(earthPosition - satellitePosition)
        if(distanceToEarth < 100000):
            print("It returns back to Earth")
            return True
    print("It doesn't return to Earth")
    return False
        

def main():
    timeStep = 86400
    
    setUpPlanets() 
    
    #To add the satellite, uncomment one of the function below
    
    #addSatelliteWithNormalV()
    #addSatelliteWithOptimalV(timeStep)

    
   
    
    performAllTimeSteps(5*365, 84600)
    
    calculateOrbitalPeriods(84600)
    printOrbitalPeriods(84600)
    
    
    #To run the graph function, comment out the simulation code
    
    #createGraphs(3)
    
    
    simulation = sim.Circles(listOfBodies)
    simulation.display()
    
    
        
if __name__ == '__main__':
    main()
    

