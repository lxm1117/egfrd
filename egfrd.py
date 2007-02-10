#!/usr/env python


import math
import bisect
import random
import sys

import numpy
import scipy
import scipy.optimize


from utils import *
from surface import *
import gfrdfunctions
import _gfrd

from gfrdbase import *

"""
class DistanceMatrix:

    def __init__( self, numSpecies ):
        row = [ numpy.array([], numpy.floating ), ] * numSpecies
        self.matrix = [ row, ] * numSpecies

    def __getitem__( self, si1, i1, si2, i2 ):
        return self.matrix[si1][si2][
"""



class EGFRDSingle:
    def __init__( self, sim, si, i ):

        #FIXME: si and i are fragile
        self.particle = Particle( si, i )
        self.sim = sim
        self.dt = 0.0
        self.closest = (-1, -1)

    def fire( self ):
        dt = self.sim.fireSingle( self )
        return dt

    def update( self, t ):
        print 'update'

    def isDependentOn( self, event ):
        #print event
        return False

    def setDr( self, dr ):
        self.sim.getSpeciesByIndex(self.particle.si).pool.drs[self.particle.i] = dr

    def getDr( self ):
        return self.sim.getSpeciesByIndex(self.particle.si).pool.drs[self.particle.i]


class EGFRDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):
        GFRDSimulatorBase.__init__( self )

        self.isDirty = True

        self.scheduler = _gfrd.EventScheduler()

        self.dtMax = INF
        self.dt = INF

        self.pairList = []
        self.singleList = []

        self.lastEvent = None
        self.eventCounter = 0


    def initialize( self ):

        self.scheduler.clear()

        self.initializeSingleList()
        self.initializeSingleDrs()

        #debug
        self.checkShells()

        for single in self.singleList:
            dt = self.calculateFirstPassageTime1( single )
            self.scheduler.addEvent( self.t + dt, single )

        self.scheduler.updateAllEventDependency()

        self.isDirty = False



    def step( self ):

        if self.isDirty:
            self.initialize()

        # if the same single stepped in the last n steps,
        # reinitialize everything.
        # FIXME: don't need to initialize everything.
        #        (1) recalculate dr to the closest with
        #            its dr shrunken.
        #        (2) new dr is
        #            min( dr to the second closest,
        #                 dr to the closest, with its dr
        #                 shrunken )

        nextEvent = self.scheduler.getTopEvent()[1]
        #print self.eventCounter, nextEvent, self.lastEvent
        if self.lastEvent is nextEvent:
            self.eventCounter += 1
            if self.eventCounter >= 10:
                print 'reinitialize'
                self.eventCounter = 0
                self.initialize()
                nextEvent = self.scheduler.getTopEvent()[1]
        else:
            self.eventCounter = 0

        self.lastEvent = nextEvent

        self.t = self.scheduler.getTime()

        self.scheduler.step()

        self.dt = self.scheduler.getTime() - self.t
        
        print 'maxdt', self.dtMax, 'dt', self.dt,\
              'reactions', self.reactionEvents,\
              'rejected moves', self.rejectedMoves,\
              'event counter', self.eventCounter
        print ''
        


    def calculateFirstPassageTime1( self, single ):
        
        species = self.getSpeciesByIndex( single.particle.si )
        fpgf = _gfrd.FirstPassageGreensFunction( species.D )
        dr = single.getDr()
        rnd = random.random()
        print dr
        return fpgf.drawTime( rnd, dr )

    def fireSingle( self, single ):

        species = self.getSpeciesByIndex( single.particle.si )

        displacementS = numpy.array( ( single.getDr(),
                                       numpy.random.uniform( 0.0, Pi ),
                                       numpy.random.uniform( 0.0, 2*Pi ) ) )
        displacement = sphericalToCartesian( displacementS )

        pos = species.pool.positions[single.particle.i]
        pos += displacement
        print 'displacement', length(displacement), single.getDr()

        closest, newdr = self.checkClosestShell( single )
        print 'newdr', newdr
        if newdr <= 0:
            print single.closest, closest
            raise 'Fatal newdr <= 0'
        single.setDr( newdr )
        single.closest = closest

        #debug
        self.checkShells()

        return self.calculateFirstPassageTime1( single )


    def initializeSingleList( self ):

        self.singleList = []

        for si in range( len( self.speciesList ) ):
            species = self.getSpeciesByIndex( si )
            for i in range( species.pool.size ):
                single = EGFRDSingle( self, si, i )
                single.setDr( 0.0 )
                self.singleList.append( single )


    def initializeSingleDrs( self ):

        for single in self.singleList:
            closest, dr = self.checkClosest( single )
            #species = self.getSpeciesByIndex( single.particle.si )
            print 'dr', dr
            #FIXME: take different D into account
            single.setDr( dr * .5 )
        
        #FIXME: here perhaps a second pass to get drs maximally large.


    def checkShells( self ):
        for single in self.singleList:
            dr = single.getDr()
            closest, distance = self.checkClosestShell( single )

            if dr - distance >= 1e-18:
                print single.particle, closest, dr, distance, dr - distance
                raise 'Fatal: shells overlap.'


    def checkClosest( self, single ):

        speciesIndex = single.particle.si
        particleIndex = single.particle.i

        closest = ( -1, -1 )
        minDistance = INF

        species1 = self.getSpeciesByIndex( speciesIndex )
        positions = species1.pool.positions
        position1 = positions[ particleIndex ].copy()

        if self.getReactionType2( species1, species1 ) != None \
           and len( positions ) >= 2 and species1.D != 0.0:

            # temporarily displace the particle
            positions[particleIndex] = NOWHERE

            closest, minDistance = self.checkClosestInSpecies( position1,\
                                                               positions,\
                                                               species1,
                                                               species1 )
            # don't forget to restore the particle position.
            positions[particleIndex] = position1


        for speciesIndex2 in range( speciesIndex )\
                + range( speciesIndex + 1, len( self.speciesList ) ):
            species2 = self.getSpeciesByIndex( speciesIndex2 )

            # empty
            if species2.pool.size == 0:
                continue

            # non reactive
            if self.reactionTypeList2.get( ( species1, species2 ), None )\
                   == None:
                continue
            
            # both of species 1 and 2 are immobile
            if species1.D == 0.0 and species1.D == species2.D:
                continue
                    
            positions = species2.pool.positions

            closest, distance = self.checkClosestInSpecies( position1,
                                                               positions,
                                                               species1,
                                                               species2 )

            if minDistance > distance:
                minDistance = distance
                closest = ( speciesIndex2, minIndex ) 

        return closest, minDistance


    def checkClosestShell( self, single ):

        speciesIndex = single.particle.si
        particleIndex = single.particle.i

        closest = ( -1, -1 )
        minDistance = INF

        species1 = self.getSpeciesByIndex( speciesIndex )
        positions = species1.pool.positions
        position1 = positions[ particleIndex ].copy()

        if self.getReactionType2( species1, species1 ) != None \
           and len( positions ) >= 2 and species1.D != 0.0:

            # temporarily displace the particle
            positions[particleIndex] = NOWHERE

            minIndex, minDistance = self.checkClosestShellInSpecies( position1,
                                                                     positions,
                                                                     species1,
                                                                     species1 )
            closest = ( speciesIndex, minIndex )

            # don't forget to restore the particle position.
            positions[particleIndex] = position1


        for speciesIndex2 in range( speciesIndex )\
                + range( speciesIndex + 1, len( self.speciesList ) ):
            species2 = self.getSpeciesByIndex( speciesIndex2 )

            # empty
            if species2.pool.size == 0:
                continue

            # non reactive
            if self.reactionTypeList2.get( ( species1, species2 ), None )\
                   == None:
                continue
            
            # both of species 1 and 2 are immobile
            if species1.D == 0.0 and species1.D == species2.D:
                continue
                    
            positions = species2.pool.positions

            minIndex, distance = self.checkClosestShellInSpecies( position1,
                                                                  positions,
                                                                  species1,
                                                                  species2 )
            if minDistance > distance:
                minDistance = distance
                closest = ( speciesIndex2, minIndex ) 

        return closest, minDistance


    def checkClosestInSpecies( self, position1, positions2,
                               species1, species2 ):

        # calculate distances
        distanceList = self.distanceSqArray( position1, positions2 )
        distanceList = numpy.sqrt( distanceList )
        
        # find the closest particle ignoring their shells.
        minIndex = distanceList.argmin()
        
        radius12 = species1.radius + species2.radius
        minDistance = distanceList[minIndex] - radius12
        
        return minIndex, minDistance

    def checkClosestShellInSpecies( self, position1, positions2,
                                    species1, species2 ):

        # calculate distances
        distanceList = self.distanceSqArray( position1, positions2 )
        distanceList = numpy.sqrt( distanceList )
        
        # find the particle of which shell is closest to this particle.
        distanceList -= species2.pool.drs
        minIndex = distanceList.argmin()
        
        radius12 = species1.radius + species2.radius
        minDistance = distanceList[minIndex] - radius12
        
        return minIndex, minDistance




'''
    def checkDistance( self, position1, positions2, species1, species2, n=1 ):

        #positions2 = species2.pool.positions

        distanceSq = self.distanceSqArray( position1, positions2 )
        sortedindices = distanceSq.argsort()

        #debug
        radius12 = species1.radius + species2.radius
        radius12sq = radius12 * radius12

        # check if particles overlap
        if distanceSq[ sortedindices[0] ] < radius12sq - 1e-20 and \
               distanceSq[ sortedindices[0] ] != 0.0:
            print position1, positions2[ sortedindices[0] ]
            print 'dr<radius', math.sqrt(distanceSq[sortedindices[0]]), radius12
            print species1.id, species2.id, sortedindices[0]
            raise "critical"

        distanceSqSorted = distanceSq.take( sortedindices[:n] )
        distances = numpy.sqrt( distanceSqSorted )
        # instead of just

        #drs = self.H * numpy.sqrt( 6.0 * species1.D * dts )
        drs = ( distances - radius12 )

        indices = sortedindices[:n]

        return indices, drs


'''
