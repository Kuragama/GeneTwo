__author__ = 'Evan'

import random
from PIL import Image
from osgeo import gdal
import numpy
from enum import Enum
from operator import itemgetter

class CrossoverType(Enum):
    SINGLEPOINT = 0
    DOUBLEPOINT = 1
    DOUBLEPOINTBIASED = 2
    UNIFORM = 3
    NONE = 4

class SelectionType(Enum):
    RANDOMELITE = 0
    CASCADINGELITE = 1
    ROULETTEELITE = 2
    ROULETTETOTAL = 3

class TerminationType(Enum):
    CONVERGENCE = 0
    OPTIMUMCONVERGENCE = 1
    NATURAL = 2

class Algorithm:

    NUMCHROMOSOMES = 2

    def __init__(self, numGenes, popSize, maxGen, crossoverType, crossoverChance, mutationSize, mutationChance, numElite, selectionType, terminationType,
                 optimumX, optimumY, convergenceConfidence, evaluationDictionary):
        self.NUMGENES = numGenes
        self.POPSIZE = popSize
        self.MAXGEN = maxGen
        self.CROSSOVERTYPE = crossoverType
        self.CROSSOVERCHANCE = crossoverChance
        self.MUTATIONSIZE = mutationSize
        self.MUTATIONCHANCE = mutationChance
        self.NUMELITE = numElite
        self.SELECTIONTYPE = selectionType
        self.TERMINATIONTYPE = terminationType
        self.OPTIMUMX = optimumX
        self.OPTIMUMY = optimumY
        self.CONVERGENCECONFIDENCE = convergenceConfidence
        self.__population = Population(self.NUMCHROMOSOMES, self.NUMGENES, self.POPSIZE)
        self.__mutator = Mutator()
        self.__evaluator = Evaluator(evaluationDictionary)

        self.__gens = 0
        self.__convergedGens = 0
        self.__bestEval = 0
        self.__bestEvalGen = 0

    def aggregateChromosomes(self, member):
        n = 0
        for i in range(self.NUMCHROMOSOMES):
            n += member.getChromosome(i) << (i * self.NUMGENES)
        return n

    def separateChromosomes(self, n):
        chromosomes = []
        mask = (1 << self.NUMGENES) - 1
        for i in range(self.NUMCHROMOSOMES):
            chromosomes.append(n & mask)
            n >>= self.NUMGENES
        return chromosomes

    def nextGen(self):
        self.__sortPopulation()
        elitePop = []
        childPop = []
        for i in range(self.NUMELITE):
            elitePop.append(self.__population.getMember(i))
        if (self.CROSSOVERTYPE == CrossoverType.NONE):
            for i in range(self.NUMELITE, self.__population.getSize()):
                childPop.append(self.__population.getMember(i))
        else:
            if (self.SELECTIONTYPE == SelectionType.ROULETTETOTAL):
                childPop = self.__createChildren(self.__population.getMembers())
            else:
                childPop = self.__createChildren(elitePop)
        self.__mutateChildren(childPop)
        self.__population.setMembers(elitePop + childPop)
        return [self.__checkTermination(), self.__bestEval]

    def __sortPopulation(self):
        evalData = []
        for i in range(self.__population.getSize()):
            eval = self.__evaluator.evaluate(self.__population.getMember(i).getChromosome(0),
                                                          self.__population.getMember(i).getChromosome(1),
                                                          self.NUMGENES)
            if (eval > self.__bestEval):
                self.__bestEval = eval
                self.__bestEvalGen = self.__gens
            evalData.append([i, eval])
        evalData = list(reversed(sorted(evalData,key=lambda x: x[1])))
        sortedPop = [None] * self.__population.getSize()
        for i in range(len(evalData)):
            sortedPop[i] = self.__population.getMember(evalData[i][0])
        for i in range(len(sortedPop)):
            self.__population.setMember(i, sortedPop[i])

    def __createChildren(self, parents):
        children = []
        if (self.SELECTIONTYPE == SelectionType.RANDOMELITE):
            while(len(children) < self.POPSIZE - self.NUMELITE):
                p1 = random.randrange(len(parents) - 1)
                p2 = random.randrange(len(parents) - 1)
                while(p2 == p1):
                    p2 = random.randrange(len(parents) - 1)
                children += self.__crossover(parents[p1], parents[p2])
            while(len(children) > self.POPSIZE - self.NUMELITE):
                children.pop(len(children) - 1)
        elif (self.SELECTIONTYPE == SelectionType.CASCADINGELITE):
            i = 0
            while(len(children) < self.POPSIZE - self.NUMELITE):
                children += self.__crossover(parents[i], parents[i + 1])
                i = (i + 2) % (len(parents) - 1)
            while(len(children) > self.POPSIZE - self.NUMELITE):
                children.pop(len(children) - 1)
        elif (self.SELECTIONTYPE == SelectionType.ROULETTEELITE or self.SELECTIONTYPE == SelectionType.ROULETTETOTAL):
            children = self.__rouletteWheelSelection(parents, self.POPSIZE - self.NUMELITE)
        return children

    def __mutateChildren(self, children):
        if (self.MUTATIONSIZE > 0):
            for child in children:
                if (random.randrange(100) < self.MUTATIONCHANCE):
                    points = []
                    for i in range(self.MUTATIONSIZE):
                        points.append(random.randrange(self.NUMGENES * self.NUMCHROMOSOMES))
                    ag = self.__mutator.pointMutation(self.aggregateChromosomes(child), points)
                    chr = self.separateChromosomes(ag)
                    child.setChromosome(0, chr[0])
                    child.setChromosome(1, chr[1])

    def __rouletteWheelSelection(self, parents, numChildren):
        evalSum = 0
        evaluations = []
        for i in range(len(parents)):
            evaluations.append(self.__evaluator.evaluate(parents[i].getChromosome(0), parents[i].getChromosome(1), self.NUMGENES))
            evalSum += evaluations[i]
        children = []
        while(len(children) < numChildren):
            stop = random.randrange(int(evalSum))
            p1 = 0
            indexSum = evaluations[0]
            while (indexSum < stop):
                indexSum += evaluations[p1 + 1]
                p1 += 1
            stop = random.randrange(int(evalSum))
            p2 = 0
            indexSum = evaluations[0]
            while (indexSum < stop):
                indexSum += evaluations[p2 + 1]
                p2 += 1
            #Repeat to avoid selecting same index twice
            while (p2 == p1):
                stop = random.randrange(int(evalSum))
                p2 = 0
                indexSum = evaluations[0]
                while (indexSum < stop):
                    indexSum += evaluations[p2 + 1]
                    p2 += 1
            #print(stop, evalSum, p1, p2, len(parents))
            children += self.__crossover(parents[p1], parents[p2])
        while(len(children) > numChildren):
                children.pop(len(children) - 1)
        return children

    def __crossover(self, parent1, parent2):
        if (random.randrange(100) <= self.CROSSOVERCHANCE):
            if (self.CROSSOVERTYPE == CrossoverType.NONE):
                return [parent1, parent2]
            elif (self.CROSSOVERTYPE == CrossoverType.SINGLEPOINT):
                child1 = Member(self.NUMCHROMOSOMES, self.NUMGENES)
                child2 = Member(self.NUMCHROMOSOMES, self.NUMGENES)
                ag = self.__mutator.pointCrossover(self.aggregateChromosomes(parent1), self.aggregateChromosomes(parent2),
                                       [random.randrange((self.NUMGENES * self.NUMCHROMOSOMES) - 1)],
                                                   self.NUMGENES * self.NUMCHROMOSOMES)
                chr1 = self.separateChromosomes(ag[0])
                chr2 = self.separateChromosomes(ag[1])
                child1.setChromosome(0, chr1[0])
                child1.setChromosome(1, chr1[1])
                child2.setChromosome(0, chr2[0])
                child2.setChromosome(1, chr2[1])
                return [child1, child2]
            elif (self.CROSSOVERTYPE == CrossoverType.DOUBLEPOINT):
                child1 = Member(self.NUMCHROMOSOMES, self.NUMGENES)
                child2 = Member(self.NUMCHROMOSOMES, self.NUMGENES)
                radices = [0, 0]
                radices[0] = random.randrange((self.NUMGENES * self.NUMCHROMOSOMES) - 1)
                radices[1] = random.randrange((self.NUMGENES * self.NUMCHROMOSOMES) - 1)
                while (radices[1] == radices[0]):
                    radices[1] = random.randrange((self.NUMGENES * self.NUMCHROMOSOMES) - 1)
                ag = self.__mutator.pointCrossover(self.aggregateChromosomes(parent1), self.aggregateChromosomes(parent2),
                                       radices, self.NUMGENES * self.NUMCHROMOSOMES)
                chr1 = self.separateChromosomes(ag[0])
                chr2 = self.separateChromosomes(ag[1])
                child1.setChromosome(0, chr1[0])
                child1.setChromosome(1, chr1[1])
                child2.setChromosome(0, chr2[0])
                child2.setChromosome(1, chr2[1])
                return [child1, child2]
            elif (self.CROSSOVERTYPE == CrossoverType.DOUBLEPOINTBIASED):
                child1 = Member(self.NUMCHROMOSOMES, self.NUMGENES)
                child2 = Member(self.NUMCHROMOSOMES, self.NUMGENES)
                chr1 = self.__mutator.pointCrossover(parent1.getChromosome(0), parent2.getChromosome(0),
                                       [random.randrange(self.NUMGENES - 1)], self.NUMGENES * self.NUMCHROMOSOMES)
                chr2 = self.__mutator.pointCrossover(parent1.getChromosome(1), parent2.getChromosome(1),
                                       [random.randrange(self.NUMGENES - 1)], self.NUMGENES * self.NUMCHROMOSOMES)
                child1.setChromosome(0, chr1[0])
                child1.setChromosome(1, chr2[0])
                child2.setChromosome(0, chr1[1])
                child2.setChromosome(1, chr2[1])
                return [child1, child2]
            elif (self.CROSSOVERTYPE == CrossoverType.UNIFORM):
                child1 = Member(self.NUMCHROMOSOMES, self.NUMGENES)
                child2 = Member(self.NUMCHROMOSOMES, self.NUMGENES)
                mask = random.randrange((1 << (self.NUMGENES * self.NUMCHROMOSOMES)) - 1)
                ag = self.__mutator.uniformCrossover(self.aggregateChromosomes(parent1), self.aggregateChromosomes(parent2),
                                                       mask)
                chr1 = self.separateChromosomes(ag[0])
                chr2 = self.separateChromosomes(ag[1])
                child1.setChromosome(0, chr1[0])
                child1.setChromosome(1, chr1[1])
                child2.setChromosome(0, chr2[0])
                child2.setChromosome(1, chr2[1])
                return [child1, child2]
        return []

    def __checkTermination(self):
        if (self.TERMINATIONTYPE == TerminationType.CONVERGENCE):
            dif = self.NUMGENES * (self.CONVERGENCECONFIDENCE / 100)
            minx = miny = maxx = maxy = 0
            m0 = self.__population.getMember(0)
            m1 = self.__population.getMember(1)
            if (m0.getChromosome(0) > m1.getChromosome(1)):
                maxx = m0.getChromosome(0)
                minx = m1.getChromosome(0)
            else:
                maxx = m1.getChromosome(0)
                minx = m0.getChromosome(0)
            if (m0.getChromosome(1) > m1.getChromosome(1)):
                maxy = m0.getChromosome(1)
                miny = m1.getChromosome(1)
            else:
                maxy = m1.getChromosome(1)
                miny = m0.getChromosome(0)
            for i in range(2, self.__population.getSize()):
                m = self.__population.getMember(i)
                if (m.getChromosome(0) > maxx):
                    maxx = m.getChromosome(0)
                elif (m.getChromosome(0) < minx):
                    minx = m.getChromosome(0)
                if (m.getChromosome(1) > maxy):
                    maxy = m.getChromosome(1)
                elif (m.getChromosome(1) < miny):
                    miny = m.getChromosome(1)
            if ((maxx - minx <= dif) and (maxy - miny <= dif)):
                self.__convergedGens += 1
                if (self.__convergedGens == self.CONVERGENCECONFIDENCE):
                    return True
            else:
                self.__convergedGens = 0
        elif (self.TERMINATIONTYPE == TerminationType.OPTIMUMCONVERGENCE):
            dif = self.NUMGENES * (self.CONVERGENCECONFIDENCE / 200)
            for i in range(self.__population.getSize()):
                if (abs(self.__population.getMember(i).getChromosome(0) - self.OPTIMUMX) >= dif or
                            abs(self.__population.getMember(i).getChromosome(1) - self.OPTIMUMY) >= dif):
                    return False
            return True
        self.__gens += 1
        if (self.__gens >= self.MAXGEN):
            return True
            self.__gens = 0
        return False

    def getEvalData(self):
        #Calculate average x, y, and fitness score
        totalX = 0
        totalY = 0
        totalFitness = 0
        bestScore = 0
        bestX = bestY = 0
        for i in range(self.__population.getSize()):
            totalX += self.__population.getMember(i).getChromosome(0)
            totalY += self.__population.getMember(i).getChromosome(1)
            score = self.__evaluator.evaluate(self.__population.getMember(i).getChromosome(0),
                                              self.__population.getMember(i).getChromosome(1), self.NUMGENES)
            if score > bestScore:
                bestScore = score
                bestX = self.__population.getMember(i).getChromosome(0)
                bestY = self.__population.getMember(i).getChromosome(1)
            totalFitness += score
        aveX = totalX / self.__population.getSize()
        aveY = totalY / self.__population.getSize()
        aveFitness = totalFitness / self.__population.getSize()
        return [aveX, aveY, aveFitness, bestScore, bestX, bestY, self.__bestEvalGen]

    def bruteForce(self, width, height):
        return self.__evaluator.bruteForce(self.NUMGENES, width, height)

    def getGenerationsPassed(self):
        return self.__gens

    def getDisplayCoords(self, minx, maxx, miny, maxy):
        coords = []
        for i in range(self.__population.getSize()):
            x = self.__evaluator.chromosomeToCoord(self.NUMGENES, minx, maxx, self.__population.getMember(i).getChromosome(0))
            y = self.__evaluator.chromosomeToCoord(self.NUMGENES, miny, maxy, self.__population.getMember(i).getChromosome(1))
            coords.append([x, y])
        return coords

class Mutator:

    """

    """
    def pointCrossover(self, parent1, parent2, radices, bits):
        for radix in radices:
            mask = (1 << bits) - 1
            lm = mask >> bits - radix
            um = mask << radix
            child1 = (parent1 & lm) | (parent2 & um)
            child2 = (parent1 & um) | (parent2 & lm)
            parent1 = child1
            parent2 = child2
        return [parent1, parent2]

    """

    """
    def uniformCrossover(self, parent1, parent2, mask):
        return [(parent1 & mask) + (parent2 & (~mask)), (parent2 & mask) + (parent1 & (~mask))]

    """

    """
    def pointMutation(self, parent, points):
        for point in points:
            mask = (1 << point)
            if (parent & mask == 0):
                parent += mask
            else:
                parent -= mask
        return parent

class Evaluator:

    def __init__(self, evaluationDictionary):
        self.EVALUATIONDICTIONARY = evaluationDictionary
        self.__images = []
        self.__widths = []
        self.__heights = []
        self.__weights = []
        for name, weight in evaluationDictionary.items():
            temp = self.__loadGeotiff(name) #self.__loadImage(name)
            self.__images.append(temp[0])
            self.__widths.append(temp[1])
            self.__heights.append(temp[2])
            self.__weights.append(weight)

    def evaluate(self, x, y, bits, coords = False):
        eval = 0
        for i in range(len(self.__images)):
            if not coords:
                lx = int(self.chromosomeToCoord(bits, 0, self.__widths[i] - 1, x))
                ly = int(self.chromosomeToCoord(bits, 0, self.__heights[i] - 1, y))
                #print(i, x, y, len(self.__images[i]), lx + ly * self.__widths[i])
                eval += (self.__images[i][ly][lx] * self.__weights[i])
                #eval += (self.__images[i][lx + ly * self.__widths[i]] * self.__weights[i])
            else:
                eval += (self.__images[i][y][x] * self.__weights[i])
                #eval += (self.__images[i][int(x) + int(y) * int(self.__widths[i])] * self.__weights[i])
        return eval

    def __loadImage(self, file):
        im = Image.open(file)
        pixels = list(im.getdata())
        width, height = im.size
        return (pixels, width, height)

    def __loadGeotiff(self, file):
        gdal.UseExceptions()
        gt = gdal.Open(file)
        scores = gt.GetRasterBand(1).ReadAsArray()
        return [scores.tolist(), scores.shape[1], scores.shape[0]]

    def coordToChromosome(self, bits, min, max, coord):
        range = max - min
        coord -= min
        #num/range = x/precision, num * precision = range * x
        return (coord *((1 << bits) - 1)) / range

    def chromosomeToCoord(self, bits, min, max, chromosome):
        range = max - min
        return ((range * chromosome) / ((1 << bits)- 1)) + min

    def bruteForce(self, bits, width, height):
        smallest = len(self.__images[0])
        smallestIndex = 0
        for i in range(1, len(self.__images)):
            if len(self.__images[i]) < smallest:
                smallest = len(self.__images)
                smallestIndex = i
        bestScore = 0
        bestScoreIndices = []
        bestScoreDisplay = []
        for x in range(self.__widths[smallestIndex]):
            for y in range(self.__heights[smallestIndex]):
                score = self.evaluate(x, y, bits, True)
                if score > bestScore:
                    bestScore = score
                    bestScoreIndices = [[x, y]]
                elif score == bestScore:
                    bestScoreIndices.append([x, y])
        for coord in bestScoreIndices:
            x = width * (coord[0] / self.__widths[smallestIndex])
            y = height * (coord[1] / self.__heights[smallestIndex])
            bestScoreDisplay.append([x, y])
        return [bestScore, bestScoreIndices, bestScoreDisplay]






"""
Storage class holding a list of members to be acted upon by the GA
"""
class Population:

    """
    Creates a new population of a given size and randomly initializes its members
    :param numChromosomes The number of chromosomes per member
    :param numGenes The number of genes (bits) per chromosome
    :param popSize The size of (number of members in) the population
    """
    def __init__(self, numChromosomes, numGenes, popSize):
        self.__members = []
        for i in range(popSize):
            self.__members.append(Member(numChromosomes, numGenes, True))

    """
    Returns the size of the population
    :return The size of (or number of members in) the population
    """
    def getSize(self):
        return len(self.__members)

    """
    Returns a list containing all of the members within the populatio
    :return The list of the population's members
    """
    def getMembers(self):
        return self.__members

    """
    Exchanges the list of the population's member with the user-defined list
    :param newMembers The new list of members in the population
    """
    def setMembers(self, newMembers):
        self.__members = newMembers

    """
    Retrieves a reference to a member within the population
    :param index The index within the population of the member to be retrieved
    :return The member at the index
    """
    def getMember(self, index):
        return self.__members[index]

    """
    Exchanges member in the population with a new member
    :param index The index within the population of the member to be replaced
    :param newMember The member to replace the current member with
    """
    def setMember(self, index, newMember):
        self.__members[index] = newMember

    """
    Appends a member to the end of the population
    :param newMember The member to be added
    """
    def addMember(self, newMember):
        self.__members.append(newMember)

    """
    Removes a member from the population
    :param index The index of the member to be removed
    """
    def removeMember(self, index):
        self.__members.remove(index)


"""
Storage class representing each member of the population.  Holds gene (bit) values as decimal chromosome values
"""
class Member:

    """
    Create a new Member with a certain number of chromosomes and genes per chromosome and the option to initialize
    chromosomes randomly
    :param numChromosomes The number of chromosomes to be stored in the member
    :param numGenes The resolution of each chromosome (number of bits, e.g. 8 genes for a byte-sized chromosome)
    :param randInit False by default.  If true, each chromosome is initialized to a random number within its number of
        genes (i.e. 0 < x < 2 ^ (numGenes - 1)
    """
    def __init__(self, numChromosomes, numGenes, randInit = False):
        self.__chromosomes = []
        self.__NUMGENES = numGenes
        if randInit:
            for i in range(numChromosomes):
                self.__chromosomes.append(random.randrange((1 << numGenes) - 1))
        else:
            self.__chromosomes = [None] * numChromosomes

    """
    Gives the value of a particular chromosome of the member, or the aggregate of its genes in decimal
    :param index The index of the chromosome to be retrieved
    :return The chromosome at the given index
    """
    def getChromosome(self, index):
        return self.__chromosomes[index]

    """
    Used to set the value or genes (in decimal) of a particular chromosome of the member
    :param index The chromosome to be set
    :param value The new value (decimal) of the chromosome
    """
    def setChromosome(self, index, value):
        self.__chromosomes[index] = value

    """
    Retrieves the list of all of the chromosomes of the Member
    :return The list of the Member's chromosomes
    """
    def getChromosomes(self):
        return self.__chromosomes

    """
    Exchanges the Member's list of chromosomes with a new list of chromosomes
    :param newChromosomes The list of chromosomes to be used in place of the Member's current list
    """
    def setChromosomes(self, newChromosomes):
        self.__chromosomes = newChromosomes

"""
m1 = Member(2, 8, True)
m2 = Member(2, 8, True)
mask = random.randrange((1 << 16) - 1)
print(m1.getChromosomes())
print("m1",bin(m1.getChromosome(0)),bin(m1.getChromosome(1)))
print("m2",bin(m2.getChromosome(0)),bin(m2.getChromosome(1)))
print("mask",bin(mask))
mut = Mutator()
alg = Algorithm(8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, {})
n1 = alg.aggregateChromosomes(m1)
n2 = alg.aggregateChromosomes(m2)
print(bin(n1), bin(n2))
children = mut.pointCrossover(n1, n2, [1, 15], 16)
for c in children:
    print(bin(c))
mut.pointMutation(children[0], [0, 1, 2, 3, 4, 5, 6, 7, 8])
for c in children:
    print(bin(c))
c1 = alg.separateChromosomes(children[0])
c2 = alg.separateChromosomes(children[1])
print("c1",bin(c1[0]),bin(c1[1]))
print("c2",bin(c2[0]),bin(c2[1]))

x = alg.chromosomeToCoord(0, 255, 127)
print(x)
print(alg.coordToChromosome(0, 255, x))"""




