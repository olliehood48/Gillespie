import random
import math

class Species:
    def __init__(self, name, initialPop):
        self.name = name
        self.initialPop = initialPop
        self.currentPop = initialPop

class TwoSpeciesPayoffMatrix:
    def __init__(self, species1, species2, CC, CD, DD, DC):
        self.species1 = species1
        self.species2 = species2
        self.CC = CC
        self.CD = CD
        self.DD = DD
        self.DC = DC
        self.matrix = [CC, CD, DD, DC]

    def describe(self):
        print(f"{self.species1.name}: current population - {self.species1.currentPop}")
        print(f"{self.species2.name}: current population - {self.species2.currentPop}")
        print(f"matrix: {self.matrix}")

stags = Species(name = 'stag', initialPop = 150)
hares = Species(name = 'hare', initialPop = 50)
lions = Species(name = 'lion', initialPop = 100)
species = [stags, hares, lions]
stagHareMatrix = TwoSpeciesPayoffMatrix(stags, hares, 4, 1, 3, 2)
stagLionMatrix = TwoSpeciesPayoffMatrix(stags, lions, 4, 1, 6, 0)
harelionMatrix = TwoSpeciesPayoffMatrix(hares, lions, 3, 2, 6, 0)
matrices = {
    ('stag', 'hare'): stagHareMatrix,
    ('stag', 'lion'): stagLionMatrix,
    ('hare', 'lion'): harelionMatrix,
}

## function to return birth and death rates from two player game based on matrix and populations
def retBirthAndDeathRates(specificMatrix, i, N, s):
    a, b, c, d = specificMatrix
    payOffC = a*( (i - 1) / (N - 1) ) + b*( (N - i)/(N - 1) )
    payOffD = c*( i/ (N - 1) ) + d*( (N - i - 1)/(N - 1) )
    fitnessC = 1 - s + s*payOffC
    fitnessD = 1 - s + s*payOffD
    fbar = (i/N)*fitnessC + (1-(i/N))*fitnessD
    birthRateC = (i/N) * ((N - i) / (N)) * (fitnessC/fbar)
    birthRateD = (i/N) * ((N - i) / (N)) * (fitnessD/fbar) ## NOTE - in  birth/death, birthRateD = deathRateC 
    return birthRateC, birthRateD

## function to choose N numbers from an array at random - used to determine which species
def chooseNAtRandomNoRepeats(probs, N):
    ## make cumulative probs to compare for which is chosen
    cumProbs = []
    for prob in range(len(probs) - 1):
        cumProbs.append(sum(probs[:prob + 1]))
    ## chosen will contain the chosen probs' positions
    chosen = []
    ## loop while chosen is less than N
    while len(chosen) < N:
        ## generate randInt to compare among cumProbs
        randInt = random.uniform(0, 1)
        ## determine where randInt lies among cumProbs
        for place, prob in enumerate(cumProbs):
            if randInt < prob:
                if place not in chosen:
                    chosen.append(place)
                    break
        if (len(cumProbs)) not in chosen:
            chosen.append(len(cumProbs))
    return chosen

## returns two species
def determine2Species(species):
    speciesPopSorted = sorted(species, key = lambda specie: specie.currentPop)
    currentPops = [specie.currentPop for specie in speciesPopSorted]
    totalPop = sum(currentPops)
    currentProps = [currentPop/totalPop for currentPop in currentPops]
    twoSpecies = chooseNAtRandomNoRepeats(currentProps, 2)
    return speciesPopSorted[twoSpecies[0]], speciesPopSorted[twoSpecies[1]]

## algorithm returns updated population based off how 
def gillAlg(s):
    ## determine which two individuals are interacting based off population proportions
    twoSpecies = determine2Species(species)
    species1 = twoSpecies[0]
    species2 = twoSpecies[1]
    if (species1.currentPop == 0) or (species2.currentPop == 0):
        print(f"population 1 or 2 is wiped out")
        return
    ## determine relevant matrix from matrices defined above
    for matrixKey in list(matrices.keys()):
        if (twoSpecies[0].name in matrixKey) and (twoSpecies[1].name in matrixKey):
            if matrixKey[0] != species1.name:
                buffer = species1
                species1 = species2
                species2 = buffer
            matrix = matrices[matrixKey].matrix
    ## determine i, N
    species1Pop = species1.currentPop
    species2Pop = species2.currentPop
    totalCombined = species1Pop + species2Pop
    ## random numbers for waiting time incrament and reference for event stastics (determining birth/death)
    r1 = random.uniform(0, 1)
    r2 = random.uniform(0, 1)
    ## determine propensity function
    birthRateC, deathRateC = retBirthAndDeathRates(matrix, species1Pop, totalCombined, s)
    propFunc = birthRateC + deathRateC
    ## find time increment for event to happen
    timeInc = -(1/propFunc)*math.log(r1) ## think it might be: math.exp( -(propFunc) * -(1/propFunc)*math.log(r1) ) 
    ## determine birth or death from event
    # find stastic to compare if C has birth or death.
    eventStat = birthRateC/(birthRateC + deathRateC)
    if r2 <= eventStat:
        twoSpecies[0].currentPop += 1
        twoSpecies[1].currentPop -= 1
        # print(f"Birth of C, in time incrament {timeInc}")
        return 
    twoSpecies[0].currentPop -= 1
    twoSpecies[1].currentPop += 1
    # print(f"Death of C, in time incrament {timeInc}")
    return 

for i in range(100):
    gillAlg(1)
    for specie in species:
        print(f"{specie.name}: CP = {specie.currentPop}")
    print()