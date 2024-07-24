import random
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

################ FUNCTIONS - variable parameters (initial conditions) below ################

# -calculates revenue at a given time
def returnAdRevenue(i, a):
    # take values k, a0 defined at bottomk*
    revenue = i*k*((1/a) - (1/a0))
    return revenue

# -returns birth rate and death rate of C
def retBirthAndDeathRates(matrix, i, N, s):
    a, b, c, d = matrix
    payOffC = a*( (i - 1) / (N - 1) ) + b*( (N - i)/(N - 1) )
    payOffD = c*( i/ (N - 1) ) + d*( (N - i - 1)/(N - 1) )
    fitnessC = 1 - s + s*payOffC
    fitnessD = 1 - s + s*payOffD
    fbar = (i/N)*fitnessC + (1-(i/N))*fitnessD
    birthRateC = (i/N) * ((N - i) / (N)) * (fitnessC/fbar)
    birthRateD = (i/N) * ((N - i) / (N)) * (fitnessD/fbar) 
    return birthRateC, birthRateD

# -performs one iteration of the Gillespie algorithm
# -updates i relative to C birth or death
def gillAlg(i, N, s, matrix):
    r1 = random.uniform(0, 1)
    r2 = random.uniform(0, 1)
    birthRateC, deathRateC = retBirthAndDeathRates(matrix, i, N, s)
    propFunc = birthRateC + deathRateC
    ## find time increment for event to happen
    timeInc =  -(1/propFunc)*math.log(r1) 
    ## determine birth or death from event
    # find stastic to compare if C has birth or death.
    eventStat = birthRateC/(birthRateC + deathRateC)
    if r2 <= eventStat:
        i += 1
        # print("Birth of C")
        return i, N, timeInc
    i -= 1
    # print("Death of C")
    return i, N, timeInc

# -performs full simulation of Gillespie algorithm, until i = 0 or i = N
# -returns by appending the 'fixTime' (fixation time)
# -'fixTime' is sum of 'timeInc's (time incraments) per birth/death, calculated in gillAlg()
def gillSim(matrix, Cfrac, N, results, s, countSoFar, totalCount):
    i = int(math.floor(N*Cfrac))  ## math.floor just in case not whole number
    # print(f"TEST: performing test for {str(int(i))} Cs in total of {str(int(N))} ({Cfrac*100}%), at selection strength {str(s)}")
    fixTime = 0
    cumFixTime = 0
    revenue = 0
    while not ((i == 0) or (i == N)):
        # print(f"TEST: i: {i} -- N: {N} -- revenue: {revenue}")
        i, N, timeInc = gillAlg(i, N, s, matrix)
        fixTime += timeInc
        cumFixTime += timeInc
        if cumFixTime > 10:
            revenue += returnAdRevenue(i, matrix[0])
            cumFixTime = 0
        
        # time.sleep(0.01)
    if i == 0:
        # print(f"TEST: D has invaded, time taken = {fixTime}s")
        # print()
        results[str(s)][str(matrix[0])]['D']['fixTimes'].append(fixTime)
        results[str(s)][str(matrix[0])]['D']['revenues'].append(revenue)
        print(f"{countSoFar}/{totalCount} - {str(math.floor((countSoFar/totalCount)*100))}% - fixated at D")
        return results
    # print(f"TEST: C has invaded, time taken = {fixTime}s")
    # print()
    results[str(s)][str(matrix[0])]['C']['fixTimes'].append(fixTime)
    results[str(s)][str(matrix[0])]['C']['revenues'].append(revenue)
    print(f"{countSoFar}/{totalCount} - {str(math.floor((countSoFar/totalCount)*100))}% - fixated at C")
    return results

# -iterates gillSim according to iterations arg
# -returns results in form {N: {cFrac: {s: {iteration: []....
def runGS(iterations, Cfraction, N, Ss, aValues):
    countSoFar = 0
    totalCount = iterations * len(Ss) * len(aValues)
    ## make new results dict for each iteration amount to be added to fullResults
    results = {}
    for s in Ss:
        for a in aValues:
           results.setdefault(str(s), {}).setdefault(str(a), {'C': {'fixTimes' : [], 'revenues': []}, 'D': {'fixTimes' : [], 'revenues': []}})
    for s in Ss:  
        for a in aValues:
            b = a/1.1
            d = 11
            c = d/1.1
            matrix = [a, b, c, d]   
            for iteration in range(iterations):
                countSoFar+=1
                results = gillSim(matrix, Cfraction, N, results, s, countSoFar, totalCount)
    return results

# -display results from fullResults dictionary
def displayResults(fullResults, iterations, Cfraction, N, Ss, aValues):
    print()
    print()
    print("         **** FULL RESULTS ****")
    for s in Ss:
        for a in aValues:
            cFixValues = fullResults[str(s)][str(a)]['C']['fixTimes']
            cRevValues = fullResults[str(s)][str(a)]['C']['revenues']
            dFixValues = fullResults[str(s)][str(a)]['D']['fixTimes']
            dRevValues = fullResults[str(s)][str(a)]['D']['revenues']
            print(f"s = {s}, a = {a}, {iterations} iterations:")
            print(f"    C fixation probability: {len(cFixValues)/(len(cFixValues) + len(dFixValues))}")
            print(f"    average fixation time: {(sum(cFixValues) + sum(dFixValues))/(len(cFixValues) + len(dFixValues))}")
            print(f"    average remaining revenue: {(sum(cRevValues) + sum(dRevValues))/(len(cRevValues) + len(dRevValues))}")
            print()
    
# -get results and display them
def runGSwithResultsOutput(iterations, Cfraction, N, Ss, aValues):
    ## get full results dictionary
    fullResults = runGS(iterations, Cfraction, N, Ss, aValues)
    ## display the results
    displayResults(fullResults, iterations, Cfraction, N, Ss, aValues)

    # print(f"TEST: {fullResults}")
    return fullResults

# -plots graphs for mean fixation time from fullResults
def plotData(dataset):
    num_rows = (len(Ss) + 1) // 2
    num_columns = 2
    num_subplots = len(Ss)

    fig, axs = plt.subplots(num_rows, num_columns, figsize=(10, 7))

    for i, s in enumerate(Ss):
        row = i // num_columns
        col = i % num_columns
        
        cFixProbs = []
        avgFixTimes = []  # Initialize list for average fixation times
        for a in aValues:
            cValues = dataset[str(s)][str(a)]['C']['fixTimes']
            dValues = dataset[str(s)][str(a)]['D']['fixTimes']
            cFixProb = len(cValues) / (len(cValues) + len(dValues))
            cFixProbs.append(cFixProb)
            avgFixTime = sum(cValues) / len(cValues)
            avgFixTimes.append(avgFixTime)  # Append average fixation time
        
        # Plot data on the appropriate subplot
        if num_rows == 1:  # If there's only one row of subplots
            axs[col].plot(aValues, cFixProbs, label=f'Iterations: {iterations}')
            axs[col].set_title(f"C Fixation Probabilities and Time with s = {s}")
            axs[col].set_xlabel("a Values")
            axs[col].set_ylabel("Fixation Probability of C")
            axs[col].legend()
            
            # Create a secondary y-axis for average fixation time
            axs2 = axs[col].twinx()
            axs2.plot(aValues, avgFixTimes, color='red', linestyle='-', label='Average Fixation Time')
            axs2.set_ylabel('Average Fixation Time', color='red')
        else:  # If there are multiple rows of subplots
            axs[row, col].plot(aValues, cFixProbs, label=f'Iterations: {iterations}')
            axs[row, col].set_title(f"C Fixation Probabilities and Time with s = {s}")
            axs[row, col].set_xlabel("a Values")
            axs[row, col].set_ylabel("Fixation Probability of C")
            axs[row, col].legend()
            
            # Create a secondary y-axis for average fixation time
            axs2 = axs[row, col].twinx()
            axs2.plot(aValues, avgFixTimes, color='red', linestyle='-', label='Average Fixation Time')
            axs2.set_ylabel('Average Fixation Time', color='red')

    # Hide the empty subplot if there's an odd number of subplots
    if num_subplots % 2 == 1:
        if num_subplots == 1:
            fig.delaxes(axs)
        else:
            fig.delaxes(axs[-1, -1])

    plt.tight_layout()
    plt.show()

    plotAvalues(dataset)

# -plot the graphs for values of a with fixProb = 1
def plotAvalues(dataset):
    refinedData = {}
    ## add simulations where C had fixation probability 1 to the refined data
    for s in Ss:
        for a in aValues:
            cFixValues = dataset[str(s)][str(a)]['C']['fixTimes']
            dFixValues = dataset[str(s)][str(a)]['D']['fixTimes']
            if len(dFixValues) == 0:
                refinedData.setdefault(str(s), {}).setdefault(str(a), {'fixValues': [], 'revenues': []})
                refinedData[str(s)][str(a)]['fixTimes'] = cFixValues
                refinedData[str(s)][str(a)]['revenues'] = dataset[str(s)][str(a)]['C']['revenues']

    sVals = list(refinedData.keys())

    num_rows = (len(sVals) + 1) // 2
    num_columns = 2

    fig, axs = plt.subplots(num_rows, num_columns, figsize=(10, 7))

    for i, s in enumerate(sVals):
        row = i // num_columns
        col = i % num_columns

        aVals = list(refinedData[str(s)].keys())
        aVals.sort()
        avgFixTimes = []
        avgRevRems = []
        for a in aVals:
            fixTimes = refinedData[str(s)][str(a)]['fixTimes']
            avgFixTime = sum(fixTimes) / len(fixTimes)
            avgFixTimes.append(avgFixTime)
            revRems = refinedData[str(s)][str(a)]['revenues']
            avgRevRem = sum(revRems) / len(revRems)
            avgRevRems.append(avgRevRem)
        axs[row, col].plot(aVals, avgFixTimes, label=f"s: {s}")

        axs[row, col].set_title(f"Average Fixation Times for s: {s}")
        axs[row, col].set_xlabel("a Values")
        axs[row, col].set_ylabel("Average Fixation Times")
        axs[row, col].legend()

        # Create a secondary y-axis for average fixation time
        axs2 = axs[row, col].twinx()
        axs2.plot(aVals, avgRevRems, color='olive', linestyle='-', label='Average Total Revenue')
        axs2.set_ylabel('Average Total Revenue', color='red')
        axs2.axhline(0, color='red', linestyle='--')  # Horizontal line at y = 0
        axs2.legend()

    # Hide any empty subplots
    for i in range(len(sVals), num_rows * num_columns):
        axs.flatten()[i].axis('off')

    plt.tight_layout()
    plt.show()

    return



################ Below is changeable/conditions ################

### only need to change these variables to run script

## test values
iterations = 100
N = 1000 ## change to 1000 for official
Cfraction = 0.8
## s = 0 implies fitnesses are equal so constant fixation probability == 0.5
Ss = [0.01, 0.02, 0.03, 0.04]
## adjust aValues to give different payoff matrices accordingly
aValues = [num for num in range(11, 18)] 
k = 50
a0 = 11






################ Actual running code ################

# GS with results
fullResults = runGSwithResultsOutput(iterations = iterations, Cfraction = Cfraction, N = N, Ss = Ss, aValues = aValues) 
# # plotting sample mean fixation time and probability
plotData(fullResults)

