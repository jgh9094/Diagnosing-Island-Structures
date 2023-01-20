#####################################################################################################
#####################################################################################################
# Will hold all of the data that we share across python scripts
#
# python3
#####################################################################################################
#####################################################################################################

import sys

# columns we are interested in grabbing
GENERATION = 'gen'
# pop level
POP_FIT_AVG = 'pop_fit_avg'
POP_FIT_MAX = 'pop_fit_max'
ISL_FIT_AVG = 'isl_fit_avg'
# coverage data
POP_SAT_COV = 'pop_sat_cov'
ISL_SAT_COV = 'isl_sat_cov'
POP_ACT_COV = 'pop_act_cov'
ISL_ACT_COV = 'isl_act_cov'
# selection scheme data
LOS_DIV = 'los_div'
SEL_PRE = 'sel_pre'

# seed experiements replicates range
REPLICATES = 100
GENERATIONS = 50000
SEED_OFFSET = 1000
RESOLUTION = 100

# return appropiate string dir name (based off run.sb file naming system)
def SetSelection(s):
    # case by case
    if s == 0:
        return 'TRUNCATION'
    elif s == 1:
        return 'TOURNAMENT'
    elif s == 2:
        return 'LEXICASE'
    else:
        sys.exit("UNKNOWN SELECTION")

# return appropiate string dir name (based off run.sb file naming system)
def SetDiagnostic(s):
    # case by case
    if s == 0:
        return 'EXPLOITATION_RATE'
    elif s == 1:
        return 'ORDERED_EXPLOITATION'
    elif s == 2:
        return 'CONTRADICTORY_OBJECTIVES'
    elif s == 3:
        return 'MULTIPATH_EXPLORATION'
    else:
        sys.exit('UNKNOWN DIAGNOSTIC')

# return the correct amount of seed ran by experiment treatment
def SetSeedSets():
    seed = []
    seed.append([x for x in range(1,101)])
    seed.append([x for x in range(101,201)])
    seed.append([x for x in range(201,301)])
    seed.append([x for x in range(301,401)])
    return seed

# set the appropiate list of variables we are checking for
def SetModelType(m):
    if m == 0:
        return 'EA'
    elif m == 1:
        return 'IS'
    elif m == 2:
        return 'NMIS'
    else:
        sys.exit("UNKNOWN VARIABLE LIST")

# set the appropiate list of variables we are checking for
def SetExperimentType(e):
    if e == 0:
        return 'BASE-EXPERIMENTS'
    else:
        sys.exit("UNKNOWN VARIABLE LIST")