#####################################################################################################
#####################################################################################################
# Will hold all of the data that we share across python scripts
#
# python3
#####################################################################################################
#####################################################################################################

import sys
import os

# columns we are interested in grabbing
GENERATION = 'gen'
# pop level
POP_FIT_AVG = 'pop_fit_avg'
POP_FIT_MAX = 'pop_fit_max'
POP_OPT_MAX = 'pop_opt_max'
ISL_FIT_AVG = 'isl_fit_avg'
# coverage data
POP_SAT_COV = 'pop_sat_cov'
ISL_SAT_COV = 'isl_sat_cov'
POP_ACT_COV = 'pop_act_cov'
ISL_ACT_COV = 'isl_act_cov'
# selection scheme data
LOS_DIV = 'los_div'
SEL_PRE = 'sel_pre'

DATA_LIST = [POP_FIT_AVG,POP_FIT_MAX,ISL_FIT_AVG,POP_SAT_COV,ISL_SAT_COV,POP_ACT_COV,ISL_ACT_COV,LOS_DIV,SEL_PRE]

# seed experiements replicates range
REPLICATES = 100
GENERATIONS = 50000
SEED_OFFSET = 1000
RESOLUTION = 100
DIMENTIONALITY = 100

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


def GetIslandDirectory(dir,exp,sel,mod,mig_int,isl_cnt,isl_size,mig_bool):
    # check if data dir exists
    if os.path.isdir(dir):
        print('Data dirctory exists=', dir)
    else:
        sys.exit('DATA DIRECTORY DOES NOT EXIST: ' + dir)

    # check if experiment folder exists
    EXP_DIR = dir + SetExperimentType(exp) + '/'
    if os.path.isdir(EXP_DIR):
        print('Experiment data folder exists', EXP_DIR)
    else:
        sys.exit('EXPERIMENT DATA DIRECTORY DOES NOT EXIST: ' + EXP_DIR)

    # check if selection scheme folder exists
    SEL_DIR = EXP_DIR + SetSelection(sel) + '/'
    if os.path.isdir(SEL_DIR):
        print('Selection scheme data folder exists', SEL_DIR)
    else:
        sys.exit('SELECTION DATA DIRECTORY DOES NOT EXIST: ' + SEL_DIR)

    # check if model folder exists
    MOD_DIR = SEL_DIR + SetModelType(mod) + '/'
    if os.path.isdir(MOD_DIR):
        print('ModeL data folder exists', MOD_DIR)
    else:
        sys.exit('MODEL DATA DIRECTORY DOES NOT EXIST: ' + MOD_DIR)

    # check if island configuration folder exists
    ISL_DIR = MOD_DIR + 'INT_' + mig_int + '__CNT_' + isl_cnt + '__SIZE_' + isl_size + '__MIG_' + mig_bool + '/'
    if os.path.isdir(ISL_DIR):
        print('Island configuration folder exists', ISL_DIR)
    else:
        sys.exit('ISLAND CONFIGURATION FOLDER DOES NOT EXISTS: ' + ISL_DIR)

    return ISL_DIR