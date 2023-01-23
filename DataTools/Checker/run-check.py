#####################################################################################################
#
# Will look for data files within each replicate directory
#
# Command Line Inputs:
#
# data_directory: directory where data is located
#      selection: selection scheme used
#     experiment: experiment ran
#          model: island structure model
#        mig_int: migration interval
#       mig_bool: are migrations happening?
#       isl_size: island size
#        isl_cnt: number of islands
#         offset: seed offset
#
# Output: list of seeds that need to be reran in terminal display
#
# python3
#####################################################################################################

######################## IMPORTS ########################
import argparse
import pandas as pd
import sys
import os

# file location for data-params.py file
sys.path.insert(1, '../')
import data_params as dp

# return the number of rows in a csv file
def CountRows(file_name):
    # create pandas data frame of entire csv
    try:
        df = pd.read_csv(file_name)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()

    if(df.shape[0] == 0):
        return 0

    gens = df[dp.GENERATION].to_list()

    return gens[-1]

# responsible for looking through the data directories for success
def CheckDir(dir,sel,exp,mod,mig_int,mig_bool,isl_size,isl_cnt,offset):

    ISL_DIR = dp.GetIslandDirectory(exp,sel,mod,mig_int,isl_cnt,isl_size,mig_bool)

    # create seed data directories and check if exist
    DIR_DNE = []
    DAT_DNE = []
    DAT_DNF = []

    print('Full data Dir=', ISL_DIR + 'DIA_XXX__SEED_XXX/')
    print('Now checking data replicates sub directories')

    SEEDS = dp.SetSeedSets()
    for i in range(len(SEEDS)):
        for s in SEEDS[i]:
            seed = str(s + offset)
            DATA_DIR =  ISL_DIR + 'DIA_' + dp.SetDiagnostic(i) + '__SEED_' + seed + '/'

            print('Sub directory:', DATA_DIR)

            # add full directory to missing list if not there
            if os.path.isdir(DATA_DIR) == False:
                DIR_DNE.append(int(seed))
                continue

            # now check if the data file exists in full data director
            if os.path.isfile(DATA_DIR + 'data.csv') == False:
                DAT_DNE.append(int(seed))
                continue

            # make sure that the data.csv file did in fact finish all generations
            if CountRows(DATA_DIR + 'data.csv') != dp.GENERATIONS:
                DAT_DNF.append(int(seed))
                continue

    # print out the list of incomplete seeds
    print('Directories that were not created:', DIR_DNE)
    print('Data files that do not exist:', DAT_DNE)
    print('Data files that did not finish:', DAT_DNF)
    print('')
    print('Total list of unfinished seeds:')

    fin = DIR_DNE + DAT_DNF + DAT_DNE
    fin.sort()
    fins = ''
    for x in fin:
        fins = fins + str(x % dp.SEED_OFFSET) + ','
    # print out the sorted unfinished seed list
    print(fins[:len(fins)-1])

# runner
def main():
    # Generate and get the arguments
    parser = argparse.ArgumentParser(description="Data aggregation script.")
    parser.add_argument("data_directory", type=str,  help="Target experiment directory.")
    parser.add_argument("selection",      type=int,  help="Selection scheme we are looking for? \n0: Truncation\n1: Tournament\n2: Espilon Lexicase")
    parser.add_argument("experiment",     type=int,  help="Experiment we are assessing? \n0: Base experiments")
    parser.add_argument("model",          type=int,  help="Model we are checking")
    parser.add_argument("mig_int",        type=str,  help="Migration interval")
    parser.add_argument("mig_bool",       type=str,  help="Migration on or off?")
    parser.add_argument("isl_size",       type=str,  help="Island size")
    parser.add_argument("isl_cnt",        type=str,  help="Island count")
    parser.add_argument("offset",         type=int,  help="Experiment seed offset.")

    # Parse all the arguments
    args = parser.parse_args()
    data_dir = args.data_directory.strip()
    print('Data directory=',data_dir)
    selection = args.selection
    print('Selection scheme=', dp.SetSelection(selection))
    experiment = int(args.experiment)
    print('Experiment=', experiment)
    model = int(args.model)
    print('Experiment=', model)
    mig_int = args.mig_int
    print('Migration interval=', mig_int)
    mig_bool = args.mig_bool
    print('Migration bool=', mig_bool)
    isl_size = args.isl_size
    print('Island size=', isl_size)
    isl_cnt = args.isl_cnt
    print('Island count=', isl_cnt)
    offset = args.offset
    print('Diagnostic=', offset)

    # Get to work!
    print("\nChecking all related data directories now!")
    CheckDir(data_dir,selection,experiment,model,mig_int,mig_bool,isl_size,isl_cnt,offset)

if __name__ == "__main__":
    main()