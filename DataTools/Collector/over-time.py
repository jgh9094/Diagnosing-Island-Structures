#####################################################################################################
#
# Will get value for each data type in data_params each 100 generations (gen % 100 == 0)
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
# dump_directory: where are we dumping everything
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

# responsible for looking through the data directories for success
def Directories(dir,sel,exp,mod,mig_int,mig_bool,isl_size,isl_cnt,offset,dump):

    # get island directory
    ISL_DIR = dp.GetIslandDirectory(dir,exp,sel,mod,mig_int,isl_cnt,isl_size,mig_bool)

    print('Full data Dir=', ISL_DIR + 'DIA_XXX__SEED_XXX/')
    print('Now checking data replicates sub directories')

    SEEDS = dp.SetSeedSets()
    # gens we are expecting
    GEN_LIST = [x for x in range(int(dp.GENERATIONS)+1) if x % dp.RESOLUTION == 0]
    # collect all data
    DF_LIST = []
    for i in range(len(SEEDS)):

        for s in SEEDS[i]:
            seed = str(s + offset)
            DATA_DIR =  ISL_DIR + 'DIA_' + dp.SetDiagnostic(i) + '__SEED_' + seed + '/'

            print('Sub directory:', DATA_DIR)

            DATA_DIR += '/data.csv'

            # check if data file even exists
            if os.path.isfile(DATA_DIR) == False:
                sys.exit('DNE: ' + DATA_DIR)

            # create pandas data frame of entire csv and grab the row
            df = pd.read_csv(DATA_DIR)
            df = df.iloc[::dp.RESOLUTION, :]

            # time to export the data
            cdf = pd.DataFrame(
                {   'Generations':          pd.Series(GEN_LIST),
                    'Structure':            pd.Series([dp.SetModelType(mod)] * len(GEN_LIST)),
                    'Selection Scheme':     pd.Series([dp.SetSelection(sel)] * len(GEN_LIST)),
                    'Diagnostic':           pd.Series([dp.SetDiagnostic(i)] * len(GEN_LIST)),
                    dp.POP_FIT_AVG:         pd.Series(df[dp.POP_FIT_AVG].tolist()),
                    dp.POP_FIT_MAX:         pd.Series(df[dp.POP_FIT_MAX].tolist()),
                    dp.ISL_FIT_AVG:         pd.Series(df[dp.ISL_FIT_AVG].tolist()),
                    dp.POP_SAT_COV:         pd.Series(df[dp.POP_SAT_COV].tolist()),
                    dp.ISL_SAT_COV:         pd.Series(df[dp.ISL_SAT_COV].tolist()),
                    dp.POP_ACT_COV:         pd.Series(df[dp.POP_ACT_COV].tolist()),
                    dp.ISL_ACT_COV:         pd.Series(df[dp.ISL_ACT_COV].tolist()),
                    dp.LOS_DIV:             pd.Series(df[dp.LOS_DIV].tolist()),
                    dp.SEL_PRE:             pd.Series(df[dp.SEL_PRE].tolist())
                })
            DF_LIST.append(cdf)

    fin_df = pd.concat(DF_LIST)

    fin_df.to_csv(path_or_buf= dump + 'over-time-' + dp.SetSelection(sel)  + '.csv', index=False)

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
    parser.add_argument("dump_directory", type=str,  help="Target experiment directory.")

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
    dump_directory = args.dump_directory
    print('Dump directory=', dump_directory)

    # Get to work!
    print("\nChecking all related data directories now!")
    Directories(data_dir,selection,experiment,model,mig_int,mig_bool,isl_size,isl_cnt,offset,dump_directory)

if __name__ == "__main__":
    main()