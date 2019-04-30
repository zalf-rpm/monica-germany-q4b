import os
from os import listdir
from os.path import isfile, join, dirname, abspath
from collections import defaultdict
import csv
import pandas as pd
import numpy as np

# User settings
crop = "silage_maize"# "winter_wheat"#
exps = {
    "start": 7,
    "end": 7,
    "run-id-converter": -6
}

script_path = dirname(abspath(__file__))
sim_data_path = dirname(script_path) + "/csv-out/"
project_data_path = dirname(script_path) + "/monica-data/projects/monica-germany/"

obs_file = crop + "_yields_1999_2017.csv"
crunch_out_file = crop + "_crunched.csv"

# 1. split cal/val dataset
lk2cal_val = defaultdict(dict) #keys: lk, yr
cal_val_years_file = "lk2yrs_cal-val_" + crop + ".csv"
with open(project_data_path + cal_val_years_file) as _:
    print("reading " + cal_val_years_file)
    reader = csv.reader(_)
    reader.next()
    for line in reader:
        if line[1] != "skip":
            lk = int(line[0])
            yr = int(line[1])
            lk2cal_val[lk][yr] = line[2]

# 2. read observations
observations = defaultdict(dict) #keys: lk, year 
with open(project_data_path + obs_file) as _:
    print("reading " + obs_file)
    reader = csv.reader(_)
    reader.next()
    for row in reader:
        if row[9] == "na":
            #skip missing data
            continue
        lk = int(row[6])
        yr = int(row[8])
        val = float(row[9])
        observations[lk][yr] = val

with open(crunch_out_file, "wb") as _:
    writer = csv.writer(_)
    header = ["run-id", "exp-id", "lk", "year", "cal-val", "obs", "sim"]
    writer.writerow(header)

for run_id in range(exps["start"], exps["end"] + 1):
    #if run_id != 40:
    #    continue
    exp_id = run_id + exps["run-id-converter"]
    # 3. read output    
    exp_dframes = []
    exp_path = sim_data_path + str(run_id) + "/"
    out_files = [join(exp_path, f) for f in listdir(exp_path) if isfile(join(exp_path, f))]
    
    for f in out_files:
        #test:
        #if "row-5.csv" not in f:
        #    continue
        print("loading " + f)    
        exp_dframes.append(pd.read_csv(f))
    exp_data = pd.concat(exp_dframes)

    lks = set(lk2cal_val.keys())
    for lk in lks:
        lk_data = []    
        #test
        #if lk != 1054:
        #    continue
        print("crunching data of lk " + str(lk))
        yrs = set(lk2cal_val[lk].keys())
        for yr in yrs:
            cal_or_val = lk2cal_val[lk][yr]
            try:
                obs_yield = observations[lk][yr]
            except:
                if cal_or_val == "val":
                    # validation years contain years with:
                    # a) pheno and yield data
                    # b) pheno or yield data
                    # in case of b, no need to worry about missing data
                    continue
                else:
                    print("missing data in calibration year! Check lk: " + str(lk) + " yr: " + str(year))
                    exit()
            sim_data = exp_data.loc[(exp_data["lk"] == lk) &
                                    (exp_data["Harvest-year"] == yr), "Yield-final"].values
            sim_yield = np.average(sim_data)            
            row = [run_id, exp_id, lk, yr, cal_or_val, obs_yield, sim_yield]
            lk_data.append(row)
    
        # 4. write output
        with open(crunch_out_file, "ab") as _:
            writer = csv.writer(_)
            for row in lk_data:
                writer.writerow(row)

print("finished!")

