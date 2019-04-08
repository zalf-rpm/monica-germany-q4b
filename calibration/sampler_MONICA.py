from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import os
import spotpy
import spotpy_setup_MONICA
import csv
from datetime import date
import numpy as np
import json
from collections import defaultdict

script_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.dirname(script_path) + "/monica-data/projects/monica-germany/"
calib_id = "7"

with open(script_path + "/" + "calibration_log.csv") as _:
    reader =csv.reader(_)
    reader.next()
    for line in reader:
        if line[0] == calib_id:
            clustering_file = line[1]
            run_cells_file = line[2]
            calibration_target = line[3]
            params_to_calibrate = line[4]
            crop = line[5]
            dumped_envs_folder = script_path + line[6]
            rep = int(line[7])
            cal_val_years_file = line[8]

'''
clustering_file = "no_clustering_lks_SM.csv"#"no_clustering_lks_WW.csv"#
run_cells_file = "sample_setup_fullcover_SM.csv"#"sample_setup_fullcover_WW.csv"#
calibration_target = "phenology"#"sowing"#"yield"#  
params_to_calibrate = "calibratethese_pheno_SM.csv"#"calibratethese_sowing_SM.csv"#"calibratethese_sowing_WW.csv"# 
crop = "silage_maize"#"winter_wheat"# 
dumped_envs_folder = script_path + "/dumped_envs/localProducer-remoteMonica/3/"
'''

cal_info = {
    "sowing": {
        "f_name": data_path +"/lk_cleaned_pheno_" + crop + "_1999_2017.csv",
        "lk_col": 0,
        "yr_col": 1,
        "var": 3,
        "val_col": 4,
    },
    "phenology": {
        "f_name": data_path +"/lk_cleaned_pheno_" + crop + "_1999_2017.csv",
        "lk_col": 0,
        "yr_col": 1,
        "var": 3,
        "val_col": 4,
    },
    "yield": {
        "f_name": data_path + "/" + crop + "_yields_1999_2017.csv",
        "lk_col": 6,
        "yr_col": 8,
        "var": "yield",
        "val_col": 9,
    },
}
DWD_2_monica_stages = {
    "silage_maize": {
        10: 1,
        12: 2,
        67: 3,
        65: 4,
        5: 5,
        19: 6,
        20: 6,
        21: 6,
        24: 7
    },
    "winter_wheat": {
        10: 1,
        12: 2,
        15: 3,
        18: 3,
        19: 5,
        21: 5,
        24: 6  
    }    
}

def make_lambda(excel):
    return lambda v, p: eval(excel)

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

#read general settings
lk_clustering = {}
with open(script_path + "/" + clustering_file) as _:
    reader = csv.reader(_)
    reader.next()
    for row in reader:
        lk_clustering[int(row[0])] = int(row[1])

sim_map = defaultdict(lambda: defaultdict(list)) #keys: cluster, lk; vals: [(row, col)]
sim_lks = set()
with open(script_path + "/" + run_cells_file) as _:
    reader = csv.reader(_)
    reader.next()
    for row in reader:
        if row[4].upper() == "Y":
            lk = int(row[1])
            row_col = (int(row[2]), int(row[3]))
            cluster = lk_clustering[lk]
            sim_map[cluster][lk].append(row_col)
            sim_lks.add(lk)

#read calibration years per lk
lk2cal_years = defaultdict(list)
with open(script_path + "/" + cal_val_years_file) as _:
    reader = csv.reader(_)
    reader.next()
    for line in reader:
        if line[1] != "skip" and line[2] == "cal":
            lk = int(line[0])
            yr = int(line[1])
            lk2cal_years[lk].append(yr)

#read observations
observations = defaultdict(lambda: defaultdict(lambda: defaultdict(list))) #keys: lk, year, var_name; vals: val(s)
cal_file = cal_info[calibration_target]
with open(cal_file["f_name"]) as _:
    reader = csv.reader(_)
    reader.next()
    for row in reader:
        lk = int(row[cal_file["lk_col"]])
        if lk not in sim_lks:
            # no need to record all the data
            continue
        if row[cal_file["val_col"]] == "na":
            #skip missing data
            continue
        val = float(row[cal_file["val_col"]])
        yr = int(row[cal_file["yr_col"]])
        if yr not in lk2cal_years[lk]:
            #skip validation years (only cal is used here)
            continue
        if RepresentsInt(cal_file["var"]):
            #pheno file
            DWD_phase = float(row[cal_file["var"]])
            var = DWD_2_monica_stages[crop][DWD_phase]
        else:
            #yield file
            var = cal_file["var"]
        observations[lk][yr][var] = val

#read params to be calibrated
params = []
with open(script_path + "/" + params_to_calibrate) as paramscsv:
    dialect = csv.Sniffer().sniff(paramscsv.read(), delimiters=';,\t')
    paramscsv.seek(0)
    reader = csv.reader(paramscsv, dialect)
    next(reader, None)  # skip the header
    for row in reader:
        p={}
        p["name"] = row[0]
        p["array"] = row[1]
        p["low"] = float(row[2])
        p["high"] = float(row[3])
        p["stepsize"] = float(row[4])
        p["optguess"] = float(row[5])
        p["minbound"] = float(row[6])
        p["maxbound"] = float(row[7])
        if len(row) == 9 and row[8] != "":
            p["derive_function"] = make_lambda(row[8])
        params.append(p)

clu_counter = 0
total_clu = len(sim_map.keys())

for clu_id, lk_clu in sim_map.iteritems():
    clu_counter += 1
    print ("starting calibration of landkreise cluster: " + str(clu_id))
    print (str(clu_counter) + " out of " + str(total_clu))
    #pick only observations from the current cluster of lk
    clu_obs = defaultdict() #keys: lk, year, var_name; vals: val(s)
    for lk in lk_clu.keys():
        clu_obs[lk] = observations[lk]
    spot_setup = spotpy_setup_MONICA.spot_setup(params, sim_map[clu_id], clu_obs, calibration_target, dumped_envs_folder)
    #rep = 25

    #sampler = spotpy.algorithms.sceua(spot_setup, dbname='SCEUA', dbformat='ram')
    #sampler.sample(rep, ngs=len(params)+1, kstop=10)#, pcento=5)

    sampler = spotpy.algorithms.rope(spot_setup,dbname='ROPE',dbformat='ram')
    sampler.sample(rep)
    
    #sampler = spotpy.algorithms.mc(spot_setup,dbname='MC',dbformat='ram')
    #sampler = spotpy.algorithms.mle(spot_setup,dbname='MLE',dbformat='ram')
    #sampler = spotpy.algorithms.lhs(spot_setup,dbname='LHS',dbformat='ram')
    #sampler = spotpy.algorithms.sceua(spot_setup,dbname='SCEUA',dbformat='ram')
    #sampler = spotpy.algorithms.demcz(spot_setup,dbname='DE-MCz',dbformat='ram')
    #sampler = spotpy.algorithms.sa(spot_setup,dbname='SA',dbformat='ram')
    
    
    best_params = sampler.status.params

    results = sampler.getdata()    
    index, maximum = spotpy.analyser.get_maxlikeindex(results)
    bestmodelrun = list(spotpy.analyser.get_modelruns(results)[index][0])
    obs_list, obs_structure = spot_setup.evaluation(obs_structure=True)

    with open(script_path + '/opt_params/'+ calib_id + "_" + crop + "_"  + calibration_target + '_optimizedparams_cl_' + str(clu_id) + '.csv', 'wb') as _:
        writer = csv.writer(_)        
        for i in range(len(best_params)):
            outrow=[]
            arr_pos = ""
            if params[i]["array"].upper() != "FALSE":
                arr_pos = params[i]["array"]        
            outrow.append(params[i]["name"]+arr_pos)
            outrow.append(best_params[i])
            writer.writerow(outrow)
        if len(params) > len(best_params):
            reminder = []
            reminder.append("Don't forget to calculate and set derived params!")
            writer.writerow(reminder)
        text='optimized parameters saved!'
        print(text)
    
    with open(script_path + '/obs_vs_sim/' + calib_id + "_" + crop + "_" + calibration_target + '_obs_sim_cl_' + str(clu_id) + '.csv', 'wb') as _:
        writer = csv.writer(_)
        header = ["lk", "date", "obs", "sim"]
        writer.writerow(header)
        for i in range(len(obs_list)):
            outrow=[]
            outrow.append(obs_structure[i][0])
            outrow.append(obs_structure[i][1])
            outrow.append(obs_list[i])
            outrow.append(bestmodelrun[i])
            writer.writerow(outrow)
        text='obs_sim file saved!'
        print(text)

print("finished!")



