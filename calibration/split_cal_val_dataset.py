import csv
from collections import defaultdict
import os
import random

'''
script to split available data into calibration and validation dataset.
'''

script_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.dirname(script_path) + "/monica-data/projects/monica-germany/"
crop = "winter_wheat"#"silage_maize"# 
clustering_file = "individual_lks_WW.csv"#"individual_lks_SM.csv"#

cal_info = {
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

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def read_observations(calibration_target):
    observations = defaultdict(lambda: defaultdict(lambda: defaultdict(list))) #keys: lk, year, var_name; vals: val(s)
    cal_file = cal_info[calibration_target]
    with open(cal_file["f_name"]) as _:
        reader = csv.reader(_)
        reader.next()
        for row in reader:
            lk = int(row[cal_file["lk_col"]])
            if row[cal_file["val_col"]] == "na":
                #skip missing data
                continue
            val = float(row[cal_file["val_col"]])
            yr = int(row[cal_file["yr_col"]])
            if RepresentsInt(cal_file["var"]):
                #pheno file
                DWD_phase = float(row[cal_file["var"]])
                var = DWD_2_monica_stages[crop][DWD_phase]
            else:
                #yield file
                var = cal_file["var"]
            observations[lk][yr][var] = val
    return observations

#read general settings
lks = set()
with open(script_path + "/" + clustering_file) as _:
    reader = csv.reader(_)
    reader.next()
    for row in reader:
        lks.add(int(row[0]))

yield_data = read_observations("yield")
pheno_data = read_observations("phenology")
summary_avail_data_pheno = defaultdict(dict) #keys: lk, yr; vals: y/n
summary_avail_data_yield = defaultdict(dict) #keys: lk, yr; vals: y/n
years = range(1999, 2018)

for lk in lks:
    lk_yield_data = yield_data[lk]
    lk_pheno_data = pheno_data[lk]

    for yr in years:
        #pheno
        if yr in lk_pheno_data.keys():
            summary_avail_data_pheno[lk][yr] = "y"
        else:
            summary_avail_data_pheno[lk][yr] = "n"
        #yield
        if yr in lk_yield_data.keys():
            summary_avail_data_yield[lk][yr] = "y"
        else:
            summary_avail_data_yield[lk][yr] = "n"

with open(script_path + "/avail_data_yrs_" + crop + ".csv", "wb") as _:
    writer = csv.writer(_)
    header = ["lk", "year", "pheno", "yield"]
    writer.writerow(header)

    for lk in lks:
        for yr in years:
            row = []
            row.append(lk)
            row.append(yr)
            row.append(summary_avail_data_pheno[lk][yr])
            row.append(summary_avail_data_yield[lk][yr])
            writer.writerow(row)

lk_yy_yn_years = defaultdict(dict) #keys:lk, yrs; vals: yy/other years

with open(script_path + "/avail_data_summary_" + crop + ".csv", "wb") as _:
    writer = csv.writer(_)
    header1 = ["pheno-yield"]
    header2 = ["lk", "y-y", "y-n", "n-y", "n-n"]
    writer.writerow(header1)
    writer.writerow(header2)

    for lk in lks:
        row = [lk]
        y_y = 0
        y_n = 0
        n_y = 0
        n_n = 0
        years_yy = []
        other_years = []
        for yr in years:
            if summary_avail_data_pheno[lk][yr] == "y":
                if summary_avail_data_yield[lk][yr] == "y":
                    y_y += 1
                    years_yy.append(yr)
                else:
                    y_n += 1
                    other_years.append(yr)
            if summary_avail_data_pheno[lk][yr] == "n":
                if summary_avail_data_yield[lk][yr] == "y":
                    n_y += 1
                    other_years.append(yr)
                else:
                    n_n += 1
                    #nn years are not used at all for cal/val
        lk_yy_yn_years[lk]["yy"] = years_yy
        lk_yy_yn_years[lk]["other"] = other_years #yn or ny
        row.append(y_y)
        row.append(y_n)
        row.append(n_y)
        row.append(n_n)
        writer.writerow(row)

#split datasets based on rule
with open(script_path + "/lk2yrs_cal-val_" + crop + ".csv", "wb") as _:    
    writer = csv.writer(_)
    header = ["lk", "yr", "cal-val"]
    writer.writerow(header)
    for lk in lks:
        years_yy = lk_yy_yn_years[lk]["yy"]
        other_years = lk_yy_yn_years[lk]["other"]
        total_years = years_yy + other_years
        
        if len(years_yy) < len(other_years):
            #not enough years for calibration
            row = [lk, "skip", "skip"]
            writer.writerow(row)
            continue
        
        else:
            n_calib_years = len(total_years)/2
            
            #draw a random sample of years for calibration
            cal_years = random.sample(years_yy, n_calib_years)

            for yr in sorted(total_years):
                if yr in cal_years:
                    row = [lk, yr, "cal"]
                else:
                    row = [lk, yr, "val"]
                writer.writerow(row)

print "finished"