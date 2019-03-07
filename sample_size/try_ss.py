import csv
import os
from collections import defaultdict
import random
import numpy as np

#0. Settings
sim_id = 1
target_var = "Anthesis-doy"#"Yield-final"
n_random_samples = 1
step_ss = 20
min_ss = 20
max_ss = 400
out_file = "ss_"+ target_var + "_SM.csv"
max_row = 320

script_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.dirname(script_path) + "/csv-out/" + str(sim_id)

data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path + "/", f))]

lk_cell_year_var = defaultdict(lambda:(defaultdict(dict)))
data_fields = {}

#1. read in data
for f in data_files:
    row_id = int(f.split(".")[0].split("-")[1])
    if row_id > max_row:
        continue
    with open(data_path + "/" + f) as _:
        print "reading file " + str(f)
        reader = csv.reader(_)
        header = reader.next()
        for i in range(len(header)):
            data_fields[header[i]] = i
        for row in reader:
            lk = row[data_fields["lk"]]
            cell = row[data_fields["row"]] + "_" + row[data_fields["col"]]
            year = int(row[data_fields["Harvest-year"]])
            try:
                my_var = float(row[data_fields[target_var]])
            except:
                #"NA"
                print "skipping row " + str(row) + " year " + str(year)
                continue
            lk_cell_year_var[lk][cell][year] = my_var

#2. sample data
print "sampling..."
all_lks = set(lk_cell_year_var.keys())
sample_sizes = range(min_ss, max_ss+1, step_ss)

#2.1. calculate lk value of target var
lk_year_val = defaultdict(dict)
for lk in all_lks:
    year_vals = defaultdict(list)
    for cell in lk_cell_year_var[lk].keys():
        for yr in lk_cell_year_var[lk][cell].keys():
            year_vals[yr].append(lk_cell_year_var[lk][cell][yr])
    for yr in year_vals.keys():
        lk_year_val[lk][yr] = np.average(year_vals[yr])

#2.2. sample
out_data = []
count_lk = 0
for lk in all_lks:
    count_lk += 1    
    print "lk " + str(count_lk) + " out of " + str(len(all_lks))
    n_cells = len(lk_cell_year_var[lk].keys())
    for ss in sample_sizes:
        for rd_s in range(0, n_random_samples):
            sampled_cells = random.sample(lk_cell_year_var[lk].keys(), min(ss, n_cells))
            year_vals = defaultdict(list)
            for cell in sampled_cells:
                for yr in lk_cell_year_var[lk][cell].keys():
                    year_vals[yr].append(lk_cell_year_var[lk][cell][yr])
            #store output
            for yr in year_vals.keys():
                perc_coverage = round(float(len(sampled_cells))/float(n_cells) * 100, 0)
                sample_val = np.average(year_vals[yr])
                abs_delta = abs(lk_year_val[lk][yr] - sample_val)
                my_out = {
                    "lk": lk,
                    "n_cells": n_cells,
                    "sample_id": rd_s,
                    "sample_size": ss,
                    "year": yr,
                    "lk_val": lk_year_val[lk][yr],
                    "sample_val": sample_val,
                    "rounded_perc_coverage": perc_coverage,
                    "abs_delta": abs_delta
                }
                out_data.append(my_out)

#write data
with open(script_path + "/" + out_file, "wb") as _:
    writer = csv.writer(_)
    header = list(out_data[0].keys())
    writer.writerow(header)
    for out in out_data:
        row = []
        for field in out.keys():
            row.append(out[field])
        writer.writerow(row)

print "finished!"
            
            


