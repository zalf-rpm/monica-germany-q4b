import csv
import os
from collections import defaultdict
import random
import numpy as np

'''
Define the sample of cells for dumping the envs (from the producer) that will be used for calibration
'''

#0. Settings
sim_id = 4
out_file = "sample_setup_"+ str(sim_id) + ".csv"
relative_sample_size = 0.05 #determined after running try_ss.py in "sample_size" folder
max_row = 9999
###

script_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.dirname(script_path) + "/csv-out/" + str(sim_id)

data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path + "/", f))]

lk_cells = defaultdict(set)

#1. read the data
data_fields = {}
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
            row_col = (row[data_fields["row"]], row[data_fields["col"]])
            lk_cells[lk].add(row_col)

#2. sample row_col for later use and write out file
total_cells = 0
total_sampled = 0
with open (script_path + "/" + out_file, "wb") as _:
    writer = csv.writer(_)
    header = ["sim_ID", "lk", "row", "col"]
    writer.writerow(header)
    for lk in lk_cells.keys():
        n_cells = len(lk_cells[lk])
        print "sampling " + str(n_cells) + " cells from lk " + str(lk)
        sample_size = int(n_cells * relative_sample_size + 1)
        sampled_cells = random.sample(lk_cells[lk], sample_size)        
        for cell in sampled_cells:
            row = [sim_id, lk, cell[0], cell[1]]
            writer.writerow(row)
        total_cells += n_cells
        total_sampled += len(sampled_cells)

print "total lk: " + str(len(lk_cells.keys()))
print "total cells: " + str(total_cells)
print "total sampled: " + str(total_sampled)
print "finished!"
