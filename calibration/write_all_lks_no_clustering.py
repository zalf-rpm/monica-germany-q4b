import os
import csv

script_path = os.path.dirname(os.path.abspath(__file__))

in_run_cells_file = "sample_setup_fullcover_SM.csv"#"sample_setup_fullcover_WW.csv"#
out_clustering_file = "no_clustering_lks_SM.csv"#"individual_lks_SM.csv"#"individual_lks_WW.csv"# "no_clustering_lks_WW.csv"#

one_lk_per_cluster = False

lk_set = set()

#take only the lk containing "use-it" cells
with open(script_path + "/" + in_run_cells_file) as _:
    reader = csv.reader(_)
    reader.next()
    for row in reader:
        if row[4].upper() == "Y":
            lk = int(row[1])
            lk_set.add(lk)

#all the lk refer to the same cluster (for testing and calibration of sowing rules)
with open(script_path + "/" + out_clustering_file, "wb") as _:
    writer = csv.writer(_)
    header = ["lk", "cluster"]
    writer.writerow(header)
    for lk in lk_set:
        row = [lk, 0]
        if one_lk_per_cluster:
            row = [lk, lk]
        writer.writerow(row)

print "finished"
            