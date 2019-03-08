#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

import time
import os
import math
import json
import csv
import copy
from StringIO import StringIO
from datetime import date, timedelta
from collections import defaultdict
import sys
import zmq

import sqlite3
import numpy as np
from pyproj import Proj, transform

import monica_io
import soil_io
import monica_run_lib as Mrunlib
import numpy as np


# chose setup for running as container (in docker) or local run 
# for local run you need a local monica server running e.g "monica-zmq-server -bi -i tcp://*:6666 -bo -o tcp://*:7777"
# if you use a local monica-zmq-server, set "archive-path-to-climate-dir" to the folder where your climate data is found

# or a local docker image:
# docker run -p 6666:6666 -p 7777:7777 --env monica_instances=3 --rm --name monica_test -v climate-data:/monica_data/climate-data zalfrpm/monica-cluster:latest
# if you use docker, set "archive-path-to-climate-dir" = "/monica_data/climate-data/" 
# and create a volume for the climate data, e.g for a network drive
# docker volume create --driver local \
#     --opt type=cifs \
#     --opt device='//network_drive_ip/archiv-daten/md/data/climate/' \
#     --opt o='username=your_username,password=your_password' \
# climate-data


USER_MODE = "localProducer-localMonica"
#USER_MODE = "localProducer-remoteMonica"
#USER_MODE = "container"

PATHS = {
    # adjust the local path to your environment
    "localProducer-localMonica": {
        "producer-path-to-climate-dir": "Z:/data/climate/dwd/csvs/", # mounted path to archive or hard drive with climate data; used only for latlon file
        "producer-path-to-data-dir": "C:/Users/stella/Documents/GitHub/monica-germany-q4b/monica-data/data/", # mounted path to archive or hard drive with data 
        "producer-path-to-projects-dir": "C:/Users/stella/Documents/GitHub/monica-germany-q4b/monica-data/projects/", # mounted path to archive or hard drive with project data 
        "monica-parameters-path": "C:/Users/stella/Documents/GitHub/monica-parameters/", # path to monica-parameters
        "monica-path-to-climate-dir": "Z:/data/climate/dwd/csvs/germany/", # mounted path to archive accessable by monica executable
    },
    "localProducer-remoteMonica": {
        "producer-path-to-climate-dir": "Z:/data/climate/dwd/csvs/", # mounted path to archive or hard drive with climate data; used only for latlon file
        "producer-path-to-data-dir": "C:/Users/stella/Documents/GitHub/monica-germany-q4b/monica-data/data/", # mounted path to archive or hard drive with data 
        "producer-path-to-projects-dir": "C:/Users/stella/Documents/GitHub/monica-germany-q4b/monica-data/projects/", # mounted path to archive or hard drive with project data 
        "monica-parameters-path": "C:/Users/stella/Documents/GitHub/monica-parameters/", # path to monica-parameters
        "monica-path-to-climate-dir": "/monica_data/climate-data/dwd/csvs/germany/", # mounted path to archive accessable by monica executable
    },
    "container": {
        "producer-path-to-climate-dir": "/monica_data/climate-data/", # needs to be mounted there
        "producer-path-to-data-dir": "/monica_data/data/", # needs to be mounted there
        "producer-path-to-projects-dir": "/monica_data/project/", # needs to be mounted there
        "monica-parameters-path": "/home/monica-parameters/", # monica parameter location in docker image
        "monica-path-to-climate-dir": "/monica_data/climate-data/",  # mounted path to archive on cluster docker image 
    }
}

CONFIGURATION = {
        "mode": USER_MODE,
        "port": "6666",
        "server": "localhost",
        "start-row": "0", 
        "end-row": "-1",
        "sim.json": "sim.json",
        "crop.json": "crop.json",
        "site.json": "site.json",
        "setups-file": "sim_setups.csv",
        "run-setups": "[2]",
        "shared_id": None,
    }

PROJECT_FOLDER = "monica-germany/"
DATA_SOIL_DB = "germany/buek1000.sqlite"
DATA_GRID_HEIGHT = "germany/dem_1000_gk5.asc" 
DATA_GRID_SLOPE = "germany/slope_1000_gk5.asc"
DATA_GRID_LAND_USE = "germany/corine2006_1000_gk5.asc"
DATA_GRID_LANDKREIS = "germany/landkreise_1000_gk3.asc"
DATA_GRID_SOIL = "germany/buek1000_1000_gk5.asc"
TEMPLATE_PATH_LATLON = "{path_to_climate_dir}/latlon_to_rowcol.json" #"{path_to_climate_dir}{climate_data}/csvs/latlon-to-rowcol.json"
TEMPLATE_PATH_HARVEST = "{path_to_projects_dir}{project_folder}ILR_SEED_HARVEST_doys_{crop_id}.csv"
TEMPLATE_PATH_CLIMATE_CSV = "row-{crow}/col-{ccol}.csv"
GEO_TARGET_GRID="epsg:31469" #proj4 -> 3-degree gauss-kruger zone 5 (=Germany) https://epsg.io/31469

DEBUG_DONOT_SEND = False
DEBUG_WRITE = False
DEBUG_ROWS = 10
DEBUG_WRITE_FOLDER = "./debug_out"
DEBUG_WRITE_CLIMATE = False

DUMP_ENVS_CALIB = True
if DUMP_ENVS_CALIB:
    DEBUG_DONOT_SEND = True #no need to run monica
    script_path = os.path.dirname(os.path.abspath(__file__))

# commandline parameters e.g "server=localhost port=6666 shared_id=2"
def run_producer(config):
    "main"

    print "config:", config

    context = zmq.Context()
    socket = context.socket(zmq.PUSH)
    #config_and_no_data_socket = context.socket(zmq.PUSH)

    # select paths 
    paths = PATHS[config["mode"]]
    # open soil db connection
    soil_db_con = sqlite3.connect(paths["producer-path-to-data-dir"] + DATA_SOIL_DB)
    # connect to monica proxy (if local, it will try to connect to a locally started monica)
    socket.connect("tcp://" + config["server"] + ":" + str(config["port"]))

    # read setup from csv file
    setups = Mrunlib.read_sim_setups(paths["producer-path-to-projects-dir"] + PROJECT_FOLDER + config["setups-file"])
    run_setups = json.loads(config["run-setups"])
    print "read sim setups: ", paths["producer-path-to-projects-dir"] + PROJECT_FOLDER + config["setups-file"]

    #transforms geospatial coordinates from one coordinate reference system to another
    # transform wgs84 into gk5
    wgs84 = Proj(init="epsg:4326") #proj4 -> (World Geodetic System 1984 https://epsg.io/4326)
    gk5 = Proj(init=GEO_TARGET_GRID) 
        

    # Load grids
    ## note numpy is able to load from a compressed file, ending with .gz or .bz2
    
    # height data for germany
    path_to_dem_grid = paths["producer-path-to-data-dir"] + DATA_GRID_HEIGHT 
    dem_metadata, _ = Mrunlib.read_header(path_to_dem_grid)
    dem_grid = np.loadtxt(path_to_dem_grid, dtype=int, skiprows=6)
    dem_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(dem_grid, dem_metadata)
    print "read: ", path_to_dem_grid
    
    # slope data
    path_to_slope_grid = paths["producer-path-to-data-dir"] + DATA_GRID_SLOPE
    slope_metadata, _ = Mrunlib.read_header(path_to_slope_grid)
    slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
    slope_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(slope_grid, slope_metadata)
    print "read: ", path_to_slope_grid

    # land use data
    path_to_corine_grid = paths["producer-path-to-data-dir"] + DATA_GRID_LAND_USE
    corine_meta, _ = Mrunlib.read_header(path_to_corine_grid)
    corine_grid = np.loadtxt(path_to_corine_grid, dtype=int, skiprows=6)
    corine_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(corine_grid, corine_meta)
    print "read: ", path_to_corine_grid

    # landkreis data
    gk3 = Proj(init="epsg:3396")
    path_to_landkreis_grid = paths["producer-path-to-data-dir"] + DATA_GRID_LANDKREIS
    landkreis_meta, _ = Mrunlib.read_header(path_to_landkreis_grid)
    landkreis_grid = np.loadtxt(path_to_landkreis_grid, dtype=int, skiprows=6)
    landkreis_gk3_interpolate = Mrunlib.create_ascii_grid_interpolator(landkreis_grid, landkreis_meta)
    print "read: ", path_to_landkreis_grid

    # soil data
    path_to_soil_grid = paths["producer-path-to-data-dir"] + DATA_GRID_SOIL
    soil_metadata, _ = Mrunlib.read_header(path_to_soil_grid)
    soil_grid = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
    print "read: ", path_to_soil_grid

    cdict = {}
    path = TEMPLATE_PATH_LATLON.format(path_to_climate_dir=paths["producer-path-to-climate-dir"])
    climate_data_to_gk5_interpolator = Mrunlib.create_climate_geoGrid_interpolator_from_json_file(path, wgs84, gk5, cdict)

    sent_env_count = 1
    start_time = time.clock()

    listOfClimateFiles = set()
    # run calculations for each setup
    for _, setup_id in enumerate(run_setups):

        if setup_id not in setups:
            continue
        start_setup_time = time.clock()

        setup = setups[setup_id]
        climate_data = setup["climate_data"]
        climate_model = setup["climate_model"]
        climate_scenario = setup["climate_scenario"]
        climate_region = setup["climate_region"]
        crop_id = setup["crop-id"]

        if DUMP_ENVS_CALIB:
            #identify sample cells            
            sample_cells = set()
            sample_cells_file = "sample_setup_partialcover.csv" if setup["landcover"] else "sample_setup_fullcover.csv"
            with open(script_path + "/calibration/" + sample_cells_file) as _:
                reader = csv.reader(_)
                header = reader.next()
                data_fields = {}
                for i in range(len(header)):
                    data_fields[header[i]] = i
                for row in reader:
                    row_col = (int(row[data_fields["row"]]), int(row[data_fields["col"]]))
                    sample_cells.add(row_col)

        #read whitelisted landkreise (i.e., those having at least some information for both phenology and yield)
        whitelist_lk = set()
        whitelist_file = ""
        if crop_id == "SM":
            whitelist_file = "whitelisted_lks_silage_maize.csv"
        elif crop_id == "WW":
            whitelist_file = "whitelisted_lks_winter_wheat.csv"
        else:
            print "no whitelist_file for crop id: " + crop_id
        with open(paths["producer-path-to-projects-dir"] + PROJECT_FOLDER + whitelist_file) as _:
            print "reading " + whitelist_file
            reader = csv.reader(_)
            header = reader.next()
            for row in reader:
                whitelist_lk.add(int(row[1]))


        #read sowing and harvest data
        #keys: lk, sow/harv, vals: doy
        lk_mgt_data = defaultdict(lambda: defaultdict(list))
        pheno_file = ""
        if crop_id == "SM":
            pheno_file = "lk_cleaned_pheno_silage_maize_1999_2017.csv"
        elif crop_id == "WW":
            pheno_file = "lk_cleaned_pheno_winter_wheat_1999_2017.csv"
        else:
            print "no pheno data for crop id: " + crop_id
        with open(paths["producer-path-to-projects-dir"] + PROJECT_FOLDER + pheno_file) as _:
            print "reading " + pheno_file
            data_fields = {}
            reader = csv.reader(_)
            header = reader.next()
            for i in range(len(header)):
                data_fields[header[i]] = i
            for row in reader:
                lk = int(float(row[data_fields["lk"]]))
                if float(row[data_fields["phase"]]) == 10:
                    lk_mgt_data[lk]["sowing"].append(float(row[data_fields["avg_doy"]]))
                if float(row[data_fields["phase"]]) == 24:
                    lk_mgt_data[lk]["harvest"].append(float(row[data_fields["avg_doy"]]))

        # read template sim.json 
        with open(setup.get("sim.json", config["sim.json"])) as _:
            sim_json = json.load(_)
        # change start and end date acording to setup
        if setup["start_year"]:
            sim_json["climate.csv-options"]["start-date"] = str(setup["start_year"]) + "-01-01"
        if setup["end_year"]:
            sim_json["climate.csv-options"]["end-date"] = str(setup["end_year"]) + "-12-31" 
        sim_json["include-file-base-path"] = paths["monica-parameters-path"]

        # read template site.json 
        with open(setup.get("site.json", config["site.json"])) as _:
            site_json = json.load(_)
        # read template crop.json
        with open(setup.get("crop.json", config["crop.json"])) as _:
            crop_json = json.load(_)

        # create environment template from json templates
        env_template = monica_io.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""
        })
        # set shared id in template
        if config["shared_id"]:
            env_template["sharedId"] = config["shared_id"]       

        # create crop rotation according to setup
        # clear crop rotation and get its template
        crop_rotation_templates = env_template.pop("cropRotation")
        env_template["cropRotation"] = []
        # get correct template
        env_template["cropRotation"] = crop_rotation_templates[crop_id]

        # we just got one cultivation method in our rotation
        worksteps_templates_dict = env_template["cropRotation"][0].pop("worksteps")

        # clear the worksteps array and rebuild it out of the setup      
        worksteps = env_template["cropRotation"][0]["worksteps"] = []
        worksteps.append(worksteps_templates_dict["sowing"][setup["sowing-date"]])
        worksteps.append(worksteps_templates_dict["harvest"][setup["harvest-date"]])

        scols = int(soil_metadata["ncols"])
        srows = int(soil_metadata["nrows"])
        scellsize = int(soil_metadata["cellsize"])
        xllcorner = int(soil_metadata["xllcorner"])
        yllcorner = int(soil_metadata["yllcorner"])

        ###############################test
        #lk_grid = np.empty([srows, scols], dtype=int)
        ##############################

        print "All Rows x Cols: " + str(srows) + "x" + str(scols) 
        for srow in xrange(0, srows):
           
            try:
                print srow,
            except Exception:
                # no out
                print srow,

            if srow < int(config["start-row"]):
                continue
            elif int(config["end-row"]) > 0 and srow > int(config["end-row"]):
                break

            for scol in xrange(0, scols):

                ##################test
                #soil_id = soil_grid[srow, scol]
                #if soil_id == -9999:
                #    lk_grid[srow, scol] = -9999
                #else:
                #    sh_gk5 = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
                #    sr_gk5 = xllcorner + (scellsize / 2) + scol * scellsize
                #    #inter = crow/ccol encoded into integer
                #    crow, ccol = climate_data_to_gk5_interpolator(sr_gk5, sh_gk5)
                #    
                #    sr_gk3, sh_gk3 = transform(gk5, gk3, sr_gk5, sh_gk5)
                #    lk_grid[srow, scol] = landkreis_gk3_interpolate(sr_gk3, sh_gk3)
                #continue
                ##############################
                

                soil_id = soil_grid[srow, scol]
                if soil_id == -9999:
                    continue
                if soil_id < 1 or soil_id > 71:
                    #print "row/col:", srow, "/", scol, "has unknown soil_id:", soil_id
                    #unknown_soil_ids.add(soil_id)
                    continue
                
                if DUMP_ENVS_CALIB:
                    if (srow, scol) not in sample_cells:
                        continue
                
                #get coordinate of clostest climate element of real soil-cell
                sh_gk5 = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
                sr_gk5 = xllcorner + (scellsize / 2) + scol * scellsize
                #inter = crow/ccol encoded into integer
                crow, ccol = climate_data_to_gk5_interpolator(sr_gk5, sh_gk5)

                sr_gk3, sh_gk3 = transform(gk5, gk3, sr_gk5, sh_gk5)

                # check if current grid cell is used for agriculture                
                if setup["landcover"]:
                    corine_id = corine_gk5_interpolate(sr_gk5, sh_gk5)
                    if corine_id not in [200, 210, 211, 212, 240, 241, 242, 243, 244]:
                        continue

                landkreis_id = landkreis_gk3_interpolate(sr_gk3, sh_gk3)
                # check if landkreis is whitelisted for the study
                if landkreis_id not in whitelist_lk:
                    print "lk " + str(landkreis_id) + " not whitelisted, skipping it"
                    continue
                #print "row/col:", srow, "/", scol, "->", landkreis_id
                
                height_nn = dem_gk5_interpolate(sr_gk5, sh_gk5)
                slope = slope_gk5_interpolate(sr_gk5, sh_gk5)
                                
                #print "scol:", scol, "crow/col:", (crow, ccol), "soil_id:", soil_id, "height_nn:", height_nn, "slope:", slope, "seed_harvest_cs:", seed_harvest_cs

                clat, _ = cdict[(crow, ccol)]
                
                # set external sowing/harvest dates               
                avg_sowing_doy = int(np.average(lk_mgt_data[landkreis_id]["sowing"]))
                earliest_sowing_doy = int(np.min(lk_mgt_data[landkreis_id]["sowing"]))
                latest_sowing_doy = int(np.max(lk_mgt_data[landkreis_id]["sowing"]))

                avg_harvest_doy = int(np.average(lk_mgt_data[landkreis_id]["harvest"]))
                latest_harvest_doy = int(np.max(lk_mgt_data[landkreis_id]["harvest"]))

                #winter crops: avoid overlapping between latest harvest and earliest sowing
                if worksteps[0]["crop"]["is-winter-crop"]:
                    if latest_harvest_doy > earliest_sowing_doy:
                        print "adjusting latest harvest and earliest sowing to avoid overlapping in lk " + str(landkreis_id)
                        correction = int((latest_harvest_doy - earliest_sowing_doy)/2 +1)
                        latest_harvest_doy -= correction
                        earliest_sowing_doy += correction
                
                #set sowing
                if setup["sowing-date"] == "fixed":
                    relative_year = worksteps[0]["date"].split("-")[0]
                    avg_date = date(2018, 12, 31) + timedelta(days=avg_sowing_doy) 
                    worksteps[0]["date"] = relative_year + "-{:02d}-{:02d}".format(avg_date.month, avg_date.day)                         
                elif setup["sowing-date"] == "auto":
                    relative_year = worksteps[0]["latest-date"].split("-")[0]
                    earliest_date = date(2018, 12, 31) + timedelta(days=earliest_sowing_doy)
                    latest_date = date(2018, 12, 31) + timedelta(days=latest_sowing_doy)
                    worksteps[0]["earliest-date"] = relative_year + "-{:02d}-{:02d}".format(earliest_date.month, earliest_date.day)
                    worksteps[0]["latest-date"] = relative_year + "-{:02d}-{:02d}".format(latest_date.month, latest_date.day)
                else:
                    print "unable to handle sowing option"
                
                #set harvest
                if setup["harvest-date"] == "fixed":
                    relative_year = worksteps[1]["date"].split("-")[0]
                    avg_date = date(2018, 12, 31) + timedelta(days=avg_harvest_doy) 
                    worksteps[1]["date"] = relative_year + "-{:02d}-{:02d}".format(avg_date.month, avg_date.day)                         
                elif setup["harvest-date"] == "auto":
                    relative_year = worksteps[1]["latest-date"].split("-")[0]
                    latest_date = date(2018, 12, 31) + timedelta(days=latest_harvest_doy)
                    worksteps[1]["latest-date"] = relative_year + "-{:02d}-{:02d}".format(latest_date.month, latest_date.day)
                else:
                    print "unable to handle harvest option"
                
                #with open("dump-" + str(c) + ".json", "w") as jdf:
                #    json.dump({"id": (str(resolution) \
                #        + "|" + str(vrow) + "|" + str(vcol) \
                #        + "|" + str(crow) + "|" + str(ccol) \
                #        + "|" + str(soil_id) \
                #        + "|" + crop_id \
                #        + "|" + str(uj_id)), "sowing": worksteps[0], "harvest": worksteps[1]}, jdf, indent=2)
                #    c += 1

                env_template["params"]["userCropParameters"]["__enable_T_response_leaf_expansion__"] = setup["LeafExtensionModifier"]

                # set soil-profile
                sp_json = soil_io.soil_parameters(soil_db_con, soil_id)
                soil_profile = monica_io.find_and_replace_references(sp_json, sp_json)["result"]
                    
                #print "soil:", soil_profile

                env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

                # setting groundwater level
                if setup["groundwater-level"]:
                    groundwaterlevel = 20
                    layer_depth = 0
                    for layer in soil_profile:
                        if layer.get("is_in_groundwater", False):
                            groundwaterlevel = layer_depth
                            #print "setting groundwaterlevel of soil_id:", str(soil_id), "to", groundwaterlevel, "m"
                            break
                        layer_depth += Mrunlib.get_value(layer["Thickness"])
                    env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                    env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [max(0, groundwaterlevel - 0.2) , "m"]
                    env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [groundwaterlevel + 0.2, "m"]
                    
                # setting impenetrable layer
                if setup["impenetrable-layer"]:
                    impenetrable_layer_depth = Mrunlib.get_value(env_template["params"]["userEnvironmentParameters"]["LeachingDepth"])
                    layer_depth = 0
                    for layer in soil_profile:
                        if layer.get("is_impenetrable", False):
                            impenetrable_layer_depth = layer_depth
                            #print "setting leaching depth of soil_id:", str(soil_id), "to", impenetrable_layer_depth, "m"
                            break
                        layer_depth += Mrunlib.get_value(layer["Thickness"])
                    env_template["params"]["userEnvironmentParameters"]["LeachingDepth"] = [impenetrable_layer_depth, "m"]
                    env_template["params"]["siteParameters"]["ImpenetrableLayerDepth"] = [impenetrable_layer_depth, "m"]

                if setup["elevation"]:
                    env_template["params"]["siteParameters"]["heightNN"] = height_nn

                if setup["slope"]:
                    env_template["params"]["siteParameters"]["slope"] = slope / 100.0

                if setup["latitude"]:
                    env_template["params"]["siteParameters"]["Latitude"] = clat

                if setup["CO2"]:
                    env_template["params"]["userEnvironmentParameters"]["AtmosphericCO2"] = float(setup["CO2"])

                if setup["O3"]:
                    env_template["params"]["userEnvironmentParameters"]["AtmosphericO3"] = float(setup["O3"])

                env_template["params"]["simulationParameters"]["UseNMinMineralFertilisingMethod"] = setup["fertilization"]
                env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = setup["irrigation"]

                env_template["params"]["simulationParameters"]["NitrogenResponseOn"] = setup["NitrogenResponseOn"]
                env_template["params"]["simulationParameters"]["WaterDeficitResponseOn"] = setup["WaterDeficitResponseOn"]
                env_template["params"]["simulationParameters"]["EmergenceMoistureControlOn"] = setup["EmergenceMoistureControlOn"]
                env_template["params"]["simulationParameters"]["EmergenceFloodingControlOn"] = setup["EmergenceFloodingControlOn"]

                env_template["csvViaHeaderOptions"] = sim_json["climate.csv-options"]
                
                
                subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(crow=str(crow),
                                                                ccol=str(ccol))

                # subpath_to_csv = climate_data + "/csvs/" \
                # + (climate_model + "/" if climate_model else "") \
                # + (climate_scenario + "/" if climate_scenario else "") \
                # + climate_region + "/row-" + str(crow) + "/col-" + str(ccol) + ".csv"
                env_template["pathToClimateCSV"] = paths["monica-path-to-climate-dir"] + subpath_to_csv
                #print env_template["pathToClimateCSV"]
                if DEBUG_WRITE_CLIMATE :
                    listOfClimateFiles.add(subpath_to_csv)

                env_template["customId"] = {
                    "setup_id": setup_id,
                    "srow": srow, "scol": scol,
                    "crow": crow, "ccol": ccol,
                    "soil_id": soil_id,
                    "crop_id": crop_id,
                    "lk_id": landkreis_id
                }

                if DUMP_ENVS_CALIB:
                    #folder structure:
                    #/calibration/dumped_envs
                    #               |--mode
                    #                   |--setup_ID
                    path_to_dump_dir = script_path + "/calibration/dumped_envs/" + str(config["mode"]) + "/" + str(setup_id)
                    if not os.path.exists(path_to_dump_dir):
                        os.makedirs(path_to_dump_dir)
                    dump_file = str(srow) + "_" + str(scol) + ".json"
                    with open(path_to_dump_dir + "/" + dump_file, "w") as _:
                        _.write(json.dumps(env_template, indent=4))

                if not DEBUG_DONOT_SEND :
                    socket.send_json(env_template)
                    print "sent env ", sent_env_count, " customId: ", env_template["customId"]

                sent_env_count += 1

                # write debug output, as json file
                if DEBUG_WRITE:
                    if not os.path.exists(DEBUG_WRITE_FOLDER):
                        os.makedirs(DEBUG_WRITE_FOLDER)
                    if sent_env_count < DEBUG_ROWS  :

                        path_to_debug_file = DEBUG_WRITE_FOLDER + "/row_" + str(sent_env_count-1) + "_" + str(setup_id) + ".json" 

                        if not os.path.isfile(path_to_debug_file):
                            with open(path_to_debug_file, "w") as _ :
                                _.write(json.dumps(env_template))
                        else:
                            print "WARNING: Row ", (sent_env_count-1), " already exists"
            #print "unknown_soil_ids:", unknown_soil_ids

            #print "crows/cols:", crows_cols
        stop_setup_time = time.clock()
        print "Setup ", (sent_env_count-1), " envs took ", (stop_setup_time - start_setup_time), " seconds"

        ################################test
        #print "saving lk grid..."
        #np.savetxt('lk_grid.asc', lk_grid, delimiter='\t', fmt="%.0f")
        ###################################

    stop_time = time.clock()

    # write summary of used json files
    if DEBUG_WRITE_CLIMATE:
        if not os.path.exists(DEBUG_WRITE_FOLDER):
            os.makedirs(DEBUG_WRITE_FOLDER)

        path_to_climate_summary = DEBUG_WRITE_FOLDER + "/climate_file_list" + ".csv"
        with open(path_to_climate_summary, "w") as _:
            _.write('\n'.join(listOfClimateFiles))

    try:
        print "sending ", (sent_env_count-1), " envs took ", (stop_time - start_time), " seconds"
        #print "ran from ", start, "/", row_cols[start], " to ", end, "/", row_cols[end]
        print "exiting run_producer()"
    except Exception:
        raise

if __name__ == "__main__":

    config = copy.deepcopy(CONFIGURATION)

    # read commandline args only if script is invoked directly from commandline
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k in config:
                config[k] = v

    run_producer(config)