import json
import sys
import monica_io
import zmq
import csv
import os
from datetime import date, timedelta
import collections
import threading
from threading import Thread
from collections import defaultdict
from copy import deepcopy
import numpy as np

script_path = os.path.dirname(os.path.abspath(__file__))

class monica_adapter(object):
    def __init__(self, sim_map, observations, calibration_target, dumped_envs_folder):

        self.obslist = [] #observations data structures for spotpy; keys: lk, year, var_name; vals: val(s)
        self.custom_events = defaultdict(list) # lk -> list of events needed
        self.sim_map = sim_map        
        self.lk_2_template_params = dict() #data structures for parameter templates 
        self.calibration_target = calibration_target
        self.dumped_envs_folder = dumped_envs_folder
        self.obslist_structure = [] #used later on to return the sim data needed in the same order as in obslist
        self.outputs = defaultdict(lambda: defaultdict())

        #count expected outs
        self.n_expected_out = 0
        for lk, sample in self.sim_map.iteritems():
            for cell in sample:
                self.n_expected_out += 1        
        
        #fill params template
        for lk, sample in self.sim_map.iteritems():
            #load one dumped env as template for params
            row_col = sample[0]
            fname = str(row_col[0]) + "_" + str(row_col[1]) + ".json"
            with open(dumped_envs_folder + "/" + fname) as _:
                template_env = json.load(_)
                #the following holds species, cv and sowing params (uniform for all the lk)
                self.lk_2_template_params[lk] = template_env["cropRotation"][0]["worksteps"][0]  

        #load custom events (for MONICA output)
        with open(script_path + "/template_events.json") as _:
            template_events = json.load(_)
        
        #add default custom events (sowing and yield)
        for lk in self.sim_map.keys():
            self.custom_events[lk].append(template_events["sowing_yield_template"][0])
            self.custom_events[lk].append(template_events["sowing_yield_template"][1])

        #build obslist and append pheno custom events:
        def append_pheno_spec(lk, year, doy):
            template = deepcopy(template_events["date_stage_template"])
            my_date = date(year-1, 12, 31) + timedelta(days=doy)
            template[0] = template[0].replace("yyyy-mm-dd", "{:04d}-{:02d}-{:02d}".format(my_date.year, my_date.month, my_date.day))
            self.custom_events[lk].append(template[0])
            self.custom_events[lk].append(template[1])
            return template[0]

        for lk in sorted(observations.keys()):
            #customize data structure for the output
            self.outputs[lk]["dates"] = set()
            self.outputs[lk]["vals"] = []
            for year in sorted(observations[lk].keys()):
                add_yr = False
                for var in sorted(observations[lk][year].keys()):
                    if calibration_target == "sowing" and var != 1:
                        continue
                    if calibration_target == "phenology":
                        #recorded phenological phase will be compared to the 
                        #phase simulated by monica on the same day
                        self.obslist.append(var)
                        record_date = append_pheno_spec(lk, year, observations[lk][year][var])
                        self.obslist_structure.append((lk, record_date))
                    else:
                        #sowing doy or yield will be compared to the simulated ones
                        self.obslist.append(observations[lk][year][var])
                        self.obslist_structure.append((lk, year))
        
        self.context = zmq.Context()
        self.socket_producer = self.context.socket(zmq.PUSH)
        self.socket_producer.connect("tcp://cluster1:6666")
        #self.socket_producer.connect("tcp://localhost:6666")

    def run(self,args):
        return self._run(*args)

    def _run(self, vector, user_params):

        def seek_set_param(par, p_value, model_params):
            p_name = par["name"]
            array = par["array"]
            add_index = False
            if isinstance(model_params[p_name], int) or isinstance(model_params[p_name], float):
                add_index = False
            elif len(model_params[p_name]) > 1 and isinstance(model_params[p_name][1], basestring):
                add_index = True #the param contains text (e.g., units)
            if array.upper() == "FALSE":
                if add_index:
                    model_params[p_name][0] = p_value
                else:
                    model_params[p_name] = p_value
            else: #param is in an array (possibly nested)
                array = array.split("_") #nested array
                if add_index:
                    array = [0] + array
                if len(array) == 1:
                    model_params[p_name][int(array[0])] = p_value
                elif len(array) == 2:
                    model_params[p_name][int(array[0])][int(array[1])] = p_value
                elif len(array) == 3:
                    model_params[p_name][int(array[0])][int(array[1])][int(array[2])] = p_value
                else:
                    print "param array too nested, contact developers"
            
        for lk, template_params in self.lk_2_template_params.iteritems():
            species_params = template_params["crop"]["cropParams"]["species"]
            cultivar_params = template_params["crop"]["cropParams"]["cultivar"]
            #set params according to spotpy sampling
            for i in range(len(user_params)):
                if user_params[i]["name"] in template_params:
                    #sowing params
                    seek_set_param(user_params[i],
                    user_params[i]["derive_function"](vector, template_params) if "derive_function" in user_params[i] else vector[i],
                    template_params)
                elif user_params[i]["name"] in species_params:
                    seek_set_param(user_params[i],
                    user_params[i]["derive_function"](vector, species_params) if "derive_function" in user_params[i] else vector[i],
                    species_params)
                elif user_params[i]["name"] in cultivar_params:
                    seek_set_param(user_params[i],
                    user_params[i]["derive_function"](vector, cultivar_params) if "derive_function" in user_params[i] else vector[i],
                    cultivar_params)
                else:
                    print(str(user_params[i]["name"]) + " not found, please revise and restart")
                    exit()

        #launch parallel thread for the collector
        collector = Thread(target=self.collect_results)
        collector.daemon = True
        collector.start()

        #launch sims from dumped envs
        for lk, sample in self.sim_map.iteritems():
            for cell in sample:
                fname = str(cell[0]) + "_" + str(cell[1]) + ".json"
                with open(self.dumped_envs_folder + "/" + fname) as _:
                    env = json.load(_)
                #update params
                env["cropRotation"][0]["worksteps"][0] = self.lk_2_template_params[lk]
                #update events (output)
                env["events"] = self.custom_events[lk]

                self.socket_producer.send_json(env)
                #print("sent env for lk " + str(lk) + " sampled cell " + str(cell))

        #wait until the collector finishes
        collector.join()
        
        #calculate averages (from sampled cells) for lk
        lk_dates2index = defaultdict(dict)
        for lk in self.outputs.keys():
            out_arr = np.array(self.outputs[lk]["vals"])
            self.outputs[lk]["avg_vals"] = np.average(out_arr, axis=0)
            #add convenience dict to easily get the index for any given lk and date
            dates_lk = list(self.outputs[lk]["dates"])[0]
            for i in range(len(dates_lk)):
                lk_dates2index[lk][dates_lk[i]] = i

        #build the evaluation list for spotpy
        evallist = []
        for lk, record_date in self.obslist_structure:
            date_index = lk_dates2index[lk][record_date]
            target_out = self.outputs[lk]["avg_vals"][date_index]
            evallist.append(target_out)

        return evallist

        
    def collect_results(self):
        socket_collector = self.context.socket(zmq.PULL)
        socket_collector.connect("tcp://cluster1:7777")
        #socket_collector.connect("tcp://localhost:7777")
        received_results = 0
        leave = False
        while not leave:
            try:
                rec_msg = socket_collector.recv_json()
            except:
                continue            
            
            dates = []
            vals = []
            lk = rec_msg["customId"]["lk_id"]                
            for res in rec_msg["data"]:
                if "crop" in res["origSpec"]:
                    #sowing and yield data
                    if self.calibration_target == "sowing":
                        dates = res["results"][0]
                        vals = res["results"][1]
                    elif self.calibration_target == "yield":
                        dates = res["results"][2]
                        vals = res["results"][3]
                else:
                    #phenology (one date and val added at a time)
                    dates.append(res["origSpec"].replace('"',''))
                    vals.append(res["results"][0][0])
            
            self.outputs[lk]["dates"].add(tuple(dates))
            self.outputs[lk]["vals"].append(vals)

            if len(self.outputs[lk]["dates"]) > 1:
                print ("mismatch in dates after adding results from " + str((rec_msg["customId"]["srow"], rec_msg["customId"]["scol"])))
                exit()            
            
            #print ("received " + str(rec_msg["customId"]))
            received_results += 1
            #print("total received: " + str(received_results))

            if received_results == self.n_expected_out:
                leave = True
