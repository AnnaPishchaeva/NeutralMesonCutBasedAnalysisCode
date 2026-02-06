# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3
# This code was written by Alica Enderich (Febuary 2024)
# This code was modified and extended by Julia Schlägel (July 2024)
# This code was modified and extended by Anna Pishchaeva (October 2025)

import os, shutil
import math
import argparse
import copy
import numpy as np
import ctypes
import yaml
from copy import deepcopy
import ROOT
import datetime
ROOT.gROOT.SetBatch(True);

class Read_config:
    """_summary_: reads from the main config file (config_run)
    """
    def __init__(self, config_run_file):
        print("/////////////////////////////////////////////////")
        print("Read_config.py")
        if config_run_file:
            print("the main config file: ", config_run_file)
        with open(config_run_file, "r", encoding="utf-8") as config_run_yml:
            self.config_run = yaml.safe_load(config_run_yml)
        print("/////////////////////////////////////////////////")

    def get_my_dir_data(self):
        my_dir_data = self.config_run["common"]["my_dir_data"]
        return my_dir_data

    def get_my_dir(self):
        my_dir = self.config_run["common"]["my_dir"]
        return my_dir

    def get_input_file(self, typ):
        input_file_name = self.config_run["common"]["input_file_{0}".format(typ)]
        input_file = os.path.join(self.get_my_dir_data(), input_file_name)
        return input_file

    def get_filename_comp(self, obj):
        comp_name = self.config_run["compare"]["comparison_file_{0}".format(obj)]
        filename_comp = os.path.join(self.get_my_dir_data(), comp_name)
        return filename_comp

    def get_energy_of_other_file(self):
        energy_of_other_file = self.config_run["compare"]["energy_of_other_file"]
        return energy_of_other_file

    def get_system_of_other_file(self):
        system_of_other_file = self.config_run["compare"]["system_of_other_file"]
        return system_of_other_file

    def get_obj_comp(self):
        comp_name = self.config_run["compare"]["obj_name"]
        return comp_name

    def get_path_in_filename_comp(self, obj):
        path = self.config_run["compare"]["dir_to_histo_{0}".format(obj)]
        return path

    def get_histoname_in_filename_comp(self, obj):
        histoname_comp = self.config_run["compare"]["histo_{0}_name".format(obj)]
        return histoname_comp

    def get_comp_system(self):
        system_comp = self.config_run["compare"]["system"]
        return system_comp

    def get_icentrality_comp(self):
        icentrality = self.config_run["compare"]["icentrality"]
        return icentrality

    def get_config_file(self, meson, typ, decay):
        config_file = self.config_run["common"]["name_config_{0}_{1}_{2}".format(meson,typ,decay)]
        return config_file

    def get_output_folder(self):
        output_folder_name = self.config_run["common"]["output_folder"]
        energy = self.get_energy()
        if " " in energy:
            energy_list = energy.split(" ")
            energy = energy_list[0]+"_"+energy_list[1]
        output_folder = str(self.get_system()) + "_" +  str(energy)+ "_" + output_folder_name
        return output_folder

    def get_start_pt_array(self):
        start_pt_array = self.config_run["common"]["start_pt_array"]
        return start_pt_array

    def get_meson_array(self):
        meson_array =self.config_run["common"]["meson_array"]
        return meson_array

    def get_type_array(self):
        type_array = self.config_run["common"]["type_array"]
        return type_array

    def get_decay_array(self):
        decay_array = self.config_run["common"]["decay_array"]
        return decay_array

    def get_date(self):
        date = self.config_run["common"]["date"]
        return date

    def is_diff_cent(self):
        is_diff_cent = self.config_run["common"]["diff_cent"]
        return is_diff_cent

    def is_diff_RV0(self):
        is_diff_RV0 = self.config_run["common"]["diff_RV0"]
        return is_diff_RV0

    def is_compare_RV0(self):
        is_compare_RV0 = self.config_run["common"]["compare_RV0"]
        return is_compare_RV0

    def is_compare_occupancy(self):
        is_compare_occupancy = self.config_run["common"]["compare_occupancy"]
        return is_compare_occupancy

    def get_all_cent(self):
        all_cent = self.config_run["common"]["all_cent"]
        return all_cent

    def get_suffix(self):
        suffix = self.config_run["common"]["suffix"]
        return suffix

    def get_energy(self):
        energy = self.config_run["common"]["energy"]
        return energy

    def get_system(self):
        system = self.config_run["common"]["system"]
        return system

    def is_inv_mass_raw_yield(self):
        inv_mass_raw_yield = self.config_run["statements"]["inv_mass_raw_yield"]
        return inv_mass_raw_yield

    def is_eff_acc_corr_yield(self):
        eff_acc_corr_yield = self.config_run["statements"]["eff_acc_corr_yield"]
        return eff_acc_corr_yield

    def is_comp_dalitz_pcm(self):
        comp_dalitz_pcm = self.config_run["statements"]["comp_dalitz_pcm"]
        return comp_dalitz_pcm

    def is_comp_eta_pi0(self):
        comp_eta_pi0 = self.config_run["statements"]["comp_eta_pi0"]
        return comp_eta_pi0

    def is_comparison_to_other_data(self):
        comparison_to_other_data = self.config_run["compare"]["comparison_to_other_data"]
        return comparison_to_other_data

    def is_ratio_comp(self):
        is_ratio_comp = self.config_run["compare"]["is_ratio"]
        return is_ratio_comp






