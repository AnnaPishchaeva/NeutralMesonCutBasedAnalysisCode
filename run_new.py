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
from Analysis import Analysis
from Read_config import Read_config
#_________________________________________________________________________________________
def run_new():
    """_summary_: reads from the main config file (config_run) using Read_config and starts the analysis using Analysis class
    """
    print("run_new() has started")
    config_run_file = "configs/config_run.yml"
    r_config = Read_config(config_run_file)
    if r_config.is_inv_mass_raw_yield():
        for is_start_pt in r_config.get_start_pt_array():
            for typ in r_config.get_type_array():
                input_file = r_config.get_input_file(typ)
                for meson in r_config.get_meson_array():
                    folder = os.path.join(r_config.get_my_dir_data(), r_config.get_output_folder(), "{0}".format(meson))
                    for decay in r_config.get_decay_array():
                        config_file = r_config.get_config_file(meson, typ, decay)
                        analysis = Analysis(r_config.is_diff_cent(), r_config.is_compare_RV0(), r_config.is_compare_occupancy(), r_config.is_diff_RV0(), r_config.get_all_cent(), is_start_pt, r_config.get_date(), folder, r_config.get_suffix(), r_config.get_energy(), r_config.get_system())
                        analysis.inv_mass_raw_yield(input_file, config_file, typ, decay, meson)

    if r_config.is_eff_acc_corr_yield():
        for is_start_pt in r_config.get_start_pt_array():
            input_file_data =  r_config.get_input_file("data")
            input_file_mc = r_config.get_input_file("mc")
            for meson in r_config.get_meson_array():
                folder = os.path.join(r_config.get_my_dir_data(), r_config.get_output_folder(), "{0}".format(meson))
                filename_comp = r_config.get_filename_comp(meson)
                for decay in r_config.get_decay_array():
                    config_file_data = r_config.get_config_file(meson, "data", decay)
                    config_file_mc = r_config.get_config_file(meson, "mc", decay)
                    analysis = Analysis(r_config.is_diff_cent(), r_config.is_compare_RV0(), r_config.is_compare_occupancy(), r_config.is_diff_RV0(), r_config.get_all_cent(), is_start_pt, r_config.get_date(), folder, r_config.get_suffix(), r_config.get_energy(), r_config.get_system())
                    analysis.eff_acc_corr_yield(input_file_mc, input_file_data, config_file_mc, config_file_data, decay, meson)

    if r_config.is_comp_dalitz_pcm():
        for is_start_pt in r_config.get_start_pt_array():
            for meson in r_config.get_meson_array():
                folder = os.path.join(r_config.get_my_dir_data(), r_config.get_output_folder(), "{0}".format(meson))
                config_file_pcm = r_config.get_config_file(meson, "data", "PCM")
                config_file_dalitz = r_config.get_config_file(meson, "data", "Dalitz")
                analysis = Analysis(r_config.is_diff_cent(), r_config.is_compare_RV0(), r_config.is_compare_occupancy(), r_config.is_diff_RV0(), r_config.get_all_cent(), is_start_pt, r_config.get_date(), folder, r_config.get_suffix(), r_config.get_energy(), r_config.get_system())
                analysis.comp_dalitz_pcm(config_file_pcm, config_file_dalitz, meson)

    if r_config.is_comp_eta_pi0():
        for is_start_pt in r_config.get_start_pt_array():
            for decay in r_config.get_decay_array():
                folder = os.path.join(r_config.get_my_dir_data(), r_config.get_output_folder(), "pi0_eta_comparison")
                config_file_pi0 = r_config.get_config_file("pi0", "data", decay)
                config_file_eta = r_config.get_config_file("eta", "data", decay)
                analysis = Analysis(r_config.is_diff_cent(), r_config.is_compare_RV0(), r_config.is_compare_occupancy(), r_config.is_diff_RV0(), r_config.get_all_cent(), is_start_pt, r_config.get_date(), folder, r_config.get_suffix(), r_config.get_energy(), r_config.get_system())
                analysis.comp_eta_pi0(config_file_pi0, config_file_eta, decay)

    if r_config.is_comparison_to_other_data():
        energy_of_other_file = r_config.get_energy_of_other_file()
        system_of_other_file = r_config.get_system_of_other_file()
        for is_start_pt in r_config.get_start_pt_array():
            obj_comp = r_config.get_obj_comp()
            if r_config.get_comp_system() == "eta_pi0":
                for decay in r_config.get_decay_array():
                    folder = os.path.join(r_config.get_my_dir_data(), r_config.get_output_folder(), "{0}".format("pi0"))
                    filename_comp = r_config.get_filename_comp(decay)
                    path_comp = r_config.get_path_in_filename_comp(decay)
                    histname_comp = r_config.get_histoname_in_filename_comp(decay)
                    config = r_config.get_config_file("pi0", "data", decay)
                    analysis = Analysis(r_config.is_diff_cent(), r_config.is_compare_RV0(), r_config.is_compare_occupancy(), r_config.is_diff_RV0(), r_config.get_all_cent(), is_start_pt, r_config.get_date(), folder, r_config.get_suffix(), r_config.get_energy(), r_config.get_system())
                    analysis.comp_to_other_file(energy_of_other_file, system_of_other_file, decay, "pi0", filename_comp, obj_comp, path_comp, histname_comp, r_config.is_ratio_comp(), config)
            if r_config.get_comp_system() == "Dalitz_PCM":
                for meson in r_config.get_meson_array():
                    folder = os.path.join(r_config.get_my_dir_data(), r_config.get_output_folder(), "{0}".format(meson))
                    filename_comp = r_config.get_filename_comp(meson)
                    path_comp = r_config.get_path_in_filename_comp(meson)
                    histname_comp = r_config.get_histoname_in_filename_comp(meson)
                    config = r_config.get_config_file(meson, "data", "PCM")
                    analysis = Analysis(r_config.is_diff_cent(), r_config.is_compare_RV0(), r_config.is_compare_occupancy(), r_config.is_diff_RV0(), r_config.get_all_cent(), is_start_pt, r_config.get_date(), folder, r_config.get_suffix(), r_config.get_energy(), r_config.get_system())
                    analysis.comp_to_other_file(energy_of_other_file, system_of_other_file, "PCM", meson, filename_comp, obj_comp, path_comp, histname_comp, r_config.is_ratio_comp(), config)
run_new()

