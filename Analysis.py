# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3
# This code was written by Alica Enderich (Febuary 2024)
# This code was modified and extended by Julia Schlägel (July 2024)
# This code was modified and extended by Anna Pishchaeva (October 2025)

import os, shutil
import math
import argparse
import pickle
import copy
import numpy as np
import ctypes
import yaml
from copy import deepcopy
import ROOT
import datetime
ROOT.gROOT.SetBatch(True);
from ROOT import TFile, THashList, TF1, TString, TCanvas, TLegend, TPaveText
from Info_decay import Info_decay
from Write_in_file import Write_in_file
from Write_pdf import Write_pdf
from Get_histo import Get_histo
from Calculation import extract_inv_mass_yield, calc_eff_acc, calc_corr_yield
from Make_paths import make_paths
#_________________________________________________________________________________________
class Analysis:
    def __init__(self, is_diff_cent, is_compare_RV0, is_compare_occupancy, is_diff_RV0, all_cent, is_start_pt, date, folder, suffix, energy, system):
        print("#################################################")
        print("Analysis.py")
        self.info_decay = Info_decay(is_diff_cent, is_compare_RV0, is_compare_occupancy, is_diff_RV0, all_cent, is_start_pt, date, suffix, energy, system) # object that contains all of the decay info
        self.w_pdf = Write_pdf(folder, self.info_decay)
        make_paths(folder)
        self.folder = folder
        print("folder = ", self.folder)
        print("info decay object Info_decay.py is initialized")
        print("write pdf object from Write_pdf.py is initialized")
        print("the folders for the outpur file is created {using make_paths(folder)}")
        print("#################################################")


    def inv_mass_raw_yield(self, input_file, config_file, typ, decay, meson):
        """
            _arg_    : for one type of meson, for one type of dateset, for one type of decay
            _summary_: extracts inv mass of the given meson of the given dataset (data/mc) and decay (pcm/dalitz);
                       calculates inv yield; writes the produced histos in the root, pkl, pdf files
           _returns_:  nothing
        """
        self.info_decay.set_typ(typ)
        self.info_decay.set_decay(decay)
        self.info_decay.set_meson(meson)
        self.info_decay.set_config_file(config_file)
        self.info_decay.set_input_file(input_file)
        g_histo = Get_histo(self.info_decay) # object that extracts histos from the file
        self.w_pdf.set_g_histo(g_histo)
        for isys in range(0, self.info_decay.get_number_ss()): #loop over centralities
            if self.info_decay.get_is_compare_occupancy() or self.info_decay.get_is_compare_RV0():
                break
            self.info_decay.set_subsystem(isys)
            self.info_decay.set_subsystem_params(isys)
            ssname = self.info_decay.get_ssname()
            #self.info_decay.set_list_ss()
            self.info_decay.set_centrality(ssname)
            cent1, cent2 = self.info_decay.get_cent()
            if (self.info_decay.get_all_cent() == False) and (cent1 == 0) and (cent2 == 100) and (self.info_decay.get_number_ss()!=1): #if we don't want ot include 0-100%
                continue
            self.info_decay.set_RV0(ssname)
            self.info_decay.set_occupancy(ssname)
            occup1, occup2 = self.info_decay.get_occupancy()

            outname = os.path.join(self.folder,"inv_mass_analysis", "{0}".format(self.info_decay.get_decay()),"{0}_{1}_{2}_{3}_{4}_ptspectrum_{5}_{6}_{7}_cent{8}_{9}_RV0_{10}_{11}_occup_{12}_{13}.root".format(self.info_decay.get_date(), self.info_decay.get_decay(), self.info_decay.get_period(), self.info_decay.get_meson(), self.info_decay.get_typ(), self.info_decay.get_system(), self.info_decay.get_energy(), self.info_decay.get_suffix(), cent1, cent2, self.info_decay.get_V0_radius_min(), self.info_decay.get_V0_radius_max(), occup1, occup2));
            print("Analysis.py, inv_mass_raw_yield, output file name = ", outname, flush=True);

            w_file = Write_in_file(outname) # object to write histos into file

            outlist_inv_mass_histos = extract_inv_mass_yield(self.info_decay, isys) #calculating the invariant mass and the raw yield

            w_file.writing_in_file(self.info_decay.get_decay(), outlist_inv_mass_histos)

            g_histo.set_input_file(outname, self.info_decay)
            self.w_pdf.set_outname(outname)
            self.w_pdf.write_fitted_inv_mass_pdf(self.info_decay, isys, True)  # pdf output of all fitted histograms (inv mass with subtracted mass) for each combination of cut and fit
            if self.info_decay.get_typ() == "data":
                self.w_pdf.write_same_mixed_pdf(self.info_decay, isys, True) # pdf output of all same and mixed scaled histograms (of inv mass) for each cut
                self.w_pdf.write_histo_one_cent_pdf(self.info_decay, "Significance", isys, g_histo.get_significance_list_one_cent(isys)) # pdf output of raw yield for one centrality
                w_file.writing_in_pkl(self.info_decay, g_histo.get_significance_one_cent(isys), "significance")
            self.w_pdf.write_parameter_one_cent_pdf(self.info_decay) # pdf output of parameters of the fit for each combination of cut and fit one cent
            self.w_pdf.write_histo_one_cent_pdf(self.info_decay, "RawYield", isys, g_histo.get_yield_list_one_cent(isys)) # pdf output of raw yield for one centrality
            self.w_pdf.write_histo_one_cent_pdf(self.info_decay, "StatUncert", isys, g_histo.get_stat_uncert_list_one_cent(isys)) # pdf output of raw yield for one centrality
            w_file.writing_in_pkl(self.info_decay, g_histo.get_yield_one_cent(isys), "raw_yield") #writing down the raw yeild into a pkl file (one file - one centrality)
            w_file.writing_in_pkl(self.info_decay, g_histo.get_stat_uncert_one_cent(isys), "stat_uncert")
            if self.info_decay.get_number_ss()!=1:
                g_histo.set_list_parameters_all_cent()
                g_histo.set_list_raw_yield_all_cent(isys)

        if self.info_decay.get_number_ss()!=1:
            if not self.info_decay.get_is_compare_occupancy() and not self.info_decay.get_is_compare_RV0():
                self.w_pdf.write_parameter_all_syst_pdf(self.info_decay, g_histo.get_list_parameters_all_cent(), self.info_decay.get_plotting_sys()) # pdf output of parameters of the fit for each combination of cut and fit all subsystems
                if self.info_decay.get_is_diff_cent():
                    ratio_plot = False
                if self.info_decay.get_is_diff_RV0():
                    ratio_plot = True
                self.w_pdf.write_yield_all_syst_pdf(self.info_decay, "RawYield", g_histo.get_list_raw_yield_all_cent(), self.info_decay.get_plotting_sys(), ratio_plot) # pdf output of raw yield all subsystems
            #self.w_pdf.write_yield_all_cent_Rcp_pdf(self.info_decay, "RawYieldRcp", g_histo.get_list_raw_yield_all_cent())
            if (self.info_decay.get_plotting_sys() == "diff_cent_and_diff_RV0s") or (self.info_decay.get_plotting_sys() == "diff_cent_and_diff_occupancies"):
                for isys in range(0, self.info_decay.get_number_ss()):
                    self.info_decay.set_subsystem(isys)
                    ssname = self.info_decay.get_ssname()
                    self.info_decay.set_centrality(ssname)
                    cent1, cent2 = self.info_decay.get_cent()
                    if (self.info_decay.get_all_cent() == False) and (cent1 == 0) and (cent2 == 100) and (self.info_decay.get_number_ss()!=1): #if we don't want ot include 0-100%
                        continue
                    raw_list_diff_syst = []
                    if self.info_decay.get_is_compare_RV0():
                        RV0_list = self.info_decay.get_RV0_list()
                        length = len(RV0_list)-1
                        occup1, occup2 = self.info_decay.get_occupancy()
                    elif self.info_decay.get_is_compare_occupancy():
                        occup_list = self.info_decay.get_occupancy_list()
                        length = len(occup_list)-1
                        self.info_decay.set_RV0(ssname)
                    for i in range (length):
                        if self.info_decay.get_is_compare_RV0():
                            g_histo.set_input_pkl_file(f"pkl/raw_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_{self.info_decay.get_typ()}_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{RV0_list[i]}_{RV0_list[i+1]}_occup_{occup1}_{occup2}.pkl")
                        elif self.info_decay.get_is_compare_occupancy():
                            g_histo.set_input_pkl_file(f"pkl/raw_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_{self.info_decay.get_typ()}_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup_list[i]}_{occup_list[i+1]}.pkl")
                            print(f"pkl/raw_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_{self.info_decay.get_typ()}_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup_list[i]}_{occup_list[i+1]}.pkl")
                        raw_sys = g_histo.get_histo_from_pkl()
                        raw_list_diff_syst.append(raw_sys)
                    if self.info_decay.get_is_compare_RV0():
                        g_histo.set_input_pkl_file(f"pkl/raw_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_{self.info_decay.get_typ()}_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{RV0_list[0]}_{RV0_list[length]}_occup_{occup1}_{occup2}.pkl")
                    elif self.info_decay.get_is_compare_occupancy():
                         g_histo.set_input_pkl_file(f"pkl/raw_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_{self.info_decay.get_typ()}_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup_list[0]}_{occup_list[length]}.pkl")
                    raw_sys = g_histo.get_histo_from_pkl()
                    raw_list_diff_syst.append(raw_sys)
                    self.w_pdf.write_yield_all_syst_pdf(self.info_decay, "RawYield", raw_list_diff_syst, self.info_decay.get_plotting_sys(), True)

    def eff_acc_corr_yield(self, input_file_mc, input_file_data, config_file_mc, config_file_data, decay, meson):
        """
           _arg_    : for one type of meson, for both types of dateset, for one type of decay
           _summary_: calculates efficiency and acceptance using mc file; calculates corr yield of the data inv yield;
                      writes the produced histos in the root, pkl, pdf files; compares the corr yield and eff_acc of this dateset with ones of another dateset
           _returns_: nothing
        """
        eff_acc_list = []
        corr_list = []

        self.info_decay.set_meson(meson)
        self.info_decay.set_decay(decay)

        info_decay_data = copy.deepcopy(self.info_decay)
        info_decay_data.set_typ("data")
        info_decay_data.set_input_file(input_file_data)
        info_decay_data.set_config_file(config_file_data)

        info_decay_mc = copy.deepcopy(self.info_decay)
        info_decay_mc.set_typ("mc")
        info_decay_mc.set_input_file(input_file_mc)
        info_decay_mc.set_config_file(config_file_mc)

        g_histo = Get_histo(info_decay_data) # object that extracts histos from the file
        self.w_pdf.set_g_histo(g_histo)

        for isys in range(0, info_decay_data.get_number_ss()):
            if self.info_decay.get_is_compare_occupancy() or self.info_decay.get_is_compare_RV0():
                break
            corr_list_one_cent =[]
            corr_scaled_list_one_cent = []
            info_decay_data.set_subsystem(isys)
            ssname_data = info_decay_data.get_ssname()
            info_decay_data.set_centrality(ssname_data)
            info_decay_data.set_RV0(ssname_data)
            info_decay_data.set_occupancy(ssname_data)
            #info_decay_data.set_list_ss()
            info_decay_mc.set_subsystem(isys)
            ssname_mc = info_decay_mc.get_ssname()
            info_decay_mc.set_centrality(ssname_mc)
            info_decay_mc.set_RV0(ssname_mc)
            info_decay_mc.set_occupancy(ssname_mc)
            #info_decay_mc.set_list_ss()
            cent1, cent2 = info_decay_data.get_cent()
            occup1, occup2 = info_decay_data.get_occupancy()
            if (self.info_decay.get_all_cent() == False) and (cent1 == 0) and (cent2 == 100) and (info_decay_data.get_number_ss()!=1): #if we don't want ot include 0-100%
                continue
            filename_rootfile = os.path.join(self.folder,"inv_yield_calculation", "{0}".format(info_decay_data.get_decay()),"{0}_{1}_results_cent_{2}_{3}_RV0_{4}_{5}_occup_{6}_{7}.root".format(self.info_decay.get_date(), ssname_data, cent1, cent2, info_decay_data.get_V0_radius_min(), info_decay_data.get_V0_radius_max(), occup1, occup2))
            w_file = Write_in_file(filename_rootfile) # object to write histos into file
        #extracting and saving raw yeilds of one cetrality on data and mc in rootfiles

            g_histo.set_input_pkl_file(f'pkl/raw_yield_{info_decay_data.get_system()}_{info_decay_data.get_energy()}_{info_decay_data.get_meson()}_data_{info_decay_data.get_decay()}_{cent1}_{cent2}_R_{info_decay_data.get_V0_radius_min()}_{info_decay_data.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl')
            raw_data = g_histo.get_histo_from_pkl()
            w_file.writing_in_one_histo("{0}".format(info_decay_data.get_decay()), raw_data, "h1_yield_data")

            g_histo.set_input_pkl_file(f'pkl/raw_yield_{info_decay_mc.get_system()}_{info_decay_mc.get_energy()}_{info_decay_mc.get_meson()}_mc_{info_decay_mc.get_decay()}_{cent1}_{cent2}_R_{info_decay_mc.get_V0_radius_min()}_{info_decay_mc.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl')
            raw_mc = g_histo.get_histo_from_pkl()
            w_file.writing_in_one_histo("{0}".format(info_decay_data.get_decay()), raw_mc, "h1_yield_mc")

        #calculating the efficiency*acceptance using mc
            eff_acc, eff_acc_scaled, teff_acc = calc_eff_acc(raw_mc, info_decay_mc, w_file, ssname_mc)
            self.w_pdf.write_eff_acc_one_cent_pdf(info_decay_mc, eff_acc, "Efficiency_Acceptance", ssname_mc)
            self.w_pdf.write_eff_acc_one_cent_pdf(info_decay_mc, eff_acc_scaled, "Efficiency_Acceptance_scaled_by_BR", ssname_mc)
            self.w_pdf.write_eff_acc_one_cent_pdf(info_decay_mc, teff_acc, "T_Efficiency_Acceptance_shifted", ssname_mc, True)
            self.w_pdf.write_eff_vs_Teff_pdf(info_decay_mc, eff_acc_scaled, teff_acc, isys)
            w_file.writing_in_pkl(info_decay_mc, eff_acc, "eff_acc")
            w_file.writing_in_pkl(info_decay_mc, eff_acc_scaled, "eff_acc_scaled")
            w_file.writing_in_one_histo("{0}".format(info_decay_mc.get_decay()), eff_acc, "h_eff_acc")
            w_file.writing_in_one_histo("{0}".format(info_decay_mc.get_decay()), eff_acc_scaled, "h_eff_acc_scaled")
            eff_acc_list.append(copy.deepcopy(eff_acc))

        #calculating the corrected yield of data
            pyth, corr, corr_scaled = calc_corr_yield(raw_data, eff_acc, info_decay_mc, info_decay_data)
            corr_list_one_cent.append(corr)
            corr_scaled_list_one_cent.append(corr_scaled)
            self.w_pdf.write_histo_one_cent_pdf(info_decay_data, "CorrectedYield", isys, corr_list_one_cent)
            self.w_pdf.write_histo_one_cent_pdf(info_decay_data, "CorrectedYield", isys, corr_scaled_list_one_cent, True)
            self.w_pdf.write_corr_yield_fit_fun_pdf(corr_scaled, info_decay_data, "this_data")
            self.w_pdf.write_pyth_one_cent_pdf(pyth, info_decay_mc)
            self.w_pdf.write_pyth_vs_corr_pdf(pyth, corr, info_decay_mc)
            w_file.writing_in_one_histo("{0}".format(info_decay_data.get_decay()), pyth)
            w_file.writing_in_one_histo("{0}".format(info_decay_data.get_decay()), corr)
            w_file.writing_in_pkl(info_decay_data, corr, "corr_yield")
            corr_list.append(copy.deepcopy(corr))

        #comapring data parameters and stat. uncert. of the fit and mc parameters of the fit
            outname_param_data = os.path.join(self.folder,"inv_mass_analysis", "{0}".format(info_decay_data.get_decay()),"{0}_{1}_{2}_{3}_{4}_ptspectrum_{5}_{6}_{7}_cent{8}_{9}_RV0_{10}_{11}_occup_{12}_{13}.root".format(info_decay_data.get_date(), info_decay_data.get_decay(), info_decay_data.get_period(), info_decay_data.get_meson(), info_decay_data.get_typ(), info_decay_data.get_system(), info_decay_data.get_energy(), info_decay_data.get_suffix(), cent1, cent2, info_decay_data.get_V0_radius_min(), info_decay_data.get_V0_radius_max(), occup1, occup2))
            outname_param_mc = os.path.join(self.folder,"inv_mass_analysis", "{0}".format(info_decay_mc.get_decay()),"{0}_{1}_{2}_{3}_{4}_ptspectrum_{5}_{6}_{7}_cent{8}_{9}_RV0_{10}_{11}_occup_{12}_{13}.root".format(info_decay_mc.get_date(), info_decay_mc.get_decay(), info_decay_mc.get_period(), info_decay_mc.get_meson(), info_decay_mc.get_typ(), info_decay_mc.get_system(), info_decay_mc.get_energy(), info_decay_mc.get_suffix(), cent1, cent2, info_decay_data.get_V0_radius_min(), info_decay_data.get_V0_radius_max(), occup1, occup2))
            g_histo_data = Get_histo(info_decay_data) # object that extracts histos from the file
            g_histo_data.set_input_file(outname_param_data, info_decay_data)
            g_histo_mc = Get_histo(info_decay_mc) # object that extracts histos from the file
            g_histo_mc.set_input_file(outname_param_mc, info_decay_mc)
            parameter_list_all = []
            parameter_list_data = g_histo_data.get_parameter_list()
            parameter_list_mc = g_histo_mc.get_parameter_list()
            parameter_list_all.append(parameter_list_data)
            parameter_list_all.append(parameter_list_mc) #parameter_list_all = [data, mc]
            self.w_pdf.write_histo_one_vs_another_and_ratio_pdf(info_decay_data, parameter_list_data[1], parameter_list_mc[1], "mean", "mc_data", info_decay_mc.get_period())
            self.w_pdf.write_histo_one_vs_another_and_ratio_pdf(info_decay_data, parameter_list_data[3], parameter_list_mc[3], "sigma", "mc_data", info_decay_mc.get_period())
            self.w_pdf.write_histo_one_vs_another_and_ratio_pdf(info_decay_data, parameter_list_data[2], parameter_list_mc[2], "lambda", "mc_data", info_decay_mc.get_period())
            self.w_pdf.write_parameter_all_syst_pdf(info_decay_data, parameter_list_all, "mc_and_data") # pdf output of parameters of the fit for each combination of cut and fit all cent
            stat_uncert_mc_data_list = []
            stat_uncert_mc_data_list.append(g_histo_data.get_stat_uncert_one_cent(isys))
            stat_uncert_mc_data_list.append(g_histo_mc.get_stat_uncert_one_cent(isys))
            self.w_pdf.write_yield_all_syst_pdf(info_decay_data, "StatUncert", stat_uncert_mc_data_list, info_decay_data.get_plotting_sys(), False, False, True, info_decay_mc.get_period())

        if info_decay_data.get_number_ss()!=1:
            if not self.info_decay.get_is_compare_occupancy() and not self.info_decay.get_is_compare_RV0():
                if info_decay_data.get_is_diff_cent():
                    ratio_plot = False
                if info_decay_data.get_is_diff_RV0():
                    ratio_plot = True
                self.w_pdf.write_yield_all_syst_pdf(info_decay_data, "CorrectedYield", corr_list, info_decay_data.get_plotting_sys(), ratio_plot)
                #self.w_pdf.write_yield_all_cent_Rcp_pdf(info_decay_data, "CorrectedYieldRcp", corr_list)
                self.w_pdf.write_eff_acc_all_syst_pdf(info_decay_mc, eff_acc_list, "Efficiency_Acceptance", info_decay_mc.get_plotting_sys(), ratio_plot)
            if (self.info_decay.get_plotting_sys() == "diff_cent_and_diff_RV0s") or (self.info_decay.get_plotting_sys() == "diff_cent_and_diff_occupancies"):
                for isys in range(0, info_decay_data.get_number_ss()):
                    corr_list_diff_sys = []
                    eff_acc_list_diff_sys = []
                    info_decay_data.set_subsystem(isys)
                    ssname = info_decay_data.get_ssname()
                    info_decay_data.set_centrality(ssname)
                    info_decay_mc.set_subsystem(isys)
                    ssname_mc = info_decay_data.get_ssname()
                    info_decay_mc.set_centrality(ssname_mc)
                    cent1, cent2 = info_decay_data.get_cent()
                    if (self.info_decay.get_all_cent() == False) and (cent1 == 0) and (cent2 == 100) and (info_decay_data.get_number_ss()!=1): #if we don't want ot include 0-100%
                        continue
                    if self.info_decay.get_is_compare_RV0():
                        RV0_list = info_decay_data.get_RV0_list()
                        length = len(RV0_list)-1
                        occup1, occup2 = info_decay_data.get_occupancy()
                    if self.info_decay.get_is_compare_occupancy():
                        occup_list = info_decay_data.get_occupancy_list()
                        length = len(occup_list)-1
                        info_decay_data.set_RV0(ssname)
                        info_decay_mc.set_RV0(ssname_mc)
                    for i in range (length):
                        if self.info_decay.get_is_compare_RV0():
                            g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{RV0_list[i]}_{RV0_list[i+1]}_occup_{occup1}_{occup2}.pkl")
                        elif self.info_decay.get_is_compare_occupancy():
                            g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{info_decay_data.get_V0_radius_min()}_{info_decay_data.get_V0_radius_max()}_occup_{occup_list[i]}_{occup_list[i+1]}.pkl")
                        corr_sys = g_histo.get_histo_from_pkl()
                        corr_list_diff_sys.append(corr_sys)
                        if self.info_decay.get_is_compare_RV0():
                            g_histo.set_input_pkl_file(f"pkl/eff_acc_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_mc_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{RV0_list[i]}_{RV0_list[i+1]}_occup_{occup1}_{occup2}.pkl")
                        elif self.info_decay.get_is_compare_occupancy():
                            g_histo.set_input_pkl_file(f"pkl/eff_acc_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_mc_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{info_decay_data.get_V0_radius_min()}_{info_decay_data.get_V0_radius_max()}_occup_{occup_list[i]}_{occup_list[i+1]}.pkl")
                        eff_acc_sys = g_histo.get_histo_from_pkl()
                        eff_acc_list_diff_sys.append(eff_acc_sys)
                    if self.info_decay.get_is_compare_RV0():
                            g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{RV0_list[0]}_{RV0_list[length]}_occup_{occup1}_{occup2}.pkl")
                    elif self.info_decay.get_is_compare_occupancy():
                            g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{info_decay_data.get_V0_radius_min()}_{info_decay_data.get_V0_radius_max()}_occup_{occup_list[0]}_{occup_list[length]}.pkl")
                    corr_sys = g_histo.get_histo_from_pkl()
                    corr_list_diff_sys.append(corr_sys)
                    if self.info_decay.get_is_compare_RV0():
                            g_histo.set_input_pkl_file(f"pkl/eff_acc_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_mc_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{RV0_list[0]}_{RV0_list[length]}_occup_{occup1}_{occup2}.pkl")
                    elif self.info_decay.get_is_compare_occupancy():
                            g_histo.set_input_pkl_file(f"pkl/eff_acc_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_mc_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{info_decay_data.get_V0_radius_min()}_{info_decay_data.get_V0_radius_max()}_occup_{occup_list[0]}_{occup_list[length]}.pkl")
                    eff_acc_sys = g_histo.get_histo_from_pkl()
                    eff_acc_list_diff_sys.append(eff_acc_sys)
                    self.w_pdf.write_yield_all_syst_pdf(info_decay_data, "CorrectedYield", corr_list_diff_sys, self.info_decay.get_plotting_sys(), True)
                    self.w_pdf.write_eff_acc_all_syst_pdf(info_decay_mc, eff_acc_list_diff_sys, "Efficiency_Acceptance", self.info_decay.get_plotting_sys(), True)

    def comp_dalitz_pcm(self, config_file_pcm, config_file_dalitz, meson):
        """
           _arg_    : for one type of meson, for both types of dateset, for both types of decay
           _summary_: compares pcm eff_acc, corr yield with dalitz eff_acc, corr yield; writes the produced histos in the root, pkl, pdf files;
           _returns_: nothing
        """
        self.info_decay.set_typ("data")
        self.info_decay.set_decay("PCM")
        self.info_decay.set_meson(meson)
        self.info_decay.set_config_file(config_file_pcm)

        g_histo = Get_histo(self.info_decay) # object that extracts histos from the file
        self.w_pdf.set_g_histo(g_histo)

        for isys in range(0, self.info_decay.get_number_ss()):
            self.info_decay.set_subsystem(isys)
            ssname = self.info_decay.get_ssname()
            self.info_decay.set_centrality(ssname)
            cent1, cent2 = self.info_decay.get_cent()
            self.info_decay.set_RV0(ssname)
            self.info_decay.set_occupancy(ssname)
            occup1, occup2 = self.info_decay.get_occupancy()
            if (self.info_decay.get_all_cent() == False) and (cent1 == 0) and (cent2 == 100) and (self.info_decay.get_number_ss()!=1): #if we don't want ot include 0-100%
                continue

            g_histo.set_input_pkl_file(f"pkl/eff_acc_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_mc_PCM_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
            eff_acc_pcm = g_histo.get_histo_from_pkl()
            g_histo.set_input_pkl_file(f"pkl/eff_acc_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_mc_Dalitz_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
            eff_acc_dalitz = g_histo.get_histo_from_pkl()
            self.w_pdf.write_histo_one_vs_another_pdf(self.info_decay, eff_acc_pcm, eff_acc_dalitz, "eff_acc", "pcm_dalitz")

            g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_data_PCM_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
            corr_pcm = g_histo.get_histo_from_pkl()
            g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_data_Dalitz_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
            corr_dalitz = g_histo.get_histo_from_pkl()
            self.w_pdf.write_histo_one_vs_another_pdf(self.info_decay, corr_pcm, corr_dalitz, "corr_yield", "pcm_dalitz")


    def comp_eta_pi0(self, config_file_pi0, config_file_eta, decay):
        """
           _arg_    : for both types of meson, for both types of dateset, for one type of decay
           _summary_: compares pi0 eff_acc, corr yield with eta eff_acc, corr yield;
                      compares the eta/pi0 ratio of corr yield of this dateset with ones of another dateset
                      writes the produced histos in the root, pkl, pdf files;
           _returns_: nothing
        """

        self.info_decay.set_typ("data")
        self.info_decay.set_decay(decay)
        self.info_decay.set_meson("pi0")
        self.info_decay.set_config_file(config_file_pi0)

        g_histo = Get_histo(self.info_decay) # object that extracts histos from the file
        self.w_pdf.set_g_histo(g_histo)

        for isys in range(0, self.info_decay.get_number_ss()):
            self.info_decay.set_subsystem(isys)
            ssname = self.info_decay.get_ssname()
            self.info_decay.set_centrality(ssname)
            cent1, cent2 = self.info_decay.get_cent()
            self.info_decay.set_RV0(ssname)
            self.info_decay.set_occupancy(ssname)
            occup1, occup2 = self.info_decay.get_occupancy()
            if (self.info_decay.get_all_cent() == False) and (cent1 == 0) and (cent2 == 100) and (self.info_decay.get_number_ss()!=1): #if we don't want ot include 0-100%
                continue

            g_histo.set_input_pkl_file(f"pkl/eff_acc_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_pi0_mc_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
            eff_acc_pi0 = g_histo.get_histo_from_pkl()
            g_histo.set_input_pkl_file(f"pkl/eff_acc_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_eta_mc_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
            eff_acc_eta = g_histo.get_histo_from_pkl()
            self.w_pdf.write_histo_one_vs_another_pdf(self.info_decay, eff_acc_pi0, eff_acc_eta, "eff_acc", "pi0_eta")
            self.w_pdf.write_histo_one_vs_another_ratio_pdf(self.info_decay, eff_acc_pi0, eff_acc_eta, "eff_acc", "pi0_eta")

            g_histo.set_input_pkl_file(f"pkl/eff_acc_scaled_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_pi0_mc_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
            eff_acc_scaled_pi0 = g_histo.get_histo_from_pkl()
            g_histo.set_input_pkl_file(f"pkl/eff_acc_scaled_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_eta_mc_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
            eff_acc_scaled_eta = g_histo.get_histo_from_pkl()
            self.w_pdf.write_histo_one_vs_another_pdf(self.info_decay, eff_acc_scaled_pi0, eff_acc_scaled_eta, "eff_acc_scaled", "pi0_eta")
            self.w_pdf.write_histo_one_vs_another_ratio_pdf(self.info_decay, eff_acc_scaled_pi0, eff_acc_scaled_eta, "eff_acc_scaled", "pi0_eta")

            g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_pi0_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
            corr_pi0 = g_histo.get_histo_from_pkl()
            g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_eta_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
            corr_eta = g_histo.get_histo_from_pkl()
            self.w_pdf.write_histo_one_vs_another_pdf(self.info_decay, corr_pi0, corr_eta, "corr_yield", "pi0_eta")
            self.w_pdf.write_histo_one_vs_another_ratio_pdf(self.info_decay, corr_pi0, corr_eta, "corr_yield", "pi0_eta")

    def comp_to_other_file(self, energy_of_other_file, system_of_other_file, decay, meson, filename_comp, obj_comp, path_comp, histname_comp, is_ratio_comp, config):

        self.info_decay.set_typ("data")
        self.info_decay.set_decay(decay)
        self.info_decay.set_meson(meson)
        self.info_decay.set_config_file(config)
        corr_yield_from_other_file_list = []
        err_stat_list = []
        err_syst_list = []

        for icentrality in range(0, self.info_decay.get_number_ss()-1):
            self.info_decay.set_subsystem(icentrality)
            ssname = self.info_decay.get_ssname()
            self.info_decay.set_centrality(ssname)
            self.info_decay.set_RV0(ssname)
            self.info_decay.set_occupancy(ssname)
            cent1, cent2 = self.info_decay.get_cent()
            occup1, occup2 = self.info_decay.get_occupancy()

            g_histo = Get_histo(self.info_decay)
            self.w_pdf.set_g_histo(g_histo)
            g_histo.set_input_file(filename_comp, self.info_decay)

            if obj_comp == "corr_yield":
                corr_yield_list = []
                g_histo.set_input_pkl_file(f'pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl')
                corr = g_histo.get_histo_from_pkl()
                corr_scaled_pi = corr.Clone()
                corr_scaled_pi.Scale(2*np.pi)
                corr_yield_list.append(corr_scaled_pi)
                hist_notrebin, corr_yield_from_other_file, h_err_stat, h_err_syst = g_histo.get_corr_yield_from_other_file(path_comp, icentrality, self.info_decay.get_bin_var()) #0 - 0-20%, 1 - 20-40%
                corr_yield_from_other_file_list.append(copy.deepcopy(corr_yield_from_other_file))
                err_stat_list.append(copy.deepcopy(h_err_stat))
                err_syst_list.append(copy.deepcopy(h_err_syst))
                output_file_comp = ROOT.TFile("test.root", "RECREATE")
                corr_yield_from_other_file.Write()
                corr_yield_list.append(corr_yield_from_other_file)
                self.w_pdf.write_comp_with_other_data_pdf(self.info_decay, corr_yield_list, obj_comp)
                self.w_pdf.write_histo_one_vs_another_ratio_pdf(self.info_decay, corr_scaled_pi, corr_yield_from_other_file, "corr_yield", "other_data")

            if obj_comp == "eff_acc":
                eff_acc_list = []
                g_histo.set_input_pkl_file(f"pkl/eff_acc_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_mc_PCM_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
                eff_acc = g_histo.get_histo_from_pkl()
                eff_acc_list.append(eff_acc)
                eff_acc_from_other_file = g_histo.get_eff_acc_BR_from_other_file(histname_comp, path_comp)
                eff_acc_list.append(eff_acc_from_other_file)
                self.w_pdf.write_comp_with_other_data_pdf(self.info_decay, eff_acc_list, obj_comp)

            if obj_comp == "eta_pi0_ratio":
                ratio_list = []
                g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_pi0_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
                corr_pi0 = g_histo.get_histo_from_pkl()
                g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_eta_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
                corr_eta = g_histo.get_histo_from_pkl()
                ratio = corr_eta.Clone("ratio_eta_pi0")
                ratio.Sumw2()
                ratio.Divide(corr_pi0)
                ratio_list.append(ratio)
                ratio_eta_pi0_other_file = g_histo.get_histo_from_other_file(histname_comp, path_comp)
                ratio_list.append(ratio_eta_pi0_other_file)
                self.w_pdf.write_comp_with_other_data_pdf(self.info_decay, ratio_list, obj_comp)

            if obj_comp == "RV0s_of_diff_systems":
                corr_list = []
                g_histo.set_input_pkl_file(f"pkl/corr_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
                corr_PbPb = g_histo.get_histo_from_pkl()
                corr_list.append(corr_PbPb)
                g_histo.set_input_pkl_file(f"pkl/corr_yield_O-O_5.36 TeV_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_0_100_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}.pkl")
                corr_OO = g_histo.get_histo_from_pkl()
                corr_list.append(corr_OO)
                g_histo.set_input_pkl_file(f"pkl/corr_yield_p-O_9.61 TeV_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_0_100_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}.pkl")
                corr_pO = g_histo.get_histo_from_pkl()
                corr_list.append(corr_pO)
                self.w_pdf.write_comp_with_other_data_pdf(self.info_decay, corr_list, "corr_diff_pkl")

                eff_acc_list = []
                g_histo.set_input_pkl_file(f"pkl/eff_acc_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_mc_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
                eff_acc_PbPb = g_histo.get_histo_from_pkl()
                eff_acc_list.append(eff_acc_PbPb)
                g_histo.set_input_pkl_file(f"pkl/eff_acc_O-O_5.36 TeV_{self.info_decay.get_meson()}_mc_{self.info_decay.get_decay()}_0_100_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}.pkl")
                eff_acc_OO = g_histo.get_histo_from_pkl()
                eff_acc_list.append(eff_acc_OO)
                g_histo.set_input_pkl_file(f"pkl/eff_acc_p-O_9.61 TeV_{self.info_decay.get_meson()}_mc_{self.info_decay.get_decay()}_0_100_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}.pkl")
                eff_acc_pO = g_histo.get_histo_from_pkl()
                eff_acc_list.append(eff_acc_pO)
                self.w_pdf.write_comp_with_other_data_pdf(self.info_decay, eff_acc_list, "eff_acc_diff_pkl")

                raw_list_data = []
                g_histo.set_input_pkl_file(f"pkl/raw_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
                raw_PbPb_data = g_histo.get_histo_from_pkl()
                raw_list_data.append(raw_PbPb_data)
                g_histo.set_input_pkl_file(f"pkl/raw_yield_O-O_5.36 TeV_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_0_100_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}.pkl")
                raw_OO_data = g_histo.get_histo_from_pkl()
                raw_list_data.append(raw_OO_data)
                g_histo.set_input_pkl_file(f"pkl/raw_yield_p-O_9.61 TeV_{self.info_decay.get_meson()}_data_{self.info_decay.get_decay()}_0_100_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}.pkl")
                raw_pO_data = g_histo.get_histo_from_pkl()
                raw_list_data.append(raw_pO_data)
                self.w_pdf.write_comp_with_other_data_pdf(self.info_decay, raw_list_data, "raw_diff_data_pkl")

                raw_list_mc = []
                g_histo.set_input_pkl_file(f"pkl/raw_yield_{self.info_decay.get_system()}_{self.info_decay.get_energy()}_{self.info_decay.get_meson()}_mc_{self.info_decay.get_decay()}_{cent1}_{cent2}_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl")
                raw_PbPb_mc = g_histo.get_histo_from_pkl()
                raw_list_mc.append(raw_PbPb_mc)
                g_histo.set_input_pkl_file(f"pkl/raw_yield_O-O_5.36 TeV_{self.info_decay.get_meson()}_mc_{self.info_decay.get_decay()}_0_100_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}.pkl")
                raw_OO_mc = g_histo.get_histo_from_pkl()
                raw_list_mc.append(raw_OO_mc)
                g_histo.set_input_pkl_file(f"pkl/raw_yield_p-O_9.61 TeV_{self.info_decay.get_meson()}_mc_{self.info_decay.get_decay()}_0_100_R_{self.info_decay.get_V0_radius_min()}_{self.info_decay.get_V0_radius_max()}.pkl")
                raw_pO_mc = g_histo.get_histo_from_pkl()
                raw_list_mc.append(raw_pO_mc)
                self.w_pdf.write_comp_with_other_data_pdf(self.info_decay, raw_list_mc, "raw_diff_mc_pkl")
        #if obj_comp == "corr_yield" and Rcp_other_data == True:
        #self.w_pdf.write_yield_all_cent_Rcp_pdf(self.info_decay, "CorrectedYieldRcp", corr_yield_from_other_file_list, True, energy_of_other_file, system_of_other_file, err_stat_list, err_syst_list)



