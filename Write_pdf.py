# This code was written by Julia Schlägel (July 2024)
import os, sys, shutil
import ROOT
import yaml
from ctypes import *
import numpy as np
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency, TString
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX

from Get_histo import Get_histo
from Plot_inv_mass import Plot_inv_mass
from Plot_yield import Plot_yield
from Plot_efficiency import Plot_efficiency
from Compare import Compare
# from Analysis import Analysis
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

class Write_pdf: #in this class all central things are done that are needed in multiple macros like loading a file from a rootfile,...
    def __init__(self, folder, info_decay):
        print("/////////////////////////////////////////////////")
        print("Write_pdf.py")
        self.folder = folder
        print("folder =", self.folder)
        print("/////////////////////////////////////////////////")

    def set_outname(self, outname):
        self.outname = outname

    def set_g_histo(self, g_histo):
        self.g_histo = g_histo

    def write_fitted_inv_mass_pdf(self, info_decay, isys, is_one_pad):
        plot_fit_inv_mass = Plot_inv_mass(info_decay.get_meson(), info_decay.get_ssname()) #this object will write down everything in the output file + in the down lines it will also draw histos for inv mass
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        subtracted_list = self.g_histo.get_subtracted_list(isys, info_decay.get_scale_mix_range())
        function_list = self.g_histo.get_function_list(isys, info_decay.get_fit_func())
        output_name_with_yield = os.path.join(self.folder,"inv_mass_analysis", "{0}".format(info_decay.get_decay()),"fitted_inv_mass" ,"{0}_{1}_{2}_InvMass_Overview_with_yield_cent{3}_{4}_npt_{5}_R_{6}_{7}_occup_{8}_{9}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), cent1, cent2, len(info_decay.get_arr_pt()), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        plot_fit_inv_mass.plot_inv_mass("inv_mass_subtracted", info_decay, subtracted_list, function_list, TString(output_name_with_yield), "yield")
        if is_one_pad:
            output_name_with_yield = os.path.join(self.folder,"inv_mass_analysis", "{0}".format(info_decay.get_decay()),"fitted_inv_mass", "{0}_{1}_{2}_InvMass_Overview_with_yield_cent{3}_{4}_one_pad_npt_{5}_R_{6}_{7}_occup_{8}_{9}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), cent1, cent2, len(info_decay.get_arr_pt()), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
            plot_fit_inv_mass.plot_inv_mass("inv_mass_subtracted", info_decay, subtracted_list, function_list, TString(output_name_with_yield), "yield", is_one_pad)
        del plot_fit_inv_mass

    def write_same_mixed_pdf(self, info_decay, isys, is_one_pad):
        plot_same_mixed = Plot_inv_mass(info_decay.get_meson(), info_decay.get_ssname())
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        same_list = self.g_histo.get_same_list(isys)
        mixed_list = self.g_histo.get_mixed_list(isys, info_decay.get_scale_mix_range())
        output_name = os.path.join(self.folder, "inv_mass_analysis","{0}".format(info_decay.get_decay()), "same_mixed_inv_mass", "{0}_{1}_{2}_InvMass_Scaled_Mixed_with_yield_cent{3}_{4}_npt_{5}_R_{6}_{7}_occup_{8}_{9}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), cent1, cent2, len(info_decay.get_arr_pt()), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        plot_same_mixed.plot_inv_mass("inv_mass_same_mix", info_decay, same_list, mixed_list, TString(output_name), "no yield") #for pi0 binning: 5, 6 (for eta: 3,4)
        if is_one_pad:
            output_name = os.path.join(self.folder, "inv_mass_analysis","{0}".format(info_decay.get_decay()), "same_mixed_inv_mass", "{0}_{1}_{2}_InvMass_Scaled_Mixed_with_yield_cent{3}_{4}_one_pad_{5}_R_{6}_{7}_occup_{8}_{9}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), cent1, cent2, len(info_decay.get_arr_pt()), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
            plot_same_mixed.plot_inv_mass("inv_mass_same_mix", info_decay, same_list, mixed_list, TString(output_name), "no yield", is_one_pad) #for pi0 binning: 5, 6 (for eta: 3,4)
        del plot_same_mixed

    def write_parameter_one_cent_pdf(self, info_decay):
        plot_parameter_fit = Plot_inv_mass(info_decay.get_meson(), info_decay.get_ssname())
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        parameter_list = self.g_histo.get_parameter_list()
        outname_histo = os.path.join(self.folder, "inv_mass_analysis","{0}".format(info_decay.get_decay()),"parameters","{0}_{1}_{2}_InvMass_Fitparameters_cent{3}_{4}_npt{5}_R_{6}_{7}_occup_{8}_{9}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), cent1, cent2, len(info_decay.get_arr_pt()), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        plot_parameter_fit.plot_parameters(info_decay, parameter_list, TString(outname_histo))
        del plot_parameter_fit

    def write_parameter_all_syst_pdf(self, info_decay, list_parameters_all_cent, plotting):
        ssname_list = info_decay.get_ssname_list()
        plot_parameters_all_cent = Plot_inv_mass(info_decay.get_meson(), info_decay.get_ssname(), plotting)
        if plotting == "centralities":
            occup1, occup2 = info_decay.get_occupancy()
            outname_histo_linear = os.path.join(self.folder, "inv_mass_analysis", "{0}".format(info_decay.get_decay()), "parameters", "{0}_{1}_{2}_InvMass_Parameters_Combined_all_cent_step{3}_npt{4}_R_{5}_{6}_occup_{7}_{8}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), info_decay.get_step_cent(), len(info_decay.get_arr_pt()), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2));
        if plotting == "RV0s" or plotting == "diff_cent_and_diff_RV0s":
            cent1, cent2 = info_decay.get_cent()
            occup1, occup2 = info_decay.get_occupancy()
            outname_histo_linear = os.path.join(self.folder, "inv_mass_analysis", "{0}".format(info_decay.get_decay()), "parameters", "{0}_{1}_{2}_InvMass_Parameters_Combined_all_RV0_cent{3}_{4}_npt{5}_occup_{6}_{7}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), info_decay.get_step_cent(), len(info_decay.get_arr_pt()), cent1, cent2, occup1, occup2))
        if plotting == "mc_and_data":
            cent1, cent2 = info_decay.get_cent()
            occup1, occup2 = info_decay.get_occupancy()
            outname_histo_linear = os.path.join(self.folder, "inv_mass_analysis", "{0}".format(info_decay.get_decay()), "parameters", "{0}_{1}_{2}_InvMass_Parameters_mc_and_data_cent_{3}_{4}_npt{5}_R_{6}_{7}_occup_{8}_{9}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), cent1, cent2, len(info_decay.get_arr_pt()), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2));
        if plotting == "occupancies" or plotting == "diff_cent_and_diff_occupancies":
            cent1, cent2 = info_decay.get_cent()
            outname_histo_linear = os.path.join(self.folder, "inv_mass_analysis", "{0}".format(info_decay.get_decay()), "parameters", "{0}_{1}_{2}_InvMass_Parameters_all_occupancies_cent_{3}_{4}_npt{5}_R_{6}_{7}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), cent1, cent2, len(info_decay.get_arr_pt()), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max()));
        plot_parameters_all_cent.plot_parameters(info_decay, list_parameters_all_cent, TString(outname_histo_linear), ssname_list, True)
        del plot_parameters_all_cent

    def write_histo_one_cent_pdf(self, info_decay, type_of_histo, isys, list_one_cent, scale_pT = False):
        plot_histo_one_cent = Plot_yield(info_decay.get_decay(), info_decay.get_meson())
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        if type_of_histo !="CorrectedYield":
            plot_name = os.path.join(self.folder, "inv_mass_analysis", "{0}".format(info_decay.get_decay()),"raw_yield", "{0}_{1}_{2}_InvMass_{3}_cent{4}_{5}_R_{6}_{7}_occup_{8}_{9}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), type_of_histo, cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        else:
            if scale_pT:
                plot_name = os.path.join(self.folder, "inv_yield_calculation", "{0}".format(info_decay.get_decay()),"corr_yield", "{0}_corrected_yield_scaled_by_pT_cent{1}_{2}_R_{3}_{4}_occup_{5}_{6}.pdf".format(info_decay.get_date(), cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
            else:
                plot_name = os.path.join(self.folder, "inv_yield_calculation", "{0}".format(info_decay.get_decay()),"corr_yield", "{0}_corrected_yield_cent{1}_{2}_R_{3}_{4}_occup_{5}_{6}.pdf".format(info_decay.get_date(), cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))

        plot_histo_one_cent.plot_histo(type_of_histo, list_one_cent, TString(plot_name), info_decay, info_decay.get_plotting_sys(), False, scale_pT)
        del plot_histo_one_cent

    def write_corr_yield_fit_fun_pdf(self, corr_scaled, info_decay, which_data):
        plot_corr_fit = Plot_yield(info_decay.get_decay(), info_decay.get_meson())
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        plot_name = os.path.join(self.folder, "inv_yield_calculation", "{0}".format(info_decay.get_decay()),"corr_yield", "TCM_fit_{0}_cent{1}_{2}_R_{3}_{4}_occup_{5}_{6}.pdf".format(which_data, cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        if which_data == "this_data":
            energy = info_decay.get_energy()
        else:
            energy = "5.02 TeV"
        plot_corr_fit.plot_fit_func(corr_scaled, energy, info_decay, which_data, TString(plot_name))
        del plot_corr_fit

    def write_comp_with_other_data_pdf(self, info_decay, histo_list, name):
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        plot_comp_corr = Compare(info_decay.get_decay(), info_decay.get_meson(), cent1, cent2)
        if name == "inv_cross_section":
            info_decay.scale_to_inv_yield(histo_list[0])
        plot_name = os.path.join(self.folder, "Comparison", "{0}_Comp_{1}_{2}_cent{3}_{4}_R_{5}_{6}_occup_{7}_{8}.pdf".format(info_decay.get_date(), name, info_decay.get_decay(), cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        plot_comp_corr.comp_with_other_data(histo_list, name, TString(plot_name), info_decay.get_energy(), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max())
        del plot_comp_corr

    def write_histo_one_vs_another_pdf(self, info_decay, histo, histo_comp, name, comparison):
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        plot_histo_one_vs_another = Compare(info_decay.get_decay(), info_decay.get_meson(), cent1, cent2)
        plot_name = os.path.join(self.folder, "Comparison", "comp_{0}_{1}_cent{2}_{3}_R_{4}_{5}_occup_{6}_{7}.pdf".format(name, comparison, cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        plot_histo_one_vs_another.plot_histo_one_vs_another(histo, histo_comp, TString(plot_name), info_decay.get_system(), info_decay.get_energy(), name, comparison, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max())
        del plot_histo_one_vs_another

    def write_histo_one_vs_another_ratio_pdf(self, info_decay, histo, histo_comp, name, comparison):
        ratio = histo_comp.Clone("Ratio_corr")
        ratio.Sumw2()
        ratio.Divide(histo)
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        plot_histo_one_vs_another = Compare(info_decay.get_decay(), info_decay.get_meson(), cent1, cent2)
        plot_name = os.path.join(self.folder, "Comparison", "comp_ratio_{0}_{1}_cent{2}_{3}_R_{4}_{5}_occup_{6}_{7}.pdf".format(name, comparison, cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        plot_histo_one_vs_another.plot_ratio(ratio, TString(plot_name), info_decay.get_system(), info_decay.get_energy(), name, comparison, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max())
        del plot_histo_one_vs_another

    def write_histo_one_vs_another_and_ratio_pdf(self, info_decay, histo, histo_comp, name, comparison, period_mc = 0):
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        plot_histo_one_vs_another_w_ratio = Compare(info_decay.get_decay(), info_decay.get_meson(), cent1, cent2)
        plot_name = os.path.join(self.folder, "Comparison", "comp_{0}_{1}_cent{2}_{3}_R_{4}_{5}_occup_{6}_{7}.pdf".format(name, comparison, cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        if comparison == "mc_data":
            plot_histo_one_vs_another_w_ratio.plot_histo_one_vs_another_w_ratio(histo, histo_comp, TString(plot_name), info_decay.get_energy(), name, comparison, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), info_decay.get_system(), info_decay.get_period(), period_mc)
        else:
            plot_histo_one_vs_another_w_ratio.plot_histo_one_vs_another_w_ratio(histo, histo_comp, TString(plot_name), info_decay.get_energy(), name, comparison, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), info_decay.get_system(), info_decay.get_period())
        del plot_histo_one_vs_another_w_ratio

    def write_pyth_one_cent_pdf(self, histo, info_decay):
        plot_pyth_one_cent = Plot_yield(info_decay.get_decay(), info_decay.get_meson())
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        pythia_name = os.path.join(self.folder, "inv_yield_calculation", "{0}".format(info_decay.get_decay()), "corr_yield", "{0}_pythia_cent{1}_{2}_R_{3}_{4}_occup_{5}_{6}.pdf".format(info_decay.get_date(), cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        plot_pyth_one_cent.plot_pythia(histo, info_decay, TString(pythia_name))
        del plot_pyth_one_cent

    def write_pyth_vs_corr_pdf(self, pyth, corr, info_decay):
        plot_pyth_vs_corr = Plot_yield(info_decay.get_decay(), info_decay.get_meson())
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        plot_name = os.path.join(self.folder, "inv_yield_calculation", "{0}".format(info_decay.get_decay()), "corr_yield", "{0}_comp_pythia_corr_yield_cent{1}_{2}_R_{3}_{4}_occup_{5}_{6}.pdf".format(info_decay.get_date(), cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        plot_pyth_vs_corr.plot_pyth_vs_corr(pyth, corr, info_decay, TString(plot_name), False)
        del plot_pyth_vs_corr

    def write_yield_all_syst_pdf(self, info_decay, type_of_histo, yield_list, plotting, ratio = False, scale_pT = False, is_mc_and_data = False, period_mc = 0):
        plot_yield = Plot_yield(info_decay.get_decay(), info_decay.get_meson())
        if type_of_histo != "CorrectedYield":
            folder2 = "inv_mass_analysis"
            folder3 = "raw_yield"
        else:
            folder2 = "inv_yield_calculation"
            folder3 = "corr_yield"
        if plotting == "centralities":
            V0_radius_min = info_decay.get_V0_radius_min()
            V0_radius_max = info_decay.get_V0_radius_max()
            occup1, occup2 = info_decay.get_occupancy()
            outname_yield = os.path.join(self.folder, "{0}".format(folder2), "{0}".format(info_decay.get_decay()), "{0}".format(folder3), "{0}_{1}_{2}_InvMass_{3}_all_cent_step{4}_R_{5}_{6}_occup_{7}_{8}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), type_of_histo, info_decay.get_step_cent(), V0_radius_min, V0_radius_max, occup1, occup2))
        if is_mc_and_data:
            V0_radius_min = info_decay.get_V0_radius_min()
            V0_radius_max = info_decay.get_V0_radius_max()
            cent1, cent2 = info_decay.get_cent()
            occup1, occup2 = info_decay.get_occupancy()
            outname_yield = os.path.join(self.folder, "{0}".format(folder2), "{0}".format(info_decay.get_decay()), "{0}".format(folder3), "{0}_{1}_{2}_InvMass_{3}_cent{4}_{5}_mc_and_data_R_{6}_{7}_occup_{8}_{9}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), type_of_histo, cent1, cent2, V0_radius_min, V0_radius_max, occup1, occup2))
        if plotting == "RV0s" or plotting == "diff_cent_and_diff_RV0s":
            cent1, cent2 = info_decay.get_cent()
            occup1, occup2 = info_decay.get_occupancy()
            outname_yield = os.path.join(self.folder, "{0}".format(folder2), "{0}".format(info_decay.get_decay()), "{0}".format(folder3), "{0}_{1}_{2}_InvMass_{3}_all_RV0_cent{4}_{5}_occup_{6}_{7}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), type_of_histo, cent1, cent2, occup1, occup2))
        if plotting == "occupancies" or plotting == "diff_cent_and_diff_occupancies":
            cent1, cent2 = info_decay.get_cent()
            outname_yield = os.path.join(self.folder, "{0}".format(folder2), "{0}".format(info_decay.get_decay()), "{0}".format(folder3), "{0}_{1}_{2}_InvMass_{3}_all_occup_cent{4}_{5}_R_{6}_{7}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), type_of_histo, cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max()))
        if is_mc_and_data:
            plot_yield.plot_histo(type_of_histo, yield_list, TString(outname_yield), info_decay, plotting, ratio, scale_pT, is_mc_and_data, period_mc)
        else:
            plot_yield.plot_histo(type_of_histo, yield_list, TString(outname_yield), info_decay, plotting, ratio, scale_pT, is_mc_and_data)
        del plot_yield

    def write_yield_all_cent_Rcp_pdf(self, info_decay, raw_or_corr_yield, yield_list, is_other_file = False, energy_of_other_file = "", system_of_other_file = "", err_stat_list = [], err_syst_list = []):
        plot_Rcp = Plot_yield(info_decay.get_decay(), info_decay.get_meson())
        if raw_or_corr_yield == "RawYieldRcp":
            folder2 = "inv_mass_analysis"
            folder3 = "raw_yield"
        else:
            folder2 = "inv_yield_calculation"
            folder3 = "corr_yield"
        if is_other_file:
            outname_yield = os.path.join(self.folder, "Comparison", "other_data_Rcp.pdf")
            energy, system, period = energy_of_other_file, system_of_other_file, ""
        else:
            outname_yield = os.path.join(self.folder, "{0}".format(folder2), "{0}".format(info_decay.get_decay()),"{0}".format(folder3),"{0}_{1}_{2}_InvMass_{3}_{4}_all_cent_step{5}.pdf".format(info_decay.get_date(), info_decay.get_decay(), info_decay.get_period(), raw_or_corr_yield, info_decay.set_fitname(0), info_decay.get_step_cent()))
            energy, system , period = info_decay.get_energy(), info_decay.get_system(), info_decay.get_period()
        plot_Rcp.plot_Rcp(info_decay.get_all_cent(), raw_or_corr_yield, yield_list, TString(outname_yield), period, info_decay.get_ssname_list(), energy, system, False, err_stat_list, err_syst_list, is_other_file)
        del plot_Rcp

    def write_eff_acc_one_cent_pdf(self, info_decay, histo, name, ssname, teff = False):
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        name_plot = os.path.join(self.folder, "inv_yield_calculation", "{0}".format(info_decay.get_decay()),"eff_acc", "{0}_cent_{1}_{2}_R_{3}_{4}_occu_{5}_{6}.pdf".format(name, cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        plot_eff_acc = Plot_efficiency()
        plot_eff_acc.plot_eff_acc(histo, TString(name_plot), name, info_decay, teff)
        del plot_eff_acc

    def write_eff_acc_all_syst_pdf(self, info_decay, histo_list, name, plotting, ratio = False):
        if plotting == "centralities":
            occup1, occup2 = info_decay.get_occupancy()
            name_plot = os.path.join(self.folder, "inv_yield_calculation", "{0}".format(info_decay.get_decay()),"eff_acc", "{0}_all_cent_step{1}_R_{2}_{3}_occup_{4}_{5}.pdf".format(name, info_decay.get_step_cent(), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        if plotting == "RV0s" or plotting == "diff_cent_and_diff_RV0s":
            cent1, cent2 = info_decay.get_cent()
            occup1, occup2 = info_decay.get_occupancy()
            name_plot = os.path.join(self.folder, "inv_yield_calculation", "{0}".format(info_decay.get_decay()),"eff_acc", "{0}_all_RV0s_cent_step{1}_cent_{2}_{3}_occup_{4}_{5}.pdf".format(name, info_decay.get_step_cent(), cent1, cent2, occup1, occup2))
        if plotting == "occupancies" or plotting == "diff_cent_and_diff_occupancies":
            cent1, cent2 = info_decay.get_cent()
            name_plot = os.path.join(self.folder, "inv_yield_calculation", "{0}".format(info_decay.get_decay()),"eff_acc", "{0}_all_occupancies_step{1}_R_{2}_{3}_cent_{4}_{5}.pdf".format(name, info_decay.get_step_cent(), info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), cent1, cent2))
        plot_eff_acc = Plot_efficiency()
        plot_eff_acc.plot_eff_acc_all_syst(histo_list, TString(name_plot), info_decay, name, plotting, ratio)
        del plot_eff_acc

    def write_eff_vs_Teff_pdf(self, info_decay, eff_acc, teff_acc, isys):
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        name_plot = os.path.join(self.folder, "inv_yield_calculation", "{0}".format(info_decay.get_decay()),"eff_acc", "Comp_eff_vs_Teff_cent{0}_{1}_R_{2}_{3}_occup_{4}_{5}.pdf".format(cent1, cent2, info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max(), occup1, occup2))
        plot_eff_acc = Plot_efficiency()
        plot_eff_acc.comp_eff_Teff(eff_acc, teff_acc, TString(name_plot), info_decay)
        del plot_eff_acc

