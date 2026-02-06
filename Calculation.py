# This code was written by Alica Enderich (Febuary, 2024)

import ROOT
from ROOT import TH1D, TH2D, TH3D, TGraph, TGraphErrors, TLatex, TPaveText, TLegend, TCanvas, TPad, TString, TLine
from ROOT import gROOT, gSystem, gStyle, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kFullCircle
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle, kTRUE

from Pair_analyzer import Pair_analyzer
from Efficiency import Efficiency
from Corrected_yield import Corrected_yield
from Write_in_file import Write_in_file

def extract_inv_mass_yield(info_decay, isys):
    print("Calculation.py: in the function extract_inv_mass_yield for the subsystem ", isys)
    inv_mass_raw_yield = Pair_analyzer(info_decay)
    outlist_inv_mass_histos = inv_mass_raw_yield.analyze_ptspectrum(info_decay, isys)
    print("Calculation.py: the function extract_inv_mass_yield is complete")

    return outlist_inv_mass_histos

def calc_eff_acc(raw_mc, info_decay, w_file, ssname):
    eff_class = Efficiency()

    eff_acc = eff_class.eff_and_acc_vs_pt("false", raw_mc, info_decay)
    eff_acc_scaled = eff_class.eff_and_acc_vs_pt("true", raw_mc, info_decay)
    teff_acc = eff_class.teff_and_acc_vs_pt(raw_mc, info_decay)

    return eff_acc, eff_acc_scaled, teff_acc

def calc_corr_yield(raw_data, eff_acc, info_decay_mc, info_decay_data):
    corr_class = Corrected_yield()
    pyth = corr_class.calc_pythia(info_decay_mc)
    corr = corr_class.calc_corr_yield(info_decay_data, raw_data, eff_acc, "false")
    corr_scaled = corr_class.calc_corr_yield(info_decay_data, raw_data, eff_acc, "true")

    return pyth, corr, corr_scaled