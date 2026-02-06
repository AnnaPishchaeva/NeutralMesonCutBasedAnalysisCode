import ROOT
import yaml
import math
import os, sys, shutil
import datetime
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX
import numpy as np
from array import array
# from PhotonMomentumResolution import FitInvMassForPt
#from FitInvMassForPt import PairAnalyzer
#from PlotRawYield import PlotRawYieldInvMass

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

class Efficiency:
    def __init__(self):
        print("#################################################")
        print("Efficiency.py")
        print("#################################################")

    def shift_histo(self, hist, shift):
        nbins = hist.GetNbinsX()
        xbins = [hist.GetBinLowEdge(i)+ shift for i in range(1, nbins+2)]
        new_hist = ROOT.TH1F(hist.GetName() + "_shifted", hist.GetTitle(), nbins, array("d", xbins))
        for i in range(1, nbins + 1):
            new_hist.SetBinContent(i, hist.GetBinContent(i))
            new_hist.SetBinError(i, hist.GetBinError(i))
        return new_hist

    def eff_and_acc_vs_pt(self, scale_BR, raw_mc, info_decay):
        #Get wanted histograms:
        h1 = raw_mc
        print("First Histogram rebinned")

        h2 = info_decay.rebin_hist("hPt") #rebinning the hist that is in mc to the pt range that we specified in the config file
        self.cent1, self.cent2 = info_decay.get_cent()
        info_decay.scale_by_Nev(h2)
        pt = info_decay.get_bin_var()
        print("####### pt range used: ", pt)

        print("second histogram rebinned")
        if h1 is None:
            print(f"Failed to load histogram: hPt_Pi0")

        h_eff_acc = h1.Clone("h_E_A")

        h_eff_acc.Divide(h2)
        if scale_BR == "true":
            info_decay.scale_by_BR(h_eff_acc)

        for i in range(0, len(pt)):
            print("###########Efficiency, values and Errors" , pt[i] ,"########")
            print("Value of raw mc: ", h1.GetBinContent(i), "±", h1.GetBinError(i))
            print("Value of hPt all: ", h2.GetBinContent(i), "±", h2.GetBinError(i))
            print("value of eff x acc x BR: ", h_eff_acc.GetBinContent(i), "±", h_eff_acc.GetBinError(i))
            if h_eff_acc.GetBinError(i) != 0:
                rel = h_eff_acc.GetBinContent(i) / h_eff_acc.GetBinError(i)
                print("relaive error: ", rel)
        print("eff and acc hist type:", type(h_eff_acc))
        return h_eff_acc

    def teff_and_acc_vs_pt(self, raw_mc, info_decay):
        #h1_b = raw_mc
        h1 = info_decay.scale_by_BR(raw_mc)
        h2 = info_decay.rebin_hist("hPt")
        info_decay.scale_by_Nev(h2)

        shift = 0.1
        shifted_h1 = self.shift_histo(h1, shift)
        shifted_h2 = self.shift_histo(h2, shift)

        h_eff_acc = ROOT.TEfficiency(shifted_h1, shifted_h2) #passed, total histograms
        h_eff_acc.SetStatisticOption(ROOT.TEfficiency.kFCP) #returns the probability of an event being larger than passed events; total
        #h_eff_acc.SetStatisticOption(ROOT.TEfficiency.kBBayesian)

        print("eff and acc hist type:", type(h_eff_acc))
        return h_eff_acc