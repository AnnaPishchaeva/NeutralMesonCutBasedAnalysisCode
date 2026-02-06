# This code was written by Julia Schlägel (July 2024)

import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX
import numpy as np

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

class Corrected_yield:
    def __init__(self):
        print("#################################################")
        print("Corrected_yield.py")
        print("#################################################")
        print("default constructor")

    def calc_corr_yield(self, info_decay, raw_data, eff_acc_hist, scale_pT):
        print("Scaling by pT:", scale_pT)
        raw = raw_data
        if raw:
            print("raw yield loaded")
        eff_acc = eff_acc_hist
        if eff_acc:
            print("eff_acc_loaded")
        deltarap = 1.6

        print(type(eff_acc))
        h_corrected = raw.Clone("h_corrected")
        h_corrected.Divide(eff_acc)
        h_corrected.Scale(1/deltarap)
        h_corrected.Scale(1/(2* np.pi))

        if scale_pT == "true":
            info_decay.scale_by_pt(h_corrected)

        return h_corrected

    def calc_pythia(self, info_decay): #This is the PYTHIA 8 model prediction
        hist =  info_decay.get_hist("hPt")
        deltay = 1.8
        if hist:
            print("histogram imported")
        Nev = info_decay.get_Nev_cent()
        print("Number of Events: ", Nev)

        h1 = hist.Clone("pythia")

        for bin in range (1, h1.GetXaxis().GetNbins()):
            bin_content = h1.GetBinContent(bin)
            bin_width = h1.GetBinWidth(bin)
            bin_error = h1.GetBinError(bin)
            if bin_width != 0:
                new_bin_content = bin_content / (bin_width)
                h1.SetBinContent(bin, new_bin_content)
                new_bin_error = bin_error / (bin_width)
                h1.SetBinError(bin, new_bin_error)

        h1.Scale(1 / Nev)
        h1.Scale(1 / deltay)
        h1.Scale(1 / (2 * np.pi))

        #This was just made, because there seems to be a bugg: the value when the binning changes is set to 0 so at
        # pT = 5 the value of the histogram is 0
        # this part of the code just cuts out this one strange bin
        x_values = []
        y_values = []
        bin_change = h1.FindBin(5)
        print("++++++ bin change: ", bin_change, h1.GetBinCenter(bin_change))
        for bin in range(1, h1.GetNbinsX() + 1):
        #for bin in range(1, bin_change +1):
            bin_content = h1.GetBinContent(bin)
            x_values.append(h1.GetBinCenter(bin))
            y_values.append(bin_content)
        print("x_values: ", x_values, "y_values: ", y_values)
        del x_values[bin_change-1]
        del y_values[bin_change-1]
        print("x_values after removing: ", x_values, "y_values: ", y_values)
        graph = ROOT.TGraph(len(x_values), np.array(x_values), np.array(y_values))

        return graph
