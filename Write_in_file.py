# This code was written by Julia Schlägel (July 2024)

import ROOT
import yaml
import pickle
from ctypes import *
import numpy as np
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

class Write_in_file: #in this class all central things are done that are needed in multiple macros like loading a file from a rootfile,...
    def __init__(self, outfile_name):
        print("#################################################")
        print("Write_in_file.py")
        self.histo = {}
        self.output_file = ROOT.TFile(outfile_name, "RECREATE")
        print("output file = ", self.output_file)
        print("#################################################")

    def __del__(self):
        if self.output_file.IsOpen():
            print("close input root file.");
            self.output_file.Close();

    def writing_in_file(self, folder, outlist_inv_mass_histos):
        dir = self.output_file.mkdir(folder)
        dir.cd()
        for histo in outlist_inv_mass_histos:
            histo.Write()
        self.output_file.Close()

    def writing_in_one_histo(self, folder, histo, new_name=None):
        if new_name is not None:
            histo.SetName(new_name)

        if folder not in self.histo:
            self.histo[folder] = []
        self.histo[folder].append(histo)

        try:
            for folder, hist_list in self.histo.items():
                dir = self.output_file.mkdir(folder)
                dir.cd()
                for hist in hist_list:
                    hist.Write()
            self.output_file.Close()
            print(f"Histograms successfully stored in {self.output_file}")
        except Exception as e:
            print(f"Error by storing of the histos:{e}")


    def writing_in_pkl(self, info_decay, histo, obj):
        cent1, cent2 = info_decay.get_cent()
        occup1, occup2 = info_decay.get_occupancy()
        filename = f'pkl/{obj}_{info_decay.get_system()}_{info_decay.get_energy()}_{info_decay.get_meson()}_{info_decay.get_typ()}_{info_decay.get_decay()}_{cent1}_{cent2}_R_{info_decay.get_V0_radius_min()}_{info_decay.get_V0_radius_max()}_occup_{occup1}_{occup2}.pkl'
        for i in range(1, histo.GetNbinsX() + 1):
            bin_content = histo.GetBinContent(i)
            bin_center = histo.GetBinCenter(i)
            print(f"Bin {i}: Center = {bin_center}, Content = {bin_content}")
        with open(filename, 'wb') as pkl_file:
            pickle.dump(histo, pkl_file)
