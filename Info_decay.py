# This code was written by Julia Schlägel (July 2024)
import sys
import ROOT
import yaml
from ctypes import *
import numpy as np
import warnings
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

class Info_decay: #in this class all central things are done that are needed in multiple macros like loading a file from a rootfile,...
    def __init__(self, is_diff_cent, is_compare_RV0, is_compare_occupancy, is_diff_RV0, all_cent, is_start_pt, date, suffix, energy, system):
        print("************************************************")
        print("Info_decay.py")
        self.is_start_pt = is_start_pt
        self.date = date
        self.all_cent = all_cent
        self.suffix = suffix
        self.energy = energy
        self.system = system
        self.ssname_list = []
        self.is_diff_cent = is_diff_cent
        self.is_diff_RV0 = is_diff_RV0
        self.is_compare_occupancy = is_compare_occupancy
        self.is_compare_RV0 = is_compare_RV0
        if self.is_diff_cent and not self.is_compare_occupancy and not self.is_compare_RV0:
            self.plotting_sys = "centralities"
        if self.is_diff_cent and self.is_compare_RV0:
            self.plotting_sys = "diff_cent_and_diff_RV0s"
        if self.is_diff_cent and self.is_compare_occupancy:
            self.plotting_sys = "diff_cent_and_diff_occupancies"
        if self.is_diff_RV0:
            self.plotting_sys = "RV0s"

        print(sys._getframe().f_code.co_name);
        print("all_cent =", self.all_cent)
        print("is_diff_cent = ", self.is_diff_cent)
        print("is_diff_RV0 = ", self.is_diff_RV0)
        print("is_compare_RV0 = ", self.is_compare_RV0)
        print("plotting_sys = ", self.plotting_sys)
        print("************************************************")

    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input root file.");
            self.rootfile.Close();

    def set_typ(self, typ):
        print("Info_decay.py: setting type of the dataset...")
        self.typ = typ
        print("type = ", self.typ)

    def set_meson(self, meson):
        print("Info_decay.py: setting meson...")
        self.meson = meson
        print("meson = ", self.meson)

    def set_decay(self, decay):
        print("Info_decay.py: setting decay and decay channel...")
        self.decay = decay
        self.set_decay_channel(self.decay)
        print("decay = ", self.decay)
        print("decay channel= ", self.decay_channel)

    def set_input_file(self, filename):
        print("Info_decay.py: setting input file...")
        self.filename = filename
        self.rootfile = TFile.Open(filename, "READ")
        print("input = ",filename)

    def set_config_file(self, config_file):
        print("Info_decay.py: setting config file of a particle and its elements...")
        with open(config_file, "r", encoding="utf-8") as config_yml:
            self.config = yaml.safe_load(config_yml)
        self.period = self.config["common"]["period"]
        self.plotting_type = self.config["common"]["plotting_type"]
        self.background_fit = self.config["common"]["background_fit"]
        self.arr_pt = np.array(self.config["common"]["pt_bin"],dtype=float);
        self.system_of_config = self.config["common"]["system"]
        self.energy_of_config = self.config["common"]["energy"]
        self.occupancy_list = np.array(self.config["common"]["occupancy_bin"],dtype=str);
        self.RV0_list = np.array(self.config["common"]["RV0_bin"],dtype=float);
        if self.energy_of_config != self.energy:
            print("energy of the main config file: ", self.energy)
            print("energy of the particle config file: ", self.energy_of_config)
            raise ValueError("Energy of the main config doesn't coincide with the energy of the particle config")
        if self.system_of_config != self.system:
            print("system of the main config file: ", self.system)
            print("system of the particle config file: ", self.system_of_config)
            raise ValueError("System of the main config doesn't coincide with the system of the particle config")
        if len(self.occupancy_list) <=2 and self.is_compare_occupancy:
            raise ValueError("Not enough occup values in the occupancy list of the particle config file to make occupancy comparison!")
        if len(self.RV0_list) <=2 and self.is_compare_RV0():
            raise ValueError("Not enough RV0 values in the occupancy list of the particle config file to make RV0 comparison!")
        if self.is_start_pt:
            self.start_pt = self.config["common"]["start_pT_bin"]
        else:
            self.start_pt = 0
        self.set_config_parameters()

        print("the config file =  ", config_file)
        print("the plotting type = ", self.plotting_type)
        print("the background fit = ", self.background_fit)
        print("period = ", self.period)
        print("energy = ", self.energy)
        print("system = ", self.system)
        print("pT array = ", self.arr_pt);
        print("start pT =", self.start_pt)
        print("occupancy array = ", self.occupancy_list)
        print("RV0 array = ", self.RV0_list)

    def set_subsystem(self, isys):
        print("Info_decay.py: setting subsystem of the particle config file...")
        if self.is_diff_RV0:
            self.ssname = self.config[self.typ]['subsystems']["diff_RV0"][isys]['name'];
        if self.is_diff_cent:
            self.ssname = self.config[self.typ]['subsystems']["diff_cent"][isys]['name'];
        self.ssname_list.append(self.ssname)
        print("analyze subsystem: ", self.ssname)

    def set_config_parameters(self):
        print("Info_decay.py: setting config parameters...")
        if self.is_diff_RV0:
            self.subsystems = self.config[self.typ]['subsystems']["diff_RV0"]
            self.nsys = len(self.config[self.typ]['subsystems']["diff_RV0"])
        if self.is_diff_cent:
            self.subsystems = self.config[self.typ]['subsystems']["diff_cent"]
            self.nsys = len(self.config[self.typ]['subsystems']["diff_cent"])
        print("Number of subsystems: ", self.nsys);

    def set_subsystem_params(self, isys):
        print("Info_decay.py: setting ranges and fit functions from the particle config file...")
        self.initial_fit_parameters = self.subsystems[isys].get('initial_fit_parameters', [])
        self.fit_parameters_name = self.subsystems[isys].get('fit_parameters_name', [])
        self.fit_param_limit_min = self.subsystems[isys].get('fit_param_limit_min', []) #fit_limit_min
        self.fit_param_limit_max = self.subsystems[isys].get('fit_param_limit_max', []) #fit_limit_max
        print("initial_fit_parameters = ", self.initial_fit_parameters)
        print("fit_parameters_name = ", self.fit_parameters_name)
        print("fit_param_limit_min = ", self.fit_param_limit_min)
        print("fit_param_limit_max = ", self.fit_param_limit_max)

        self.fit_func = self.subsystems[isys].get('fit_func', [])
        if len(self.fit_func) != (len(self.arr_pt)-1) and len(self.fit_func) != 1:
            print("number of fit functions: ", len(self.fit_func))
            print("number of pT ranges: ", len(self.arr_pt)-1)
            raise ValueError("Number of fit functions does not equal to the number of pT ranges")
        print("fit_func: ", self.fit_func)

        self.fit_range = self.subsystems[isys].get('fit_range', [])
        if len(self.fit_range) != (len(self.arr_pt)-1) and len(self.fit_range) != 1:
            print("number of fit ranges: ", len(self.fit_range))
            print("number of pT ranges: ", len(self.arr_pt)-1)
            raise ValueError("Number of fit ranges does not equal to the number of pT ranges")
        print("fit_range: ", self.fit_range)

        self.scale_mix_range = self.subsystems[isys].get('scale_mix_range', [])
        if len(self.scale_mix_range) != (len(self.arr_pt)-1) and len(self.scale_mix_range) != 1:
            print("number of ranges for scaling mix inv mass histogram: ", len(self.scale_mix_range))
            print("number of pT ranges: ", len(self.arr_pt)-1)
            raise ValueError("Number of ranges for scaling mix inv mass histogram does not equal to the number of pT ranges")
        print("scale_mix_range: ", self.scale_mix_range)

        self.yield_range = self.subsystems[isys].get('yield_range', [])
        if len(self.yield_range) != (len(self.arr_pt)-1) and len(self.yield_range) != 1:
            print("number of ranges for raw yield calculation: ", len(self.yield_range))
            print("number of pT ranges: ", len(self.arr_pt)-1)
            raise ValueError("Number of ranges for raw yield calculation does not equal to the number of pT ranges")
        print("yield_range: ", self.yield_range)

    def set_centrality(self, ssname):
        print("Info_decay.py: setting centrality of the subsystem of the particle config file...")
        if self.is_diff_cent:
            if "_" in ssname:
                splitlist = ssname.split("_")
                centralities = splitlist[-1]
                if len(splitlist) == 4:
                    centralities = splitlist[-2] + "-" + splitlist[-1]
                if "-" in centralities:
                        self.cent1 = int(centralities.split("-")[0])
                        self.cent2 = int(centralities.split("-")[1])
                else:
                    raise ValueError("Centrality format not recognized")
            else:
                self.cent1 = 0
                self.cent2 = 100
            if self.cent2 < 100 or self.system == "pp" or self.system == "pO" or self.system == "OO":
                self.step_cent = self.cent2 - self.cent1
                print("centrality step = ", self.step_cent)
            print("centrality 1 = ",self.cent1, "centrality 2 = ", self.cent2)
        if self.is_diff_RV0:
            self.cent1 = 0
            self.cent2 = 80
            self.step_cent = self.cent2 - self.cent1
            if self.system != "Pb-Pb" and self.system !="PbPb" and self.system !="Pb--Pb":
                self.cent2 = 100
            print("centrality step = ", self.step_cent)
            print("centrality 1 = ",self.cent1, "centrality 2 = ", self.cent2)


    def set_RV0(self, ssname):
        print("Info_decay.py: setting RV0 os the subsystem...")
        if self.is_diff_cent:
            self.V0_radius_min = self.RV0_list[0]
            self.V0_radius_max = self.RV0_list[len(self.RV0_list)-1]
        if self.is_diff_RV0:
            digit_list = [str(character) for character in ssname if character.isdigit()]
            if len(digit_list)>1:
                del digit_list[0] #deletes 0 from pi0
                del digit_list[0] #deletes 0 from V0
            if len(digit_list) == 1:
                del digit_list[0] #deletes 0 from pi0 from the overall subsystem
            if len(digit_list) == 0:
                self.V0_radius_min = self.RV0_list[0]
                self.V0_radius_max = self.RV0_list[len(self.RV0_list)-1]
            if len(digit_list) == 2:
                self.V0_radius_min = int(digit_list[0])
                self.V0_radius_max = int(digit_list[1])
            if len(digit_list) == 3:
                self.V0_radius_min = int(digit_list[0])
                self.V0_radius_max = int(digit_list[1]+digit_list[2])
            if len(digit_list) == 4:
                self.V0_radius_min = int(digit_list[0]+digit_list[1])
                self.V0_radius_max = int(digit_list[2]+digit_list[3])
        print("radius V0 min: ", self.V0_radius_min, "radius V0 max: ", self.V0_radius_max)


    def set_occupancy(self, ssname):
        print("Info_decay.py: setting occupancy of the subsystem...")
        if self.is_diff_cent or self.is_diff_RV0:
            self.occupancy_min = self.occupancy_list[0]
            self.occupancy_max = self.occupancy_list[len(self.occupancy_list)-1]

    def set_decay_channel(self, decay):
        if decay == "Dalitz":
            self.decay_channel = "e^{+}e^{-}#gamma"
        if decay == "PCM":
            self.decay_channel = "#gamma#gamma"

    def get_filename(self):
        print("Info_decay.py: getting name of the initial input root file...")
        return self.filename

    def get_config(self):
        return self.config

    def get_number_ss(self):
        print("Info_decay.py: getting number of subsystems of the particle config file...")
        return self.nsys

    def get_ssname(self):
        print("Info_decay.py: getting the subsystem of the particle config file...")
        return self.ssname

    def get_ssname_list(self):
        return self.ssname_list

    def get_energy(self):
        return self.energy

    def get_system(self):
        return self.system

    def get_period(self):
        return self.period

    def get_typ(self):
        print("Info_decay.py: getting the type of the dataset...")
        return self.typ

    def get_decay(self):
        return self.decay

    def get_meson(self):
        print("Info_decay.py: getting meson...")
        return self.meson

    def get_cent(self):
        print("Info_decay.py: getting centrality...")
        return self.cent1, self.cent2

    def get_step_cent(self):
        return self.step_cent

    def get_initial_fit_parameters(self):
        print("Info_decay.py: getting initial fit parameters...")
        return self.initial_fit_parameters

    def get_fit_parameters_name(self):
        print("Info_decay.py: getting fit parameters names...")
        return self.fit_parameters_name

    def get_fit_param_limit_min(self): #get_initial_lower_limit
        print("Info_decay.py: getting fit parameters limit min...")
        return self.fit_param_limit_min

    def get_fit_param_limit_max(self):#get_initial_upper_limit
        print("Info_decay.py: getting fit parameters limit max...")
        return self.fit_param_limit_max

    def get_fit_func(self):
        print("Info_decay.py: getting fit functions...")
        if len(self.fit_func) == 1:
            warnings.warn("all of the pT ranges have the same fit_func")
        return self.fit_func

    def get_fit_range(self):
        print("Info_decay.py: getting fit ranges...")
        if len(self.fit_range) == 1:
            warnings.warn("all of the pT ranges have the same fit range")
        return self.fit_range

    def get_scale_mix_range(self):
        print("Info_decay.py: getting scaling ranges for the mix histo...")
        if len(self.scale_mix_range) == 1:
            warnings.warn("all of the pT ranges have the same range for scaling mix inv mass")
        return self.scale_mix_range

    def get_yield_range(self):
        print("Info_decay.py: getting yield ranges...")
        if len(self.yield_range) == 1:
            warnings.warn("all of the pT ranges have the same range for calculating raw yield")
        return self.yield_range

    def get_arr_pt (self):
        return self.arr_pt

    def get_RV0_list(self):
        return self.RV0_list

    def get_occupancy_list(self):
        return self.occupancy_list

    def get_start_pt(self):
        return self.start_pt

    def get_date(self):
        return self.date

    def get_suffix(self):
        return self.suffix

    def get_all_cent(self):
        print("Info_decay.py: getting all cent, a.k.a. whether the whole range of centrality is plotted...")
        return self.all_cent

    def get_V0_radius_max(self):
        return self.V0_radius_max

    def get_V0_radius_min(self):
        return self.V0_radius_min

    def get_occupancy(self):
        print("Info_decay.py: getting occupancy...")
        return self.occupancy_min, self.occupancy_max

    def get_decay_channel(self):
        return self.decay_channel

    def get_is_diff_RV0(self):
        return self.is_diff_RV0

    def get_is_diff_cent(self):
        return self.is_diff_cent

    def get_is_compare_RV0(self):
        return self.is_compare_RV0

    def get_is_compare_occupancy(self):
        return self.is_compare_occupancy

    def get_plotting_sys(self):
        return self.plotting_sys

    def get_hist (self, hist_name):
        dir_name1 = "{0}/Generated".format(self.ssname)
        if self.meson == "pi0":
            dir_name2 = "Pi0"
        if self.meson == "eta":
            dir_name2 = "Eta"
        dir = self.rootfile.Get(dir_name1).Get(dir_name2)
        self.hist_name = dir.Get(hist_name)
        if not self.hist_name:
            raise ValueError("Histogram not found")
        return self.hist_name

    def get_Nev_cent (self): #function to get the number of events from the rootfile for particular centrality
        print("Info_decay.py: getting number of events...")
        direv = self.rootfile.Get("{0}/Event".format(self.ssname))
        events = direv.Get("after").Get("hCollisionCounter")
        self.Nev = events.GetBinContent(9)
        print("number of events = ", self.Nev)
        return self.Nev

    def get_bin_var(self): #function to get the pT binning that was set in the config file
        self.bin_var = self.config["common"]["pt_bin"]
        print(self.bin_var)
        return self.bin_var

    def scale_by_Nev(self, input_hist): #function to divide a histogram by the number of events
        Nev = self.get_Nev_cent()
        input_hist.Scale(1/Nev)
        return input_hist

    def scale_to_inv_yield(self, input_hist):
        #input_hist.Scale((57.8*1e-3))
        input_hist.Scale(59.4*1e-3)
        input_hist.Scale(1e12)

    def scale_by_pt(self, input_hist): # function to divide a histogram by pT
        p_T = self.get_bin_var()
        for bin in range (1, len(p_T)):
            bin_content = input_hist.GetBinContent(bin)
            bin_center = input_hist.GetBinCenter(bin)
            bin_error = input_hist.GetBinError(bin)

            new_bin_content = bin_content / (bin_center)
            input_hist.SetBinContent(bin, new_bin_content)
            new_bin_error = bin_error / (bin_center)
            input_hist.SetBinError(bin, new_bin_error)
        return input_hist

    def scale_by_BR(self, input_hist): # function to scale a histogram with the branching ratio
        if self.meson == "eta":
            if self.decay == "PCM":
                Br = 0.3936
            if self.decay == "Dalitz":
                Br = 0.0069
        if self.meson == "pi0":
            if self.decay == "PCM":
                Br = 0.98823
            if self.decay == "Dalitz":
                Br = 0.01174
        input_hist.Scale(1/Br)
        return input_hist

    def rebin_hist(self, hist_name):
        bin_var = self.get_bin_var()
        hist = self.get_hist(hist_name)
        new_hist = TH1D("p_T1", "p_T1", len(bin_var) - 1, np.array(bin_var))
        new_hist.Sumw2()

        for bin in range(0, len(bin_var) - 1):
            b1 = bin_var[bin]
            b2 = bin_var[bin + 1] #basically value of x-axis
            bin_b1 = hist.GetXaxis().FindBin(b1 + 1e-6) #the bin nmber of the value of x-axis
            bin_b2 = hist.GetXaxis().FindBin(b2 - 1e-6)

            error = c_double(0.0)
            # content = hist.IntegralAndError(bin_b1, bin_b2, error, "")
            content = 0
            print("bin1: ", bin_b1, bin_b2)
            for i in range(bin_b1, bin_b2 +1):
                content = content+ hist.GetBinContent(i)
            error = np.sqrt(content)

            bin_center = (b1 + b2) / 2  # Take the center of the bin to fill
            bin_width = b2-b1 #basically the width of one bin
            print("in rebin function error:", error)
            print("content: ", content)
            #
            new_content = content / bin_width
            # new_hist.Fill(bin_center, new_content)  # Filling with this center of each bin
            new_hist.SetBinContent(new_hist.FindBin(bin_center), new_content)
            #new_hist.SetBinError(bin, error) # / bin_width)  # calculate error to the right bin#
            new_hist.SetBinError(new_hist.FindBin(bin_center), error)


        if new_hist is None:
            print(f"Failed to load histogram: {hist_name}")


        # Create canvas and draw histogram
        canvas = TCanvas("acc", "acc", 0, 0, 900, 900)
        pad = canvas.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        new_hist.Draw("ESame")
        ROOT.SetOwnership(canvas, False)


        return new_hist
