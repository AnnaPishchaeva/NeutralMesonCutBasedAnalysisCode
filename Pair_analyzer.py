# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3

import os, sys, shutil
import numpy as np
import math
import ctypes
import ROOT
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TMath, TH1D, TROOT, gROOT
from ctypes import *

class Pair_analyzer:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, info_decay):
        print("#################################################")
        print("Pair_analyzer.py", flush=True)
        self.xtitle = "#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})"
        self.ytitle = "#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})"
        self.meson = info_decay.get_meson()
        self.typ = info_decay.get_typ()
        self.filename = info_decay.get_filename()
        self.decay = info_decay.get_decay()
        self.ssname = info_decay.get_ssname()
        self.rootfile = TFile.Open(self.filename, "READ");
        self.rootdir = self.rootfile.Get(self.ssname);
        self.list_ev = self.rootdir.Get("Event");
        self.list_pair = self.rootdir.Get("Pair");
        self.arr_pt = np.array([0,1,2,3,4,5], dtype=float);
        self.cent1, self.cent2 = info_decay.get_cent()
        self.arr_pt = info_decay.get_arr_pt()
        self.nev = info_decay.get_Nev_cent()
        self.all_cent = False
        self.isys = 0
        self.f1total = TF1("GaussExpLinear",
               "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x)+\
               (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x)", 0.04 ,0.24);

        self.f1total.SetNpx(1000);
        print("using info_decay object, the elements of the info_decay were acquired")
        print("#################################################")

    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input root file.");
            self.rootfile.Close();
    def get_2D_histo_same(self, ssname, same_mix_dir):
        print("Pair_analyzer.py: in the function get_2D_histo_same", flush=True)
        #rootfile = TFile.Open(filename_data, "READ")
        # ss_dir = self.rootfile.Get("pi0eta-to-gammagamma-pcmpcm")
        ss_dir = self.rootfile.Get(ssname)
        pair_dir = ss_dir.Get("Pair")
        #same_dir = pair_dir.Get("same")
        same_dir = pair_dir.Get(same_mix_dir)
        hs_sparse = same_dir.Get("hs")
        histo_proj = hs_sparse.Projection(1,0)
        print("type of hMggPt: ", type(histo_proj))
        return histo_proj

    def get_2D_histo_primary(self, ssname_mc):
        print("Pair_analyzer.py: in the function get_2D_histo_primary", flush=True)
        #ss_dir = self.rootfile.Get("pi0eta-to-gammagamma-mc-pcmpcm")
        ss_dir = self.rootfile.Get(ssname_mc)
        pair_dir = ss_dir.Get("Pair")
        if self.meson == "pi0":
            meson_dir = pair_dir.Get("Pi0")
        if self.meson == "eta":
            meson_dir = pair_dir.Get("Eta")
        hs_sparse = meson_dir.Get("hs_Primary")
        histo_proj = hs_sparse.Projection(1,0)
        print("type of hs_Primary: ", type(histo_proj))
        return histo_proj

    def slice_histogram(self, h2,x0,x1,axis,isdiff):
        print("Pair_analyzer.py: in the function slice_histogram", flush=True)
        h1 = 0;
        delta = 1e-6;
        if "x" in axis.lower():
            bin0 = h2.GetYaxis().FindBin(x0 + delta);
            bin1 = h2.GetYaxis().FindBin(x1 - delta);
            h1 = h2.ProjectionX("h1prjx_{0}".format(h2.GetName()),bin0,bin1,"e");
        elif "y" in axis.lower():
            bin0 = h2.GetXaxis().FindBin(x0 + delta);
            bin1 = h2.GetXaxis().FindBin(x1 - delta);
            h1 = h2.ProjectionY("h1prjy_{0}".format(h2.GetName()),bin0,bin1,"e"); # "e" for error calculation

        if isdiff and h1.Class() == TH1D.Class():
            h1.Scale(1.,"width");
        return h1;

    def set_fit_function(self, funct):
        print("Pair_analyzer.py: in the function set_fit_function", flush=True)
        if self.typ == "data":
            if funct == "gausexplinear":
                self.f1total = TF1("GaussExpLinear",
                "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x)+\
                (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x)", 0.,1.);
                self.func = "GaussExpLinear"
            elif funct == "gausexpquadratic":
                self.f1total = TF1("GaussExpQuadratic",
                "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x +[6]*(x-[1])^2)+\
                (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x + [6]*x*x)", 0,1);
                self.func = "GaussExpQuadratic"
            elif funct == "gausexpcube":
                self.f1total = TF1("GaussExpCube",
                "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x +[6]*(x-[1])^2)+\
                (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x + [6]*x*x + [7]*x*x*x)", 0,1);
                self.func = "GaussExpCube"

        elif self.typ == "mc":
            self.f1total = TF1("fGaussExp",
                 "(x<[1]) * ([0]*(TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)) + TMath::Exp((x-[1])/[3])*(1. - TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2))))) + \
                 (x>=[1]) * ([0]*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)))", 0,1);
            self.func = "GaussExp"

        self.f1total.SetNpx(1000);

    def calc_FWHM(self, histo, params, paramsErr):
        print("Pair_analyzer.py: in the function calc_FWHM", flush=True)
        tf1_fwhm     = TF1("tf1_fwhm",   "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))))+\
            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))", 0,1);
        for i in range(4):
            tf1_fwhm.SetParameter(i, params[i])
            tf1_fwhm.SetParError(i, paramsErr[i])
        tf1_fwhm.SetNpx(1000);

        maximum = tf1_fwhm.GetMaximum()
        maximum_x = tf1_fwhm.GetMaximumX()
        half_maximum = maximum / 2
        left_x = tf1_fwhm.GetX(half_maximum, 0, maximum_x);
        right_x = tf1_fwhm.GetX(half_maximum, maximum_x, 1);
        FWHM = (right_x-left_x)/TMath.Sqrt(8*TMath.Log(2))  # in equivalents to sigma

        # Calculate PLUS-error with (FWHM+FWHM_err)
        tf1_fwhm_plus_err     = TF1("tf1_fwhm_plus_err",   "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))))+\
            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))", 0,1);
        for i in range(4):
            param = params[i] + paramsErr[i]
            tf1_fwhm_plus_err.SetParameter(i, param)

        maximum_plus = tf1_fwhm_plus_err.GetMaximum()
        maximum_x_plus = tf1_fwhm_plus_err.GetMaximumX()
        half_maximum_plus = maximum_plus / 2
        left_x_plus = tf1_fwhm_plus_err.GetX(half_maximum_plus, 0, maximum_x_plus);
        right_x_plus = tf1_fwhm_plus_err.GetX(half_maximum_plus, maximum_x_plus, 1);
        FWHM_plus = (right_x_plus-left_x_plus)/TMath.Sqrt(8*TMath.Log(2))

        # Calculate MINUS-error with (FWHM-FWHM_err)
        tf1_fwhm_minus_err     = TF1("tf1_fwhm_plus_err",   "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))))+\
            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))", 0,1);
        for i in range(4):
            param = params[i] - paramsErr[i]
            tf1_fwhm_minus_err.SetParameter(i, param)

        maximum_minus = tf1_fwhm_minus_err.GetMaximum()
        maximum_x_minus = tf1_fwhm_minus_err.GetMaximumX()
        half_maximum_minus = maximum_minus / 2
        left_x_minus = tf1_fwhm_minus_err.GetX(half_maximum_minus, 0, maximum_x_minus);
        right_x_minus = tf1_fwhm_minus_err.GetX(half_maximum_minus, maximum_x_minus, 1);
        FWHM_minus = (right_x_minus-left_x_minus)/TMath.Sqrt(8*TMath.Log(2))

        # Calculate TOTAL error
        FWHM_err = max(abs(FWHM_plus-FWHM), abs(FWHM_minus-FWHM))
        FWHM_err = FWHM_err/2
        print("FWHM calc: ", FWHM*1e3, "+/-", FWHM_err*1e3, ";", FWHM_plus*1e3, FWHM_minus*1e3)
        return tf1_fwhm, FWHM, FWHM_err

    def calc_FWHM_MC(self, histo, params, paramsErr):
        print("Pair_analyzer.py: in the function calc_FWHM_MC", flush=True)
        tf1_fwhm     = TF1("tf1_fwhm",   "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))))+\
            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))", 0,1);
        for i in range(4):
            tf1_fwhm.SetParameter(i, params[i])
            tf1_fwhm.SetParError(i, paramsErr[i])
        tf1_fwhm.SetNpx(1000);

        maximum = tf1_fwhm.GetMaximum()
        maximum_x = tf1_fwhm.GetMaximumX()
        half_maximum = maximum / 2
        left_x = tf1_fwhm.GetX(half_maximum, 0, maximum_x);
        right_x = tf1_fwhm.GetX(half_maximum, maximum_x, 1);
        FWHM = (right_x-left_x)/TMath.Sqrt(8*TMath.Log(2))  # in equivalents to sigma

        # Calculate PLUS-error with (FWHM+FWHM_err)
        tf1_fwhm_plus_err     = TF1("tf1_fwhm_plus_err",   "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))))+\
            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))", 0,1);
        for i in range(4):
            param = params[i] + paramsErr[i]
            tf1_fwhm_plus_err.SetParameter(i, param)

        maximum_plus = tf1_fwhm_plus_err.GetMaximum()
        maximum_x_plus = tf1_fwhm_plus_err.GetMaximumX()
        half_maximum_plus = maximum_plus / 2
        left_x_plus = tf1_fwhm_plus_err.GetX(half_maximum_plus, 0, maximum_x_plus);
        right_x_plus = tf1_fwhm_plus_err.GetX(half_maximum_plus, maximum_x_plus, 1);
        FWHM_plus = (right_x_plus-left_x_plus)/TMath.Sqrt(8*TMath.Log(2))

        # Calculate MINUS-error with (FWHM-FWHM_err)
        tf1_fwhm_minus_err     = TF1("tf1_fwhm_plus_err",   "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))))+\
            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))", 0,1);
        for i in range(4):
            param = params[i] - paramsErr[i]
            tf1_fwhm_minus_err.SetParameter(i, param)

        maximum_minus = tf1_fwhm_minus_err.GetMaximum()
        maximum_x_minus = tf1_fwhm_minus_err.GetMaximumX()
        half_maximum_minus = maximum_minus / 2
        left_x_minus = tf1_fwhm_minus_err.GetX(half_maximum_minus, 0, maximum_x_minus);
        right_x_minus = tf1_fwhm_minus_err.GetX(half_maximum_minus, maximum_x_minus, 1);
        FWHM_minus = (right_x_minus-left_x_minus)/TMath.Sqrt(8*TMath.Log(2))

        # Calculate TOTAL error
        FWHM_err = max(abs(FWHM_plus-FWHM), abs(FWHM_minus-FWHM))
        FWHM_err = FWHM_err/2  # in equivalents to sigma and divided by 2 for error propagation
        return tf1_fwhm, FWHM, FWHM_err

    def calculate_raw_yield(self, yield_range_for_one_pT, histo, params, params_err, dpt):
        print("Pair_analyzer.py: in the function calculate_raw_yield", flush=True)
        integral_min = histo.FindBin(params[1] - yield_range_for_one_pT[0])
        integral_max = histo.FindBin(params[1] + yield_range_for_one_pT[1])
        print("bin_width of integral: ", histo.GetBinWidth(2))
        error_integral = c_double(0.0)
        integral_histo = histo.IntegralAndError(integral_min, integral_max, error_integral) #without "width"
        print("integral_histo: ", integral_histo , error_integral)
        #rel = integral_histo / error_integral
        #print("relative error: ", rel)

        # subtract linear background
        if self.func == "GaussExpLinear":
            tf1 = TF1("tf1","[4]+[5]*x", 0,1)
            print("linear")
        if self.func == "GaussExpQuadratic":
            tf1 = TF1("tf1","[4]+[5]*x+[6]*x*x", 0,1)
            print("quadratic")
        if self.func == "GaussExpCube":
            tf1 = TF1("tf1","[4]+[5]*x+[6]*x*x+[7]*x*x*x", 0,1)
            print("cubic")
        tf1.SetParameter(5, params[5])
        tf1.SetParameter(4, params[4])
        if self.func == "GaussExpQuadratic" or self.func == "GaussExpCube":
            tf1.SetParameter(6, params[6])
            if self.func == "GaussExpCube":
                tf1.SetParameter(7, params[7])
        lin_integral_min = params[1] - yield_range_for_one_pT[0]
        lin_integral_max = params[1] + yield_range_for_one_pT[1]
        integral_background = tf1.Integral(lin_integral_min, lin_integral_max)
        integral_background = integral_background / histo.GetBinWidth(2)
        print("integral_background: ", integral_background)
        if self.func == "GaussExpLinear":
            error_fit = np.sqrt( ((lin_integral_max-lin_integral_min)*params_err[4])**2 +
                               (0.5*(lin_integral_max*lin_integral_max-lin_integral_min*lin_integral_min)*params_err[5])**2)
        if self.func == "GaussExpQuadratic":
            error_fit = np.sqrt( ((lin_integral_max-lin_integral_min)*params_err[4])**2 +
                               (0.5*(lin_integral_max*lin_integral_max-lin_integral_min*lin_integral_min)*params_err[5])**2 +
                               ((1/3)*(lin_integral_max*lin_integral_max*lin_integral_max-lin_integral_min*lin_integral_min*lin_integral_min)*params_err[6])**2)
        if self.func == "GaussExpCube":
            error_fit = np.sqrt( ((lin_integral_max-lin_integral_min)*params_err[4])**2 +
                               (0.5*(lin_integral_max*lin_integral_max-lin_integral_min*lin_integral_min)*params_err[5])**2 +
                               ((1/3)*(lin_integral_max*lin_integral_max*lin_integral_max-lin_integral_min*lin_integral_min*lin_integral_min)*params_err[6])**2+
                               ((1/4)*(lin_integral_max*lin_integral_max*lin_integral_max*lin_integral_max-lin_integral_min*lin_integral_min*lin_integral_min*lin_integral_min)*params_err[7])**2)
        error_fit = error_fit / histo.GetBinWidth(2)
        error = TMath.Sqrt(error_integral.value**2 + error_fit**2)
        signal = integral_histo - integral_background
        print("signal = ", signal)
        if signal !=0:
            stat_uncert = error/signal
        else:
            stat_uncert = -100
        print("stat. uncert = ", "error", error, "/", "signal", signal, "=", stat_uncert)
        error = error /dpt /self.nev
        raw_yield = (integral_histo - integral_background)/dpt/self.nev
        print("raw yield", raw_yield, " error of raw yield", error)
        print
        return raw_yield, error, signal, stat_uncert

    def calculate_raw_yield_MC(self, yield_range_for_one_pT, histo, params, params_err, dpt):
        print("Pair_analyzer.py: in the function calculate_raw_yield_MC", flush=True)
        integral_min = histo.FindBin(params[1] - yield_range_for_one_pT[0]) #
        integral_max = histo.FindBin(params[1] + yield_range_for_one_pT[1])
        error_integral = c_double(0.0)
        integral_histo = histo.IntegralAndError(integral_min, integral_max, error_integral) #without "width"
        error = error_integral.value / self.nev
        raw_yield = (integral_histo)/dpt/self.nev
        print("raw yield", raw_yield," error of raw yield", error, " signal ",integral_histo)
        if integral_histo != 0:
            stat_uncert = error_integral.value/integral_histo
        else:
            stat_uncert = -100
        print("stat. uncert = ", "error", error_integral.value, "/", "signal", integral_histo, "=", stat_uncert)
        #stat_uncert = 1/np.sqrt(integral_histo)
        return raw_yield, error, stat_uncert

    def calculate_significance(self, h1same, signal, yield_range_for_one_pT, params):
        print("Pair_analyzer.py: in the function calculate_significance", flush=True)
        signal_background = 0
        min = h1same.FindBin(params[1] - yield_range_for_one_pT[0])
        max = h1same.FindBin(params[1] + yield_range_for_one_pT[1])
        for ibin in range(min, max+1):
            signal_background += h1same.GetBinContent(ibin)
        significance = signal/np.sqrt(signal_background)
        print("significance = ", significance)
        return significance

    def set_2D_histograms(self, outlist):
        print("Pair_analyzer.py: in the function set_2D_histograms")

        ########################################
        #   Collision Counter, h2same, h2mix   #
        ########################################

        self.list_ev_af = self.list_ev.Get("after");
        self.list_ev_af.ls();

        h1ev   = self.list_ev_af.FindObject("hCollisionCounter");
        #nentries = h1ev.GetEntries();
        if h1ev:
            print("Object 'h1ev' found in list_ev_ss")
            #print("entries = {0} ".format(nentries));
        else:
            print("Object 'h1ev' found in list_ev_ss")


        h2same = self.get_2D_histo_same(self.ssname,"same")

        h2mix = self.get_2D_histo_same(self.ssname, "mix")

        h2same.Sumw2();
        h2mix.Sumw2();
        h2same.SetDirectory(0);
        h2mix .SetDirectory(0);
        if self.meson == "pi0":
            h2same.RebinX(1);
            h2mix .RebinX(1);
        elif self.meson == "eta":
            h2same.RebinX(2);
            h2mix .RebinX(2);
        h2same.Sumw2();
        h2mix.Sumw2();

        #outlist.append(h1ev);
        outlist.append(h2same);
        outlist.append(h2mix);

        return h2same, h2mix, h1ev

    def set_2D_histograms_MC(self, outlist):
        print("Pair_analyzer.py: in the function set_2D_histograms_MC")
        self.list_ev_ss_af = self.list_ev.Get("after");
        self.list_ev_ss_af.ls();
        h1ev   = self.list_ev_ss_af.FindObject("hCollisionCounter");

        h2mc = self.get_2D_histo_primary(self.ssname)

        h2mc.Sumw2();
        h2mc.SetDirectory(0);
        if self.meson == "pi0":
            h2mc.RebinX(1);
        elif self.meson == "eta":
            h2mc.RebinX(2);
        h2mc.Sumw2();
        #outlist.append(h1ev);
        outlist.append(h2mc);
        return h2mc, h1ev,

    def slice_histo_and_subtract_background(self, scale_mix_range_for_one_pT, isys, h2same, h2mix, pt1, pt2):
        print("Pair_analyzer.py: in the function slice_histo_and_subtract_background")
        h1same = self.slice_histogram(h2same, pt1, pt2, "x", False); # implemented "e" (error calc) in ProjectionX/ProjectionY in the used function
        h1same.SetName("h1mgg_same_pt_{0}_{1}_syst_{2}".format(pt1, pt2, isys));
        h1same.SetTitle("m_{{#gamma#gamma}}^{{same}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}} centrality {2:d}#%-{3:d}#%".format(pt1, pt2, self.cent1, self.cent2));
        h1mix  = self.slice_histogram(h2mix , pt1, pt2, "x", False);
        h1mix.SetTitle("m_{{#gamma#gamma}}^{{mix}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}} centrality {2:d}#%-{3:d}#%".format(pt1, pt2, self.cent1, self.cent2));
        h1mix.SetName("h1mgg_mix_pt_{0}_{1}_syst_{2}".format(pt1, pt2, isys));
        h1same.RebinX(2);
        # h1same.Scale(1./nev);
        h1same.Sumw2();
        h1mix.RebinX(2);
        # h1mix.Scale(1./nev);
        h1mix.Sumw2();
        h1same.SetDirectory(0);
        h1mix .SetDirectory(0);

        bin_min_same = h1same.FindBin(scale_mix_range_for_one_pT[0]);
        bin_max_same = h1same.FindBin(scale_mix_range_for_one_pT[1]);
        integral_same_event = h1same.Integral(bin_min_same, bin_max_same);
        integral_complete_same = h1same.Integral()
        bin_min_mix = h1mix.FindBin(scale_mix_range_for_one_pT[0]);
        bin_max_mix = h1mix.FindBin(scale_mix_range_for_one_pT[1]);
        integral_mixed_event = h1mix.Integral(bin_min_mix, bin_max_mix);
        integral_complete_mix = h1mix.Integral()
        h1mix_scaled = h1mix.Clone("h1mgg_mix_scaled_pt_{0}_{1}_syst_{2}_scale_{3}_{4}".format(pt1, pt2, isys, scale_mix_range_for_one_pT[0], scale_mix_range_for_one_pT[1]));
        if integral_mixed_event != 0:
            h1mix_scaled.Scale(integral_same_event / integral_mixed_event)
        else:
            print("Error: Division by zero. integral_mixed_event is zero.")
        h1mix_scaled.SetDirectory(0);
        h1subtracted = h1same.Clone("h1mgg_pt_{0}_{1}_syst_{2}_scale_{3}_{4}".format(pt1, pt2,isys, scale_mix_range_for_one_pT[0], scale_mix_range_for_one_pT[1]));
        h1subtracted.Add(h1mix_scaled, -1);
        h1subtracted.SetTitle("m_{{#gamma#gamma}}^{{sub.}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}} centrality {2:d}%-{3:d}%".format(pt1, pt2, self.cent1, self.cent2));
        if self.meson == "pi0":
            h1subtracted.GetXaxis().SetRangeUser(0, 0.3);
        if self.meson == "eta":
            h1subtracted.GetXaxis().SetRangeUser(0.2, 0.8)
        h1subtracted .SetDirectory(0);

        return h1subtracted, h1mix_scaled, h1same, h1mix

    def slice_histo(self, scale_mix_range_for_one_pT, isys, h2mc, pt1, pt2):
        print("Pair_analyzer.py: in the function slice_histo")
        h1mc = self.slice_histogram(h2mc, pt1, pt2, "x", False); # implemented "e" (error calc) in ProjectionX/ProjectionY in the used function
        h1mc.SetName("h1mgg_pt_{0}_{1}_syst_{2}_scale_{3}_{4}".format(pt1, pt2, isys, scale_mix_range_for_one_pT[0], scale_mix_range_for_one_pT[1]));
        h1mc.SetTitle("m_{{#gamma#gamma}}^{{mc}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}} centrality {2:d}%-{3:d}%".format(pt1, pt2, self.cent1, self.cent2));
        h1mc.RebinX(2);
        # h1mc.Scale(1./nev);
        h1mc.Sumw2();
        h1mc.SetDirectory(0);

        bin_min = h1mc.FindBin(scale_mix_range_for_one_pT[0]);
        bin_max = h1mc.FindBin(scale_mix_range_for_one_pT[1]);
        integral_mc = h1mc.Integral(bin_min, bin_max);
        integral_complete = h1mc.Integral()
        return h1mc

    #______________________________________________________________________
    def analyze_ptspectrum(self, info_decay, isys):
        print("Pair_analyzer.py: in the function analyze_ptspectrum for the subsystem ", isys)
        outlist = []
        npt = len(self.arr_pt);

        if self.typ == "data":
            h2same, h2mix, h1ev = self.set_2D_histograms(outlist)
        elif self.typ == "mc":
            h2mc, h1ev = self.set_2D_histograms_MC(outlist)
        ########################################
        #           h1parameter plots          #
        ########################################

        print("Pair_analyzer.py: in the function analyze_ptspectrum setting up histograms for inv. mass, parameters of the fit")
        h1chi2          = TH1F("h1chi2_param{0}".format(isys),   "chi2",                npt-1, self.arr_pt);
        h1amplitude     = TH1F("h1amplitude_param{0}".format(isys),   "amplitude",                npt-1, self.arr_pt);
        h1mean          = TH1F("h1mean_param{0}".format(isys),       "mean",                     npt-1, self.arr_pt);
        h1sigma         = TH1F("h1sigma_param{0}".format(isys),       "sigma",                    npt-1, self.arr_pt);
        h1exponential   = TH1F("h1exponential_param{0}".format(isys), "lambda",                   npt-1, self.arr_pt);
        h1offset        = TH1F("h1offset_param{0}".format(isys),      "offset",                   npt-1, self.arr_pt);
        h1linear        = TH1F("h1linear_param{0}".format(isys),      "linear coeff.",            npt-1, self.arr_pt);
        h1quadratic     = TH1F("h1quadratic_param{0}".format(isys),   "quadratic coeff.",         npt-1, self.arr_pt);
        h1cubic         = TH1F("h1cubic_param{0}".format(isys),   "cubic coeff.",         npt-1, self.arr_pt);
        h1FWHM          = TH1F("h1fwhm_param{0}".format(isys),        "fwhm/2.36",                npt-1, self.arr_pt);
        h1yield         = TH1F("h1yield_param{0}".format(isys),       "raw yield",                npt-1, self.arr_pt);
        h1significance  = TH1F("h1significance_param{0}".format(isys),       "significance",                npt-1, self.arr_pt);
        h1statuncert  = TH1F("h1statuncert_param{0}".format(isys),       "stat. uncert.",                npt-1, self.arr_pt);

        h1chi2.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1chi2.SetYTitle("#chi^{2} of fit")
        h1amplitude.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1amplitude.SetYTitle("amplitude of fit");
        h1mean.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1mean.SetYTitle("peak mean (GeV/#it{c}^{2})");
        h1sigma.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1sigma.SetYTitle("peak sigma (GeV/#it{c}^{2})");
        h1exponential.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1exponential.SetYTitle("#exponential coeff. of fit");
        h1offset.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1offset.SetYTitle("offset of fit");
        h1linear.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1linear.SetYTitle("linear coeff. of fit");
        h1quadratic.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1quadratic.SetYTitle("quadratic coeff. of fit");
        h1cubic.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1cubic.SetYTitle("cubic coeff. of fit");
        h1FWHM.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1FWHM.SetYTitle("fwhm/#sqrt{8ln(2)}");
        h1yield.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1yield.SetYTitle("raw yield");
        h1significance.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1significance.SetYTitle("significance");
        h1statuncert.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1statuncert.SetYTitle("stat. uncert. %");

        ##############################################################
        #         get all of the config parameters for calc          #
        ##############################################################

        initial_fit_parameters = info_decay.get_initial_fit_parameters()
        fit_param_limit_min = info_decay.get_fit_param_limit_min()
        fit_param_limit_max = info_decay.get_fit_param_limit_max()

        fit_func = info_decay.get_fit_func()
        fit_range = info_decay.get_fit_range()
        scale_mix_range = info_decay.get_scale_mix_range()
        yield_range = info_decay.get_yield_range()


        ########################################
        #         loop over pt slices          #
        ########################################

        print("Pair_analyzer.py: in the function analyze_ptspectrum going in the loop over pT bins")
        for i in range(info_decay.get_start_pt(), npt-1):
            pt1 = self.arr_pt[i]
            pt2 = self.arr_pt[i+1]

            if len(scale_mix_range) == 1:
                scale_mix_range_for_one_pT = scale_mix_range[0]
            else:
                scale_mix_range_for_one_pT = scale_mix_range[i]

            if self.typ == "data":
                h1prepared, h1mix_scaled, h1same, h1mix = self.slice_histo_and_subtract_background(scale_mix_range_for_one_pT, isys, h2same, h2mix, pt1, pt2)

            if self.typ == "mc":
                h1prepared = self.slice_histo(scale_mix_range_for_one_pT, isys, h2mc, pt1, pt2)
            ###########################
            #         FITTING         #
            ###########################
            print("Pair_analyzer.py: in the function analyze_ptspectrum ######FITTING######", flush=True)
            height = initial_fit_parameters[0]
            mean_init = initial_fit_parameters[1]
            sigma_init = initial_fit_parameters[2]
            lambda_init = initial_fit_parameters[3]
            if "pi0" in self.meson:
                for j in range ( h1prepared.GetXaxis().FindBin(0.125), h1prepared.GetXaxis().FindBin(0.138)):
                    if h1prepared.GetBinContent(j) > height:
                        height = h1prepared.GetBinContent(j);
                s1_bin = h1prepared.GetXaxis().FindBin(0.110)
                s1 = h1prepared.GetBinContent(s1_bin)
                s2_bin = h1prepared.GetXaxis().FindBin(0.150)
                s2 = h1prepared.GetBinContent(s2_bin)
                sub = (s1 + s2)/2
                height_sub = height - sub

            elif "eta" in self.meson:
                for j in range ( h1prepared.GetXaxis().FindBin(0.538), h1prepared.GetXaxis().FindBin(0.552)):
                        if h1prepared.GetBinContent(j) > height:
                            height = h1prepared.GetBinContent(j);
                s1_bin = h1prepared.GetXaxis().FindBin(0.47)
                s1 = h1prepared.GetBinContent(s1_bin)
                s2_bin = h1prepared.GetXaxis().FindBin(0.58)
                s2 = h1prepared.GetBinContent(s2_bin)
                sub = (s1 + s2)/2
                height_sub =  height - sub

            if height < 7: #a.k.a. if the amount of accounts inside the peak is too low, we skip this diagram
                continue

            if len(fit_func) == 1:
                fit_func_for_one_pT = fit_func[0]
            else:
                fit_func_for_one_pT = fit_func[i]
            self.set_fit_function(fit_func_for_one_pT)

            f1total = self.f1total.Clone("f1total_pt_{0}_{1}_syst_{2}_{3}".format(pt1, pt2, isys, fit_func_for_one_pT))
            if len(fit_range) == 1:
                fit_range_for_one_pT = fit_range[0]
            else:
                fit_range_for_one_pT = fit_range[i]
            f1total.SetRange(fit_range_for_one_pT[0], fit_range_for_one_pT[1])
            param_init = [height_sub, mean_init, sigma_init, lambda_init]
            for i_lim, (lower_limit, upper_limit) in enumerate(zip(fit_param_limit_min, fit_param_limit_max)):
                print(i_lim, type(lower_limit))
                f1total.SetParLimits(i_lim, lower_limit*param_init[i_lim], upper_limit*param_init[i_lim])

            f1total.SetParameter(0,height_sub); # amplitude
            f1total.SetParameter(1,mean_init); # mean
            f1total.SetParameter(2,sigma_init); # sigma
            f1total.SetParameter(3,lambda_init); # exponential
            if self.typ =="data":
                f1total.SetParameter(4,50); # offset
                f1total.SetParameter(6,-100); # quadratic

            if initial_fit_parameters[0] == 0: # amplitude 0 -> did not yet assign values in config
                f1total.SetParLimits(0, 0.9*height_sub, 1.2*height_sub);
                f1total.SetParLimits(1, 0.125,0.138);
                f1total.SetParLimits(2, 0.9 * sigma_init, 1.5*sigma_init);
                f1total.SetParLimits(3,0.9 * lambda_init, 1.6*lambda_init);

            h1fit = h1prepared.Clone("h1mgg_fitted_pt_{0}_{1}_syst_{2}_{3}_scale_{4}_{5}_fit_{6}_{7}".format(pt1, pt2, isys,fit_func_for_one_pT, scale_mix_range_for_one_pT[0], scale_mix_range_for_one_pT[1], fit_range_for_one_pT[0], fit_range_for_one_pT[1]));
            h1fit.SetDirectory(0)
            #h1fit.Fit(f1total,"RLEM","",info_decay.get_fit_min(ir), info_decay.get_fit_max(ir)); #RME or RLEM
            fitresult = h1fit.Fit(f1total,"RQS","",fit_range_for_one_pT[0], fit_range_for_one_pT[1]); #RME or RLEM

            print(fitresult)

            amplitude       = f1total.GetParameter(0);
            amplitude_err   = f1total.GetParError(0)
            mean            = f1total.GetParameter(1);
            mean_err        = f1total.GetParError(1);
            sigma           = f1total.GetParameter(2);
            sigma_err       = f1total.GetParError(2);
            exponential     = f1total.GetParameter(3);
            exponential_err = f1total.GetParError(3);

            if self.typ == "data":
                offset          = f1total.GetParameter(4);
                offset_err      = f1total.GetParError(4);
                linear          = f1total.GetParameter(5);
                linear_err      = f1total.GetParError(5);
                if self.func == "GaussExpQuadratic" or self.func == "GaussExpCube":
                    quadratic = f1total.GetParameter(6);
                    quadratic_err = f1total.GetParError(6);
                    if self.func == "GaussExpCube":
                        cubic = f1total.GetParameter(7);
                        cubic_err = f1total.GetParError(7);
                params = [amplitude, mean, sigma, exponential, offset, linear]
                params_err = [amplitude_err, mean_err, sigma_err, exponential_err, offset_err, linear_err]
                if self.func == "GaussExpQuadratic":
                    params = [amplitude, mean, sigma, exponential, offset, linear, quadratic]
                    params_err = [amplitude_err, mean_err, sigma_err, exponential_err, offset_err, linear_err, quadratic_err]
                if self.func == "GaussExpCube":
                    params = [amplitude, mean, sigma, exponential, offset, linear, quadratic, cubic]
                    params_err = [amplitude_err, mean_err, sigma_err, exponential_err, offset_err, linear_err, quadratic_err, cubic_err]
                tf1_fwhm, FWHM, FWHM_err  = self.calc_FWHM(f1total, params, params_err)
            elif self.typ == "mc":
                params = [amplitude, mean, sigma, exponential]
                params_err = [amplitude_err, mean_err, sigma_err, exponential_err]
                tf1_fwhm, FWHM, FWHM_err  = self.calc_FWHM_MC(f1total, params, params_err)
                print("FWHM", FWHM, type(FWHM))

            tf1_fwhm.SetName("f1fwhm_pt{0}_{1}_syst{2}".format(pt1, pt2, isys));

            if len(yield_range) == 1:
                yield_range_for_one_pT = yield_range[0]
            else:
                yield_range_for_one_pT = yield_range[i]
            print("parameters: ", params[0], " ", params[1], " ", params[2], " ", params[3])
            print("cent1 = ", self.cent1, " cent2 = ", self.cent2)
            print("pt1 = ", pt1, " pt2 = ", pt2, "npt = ", i)
            print("scaling of mix inv mass min = ", scale_mix_range_for_one_pT[0], " max = ", scale_mix_range_for_one_pT[1])
            print("fit function: ", fit_func_for_one_pT)
            print("fit range min = ", fit_range_for_one_pT[0], " max = ", fit_range_for_one_pT[1])
            print("raw yield range min = ", yield_range_for_one_pT[0], " max = ", yield_range_for_one_pT[1])

            if self.typ == "data":
                raw_yield, error_raw_yield, signal, stat_uncert = self.calculate_raw_yield(yield_range_for_one_pT, h1prepared, params, params_err, (pt2-pt1))
                print("Typ von raw_yield: ", type(raw_yield))
                significane = self.calculate_significance(h1same, signal, yield_range_for_one_pT, params)

            elif self.typ == "mc":
                raw_yield, error_raw_yield, stat_uncert = self.calculate_raw_yield_MC(yield_range_for_one_pT, h1prepared, params, params_err, (pt2-pt1))

                print("################Comparison of uncertainties, raw yield monte carlo:")
                print("value of raw yield: ", raw_yield, "±", error_raw_yield )

            dof = fitresult.Ndf()
            if dof != 0:
                chi2 = fitresult.Chi2()
                chi2_dof = chi2/dof
                h1chi2.SetBinContent(i+1, chi2_dof)
            h1chi2.SetBinError(i+1, 0)
            h1amplitude.SetBinContent(i+1, amplitude)
            h1amplitude.SetBinError(i+1, amplitude_err)
            h1mean.SetBinContent(i+1, mean);
            h1mean.SetBinError(i+1, mean_err);
            h1sigma.SetBinContent(i+1, sigma);
            h1sigma.SetBinError(i+1, sigma_err);
            h1exponential.SetBinContent(i+1, exponential);
            h1exponential.SetBinError(i+1, exponential_err);
            if not math.isnan(FWHM):
                print("is not nan", FWHM)
                h1FWHM.SetBinContent(i+1, FWHM)
                h1FWHM.SetBinError(i+1, FWHM_err);
            else:
                print("is nan", FWHM)
            if self.typ == "data":
                h1linear.SetBinContent(i+1, linear);
                h1linear.SetBinError(i+1, linear_err);
                h1offset.SetBinContent(i+1, offset);
                h1offset.SetBinError(i+1, offset_err);
                h1significance.SetBinContent(i+1, significane)
                h1significance.SetBinError(i+1, 0)

                if self.func == "GaussExpQuadratic" or self.func == "GaussExpCube":
                    h1quadratic.SetBinContent(i+1, quadratic);
                    h1quadratic.SetBinError(i+1, quadratic_err);
                    if self.func == "GaussExpCube":
                        h1cubic.SetBinContent(i+1, cubic);
                        h1cubic.SetBinError(i+1, cubic_err);
            h1yield.SetBinContent(i+1, raw_yield)
            h1yield.SetBinError(i+1, error_raw_yield)
            if stat_uncert!=-100:
                h1statuncert.SetBinContent(i+1, stat_uncert)
                h1statuncert.SetBinError(i+1, 0)

            if self.typ == "data":
                outlist.append(h1mix_scaled);
                outlist.append(h1same);
                outlist.append(h1mix);
            outlist.append(h1prepared);
            outlist.append(h1fit);
            outlist.append(f1total);
            print("\n \n \n")

        h1yield.SetDirectory(0)
        h1mean.SetDirectory(0)
        h1exponential.SetDirectory(0)
        h1sigma.SetDirectory(0)
        h1FWHM.SetDirectory(0)
        h1chi2.SetDirectory(0)
        h1amplitude.SetDirectory(0)
        h1linear.SetDirectory(0)
        h1offset.SetDirectory(0)
        h1quadratic.SetDirectory(0)
        if self.typ == "data":
            h1significance.SetDirectory(0)
        h1statuncert.SetDirectory(0)
        h1statuncert.Scale(100.)

        outlist.append(h1yield);
        outlist.append(h1mean);
        outlist.append(h1exponential);
        outlist.append(h1sigma);
        outlist.append(h1FWHM);
        outlist.append(h1chi2)
        outlist.append(h1amplitude);
        outlist.append(h1linear);
        outlist.append(h1offset);
        outlist.append(h1quadratic);
        if self.typ == "data":
            outlist.append(h1significance)
        outlist.append(h1statuncert)

        del h1amplitude
        del h1mean
        del tf1_fwhm
        del h1ev
        del h1exponential
        del h1FWHM
        del h1sigma
        del h1linear
        del h1offset
        del h1quadratic
        if self.typ == "data":
            del h1mix
            del h1same
            del h1significance
        del h1prepared
        del h1fit
        del f1total
        del h1yield
        del h1statuncert
        print("Pair_analyzer.py the function analyze_ptspectrum is done for this subsystem")

        return outlist;
