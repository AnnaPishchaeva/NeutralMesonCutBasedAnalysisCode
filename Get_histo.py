# This code was written by Julia Schlägel (July 2024)

import ROOT
import yaml
import os
import pickle
from ctypes import *
import numpy as np
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency, TString
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

class Get_histo: #in this class all central things are done that are needed in multiple macros like loading a file from a rootfile,...
    def __init__(self, info_decay):
        print("/////////////////////////////////////////////////")
        print("Get_histo.py")
        print("/////////////////////////////////////////////////")
        self.yield_list = []
        self.npt = len(info_decay.get_arr_pt())
        self.arr_pt = info_decay.get_arr_pt()
        self.start_pt = info_decay.get_start_pt()
        self.list_parameters_all_cent = []

    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input root file.");
            self.rootfile.Close();

    def set_input_file(self, filename, info_decay):
        self.rootfile = TFile.Open(filename, "READ")
        self.list_plot = self.rootfile.Get(info_decay.get_decay())

    def set_input_pkl_file(self, inname):
        self.filename_pkl = inname

    def set_list_raw_yield_all_cent(self, isys):
        self.yield_list.append(self.get_yield_one_cent(isys))

    def set_list_parameters_all_cent(self):
        parameter_list = self.get_parameter_list()
        parameter_list_cut = self.get_list_parameters_comparison(parameter_list)
        self.list_parameters_all_cent.append(parameter_list_cut)

    def get_list_raw_yield_all_cent(self):
        return self.yield_list

    def get_list_parameters_all_cent(self):
        return self.list_parameters_all_cent

    def get_list_parameters_comparison(self, parameter_list):
        parameter_comparison = []
        for i_parameter in range(6):
            parameter_comparison.append(parameter_list[i_parameter]);
        return parameter_comparison

    def get_parameter_list(self):
        self.parameter_list = []
        list_keys = self.list_plot.GetListOfKeys()
        for key in list_keys:
            histo = key.ReadObj()
            histo_name = histo.GetName()
            if "param" in histo_name:
                histo.SetDirectory(0)
                self.parameter_list.append(histo)
        return self.parameter_list

    def get_yield_list_one_cent(self, isys):
        self.yield_list_one_cent = []
        histo = self.list_plot.Get("h1yield_param{0}".format(isys))
        histo.SetDirectory(0)
        self.yield_list_one_cent.append(histo)
        return self.yield_list_one_cent

    def get_yield_one_cent(self, isys):
        histo = self.list_plot.Get("h1yield_param{0}".format(isys))
        histo.SetDirectory(0)
        return histo

    def get_significance_one_cent(self, isys):
        histo = self.list_plot.Get("h1significance_param{0}".format(isys))
        histo.SetDirectory(0)
        return histo

    def get_significance_list_one_cent(self, isys):
        self.significance_list_one_cent = []
        histo = self.list_plot.Get("h1significance_param{0}".format(isys))
        histo.SetDirectory(0)
        self.significance_list_one_cent.append(histo)
        return self.significance_list_one_cent

    def get_stat_uncert_one_cent(self, isys):
        histo = self.list_plot.Get("h1statuncert_param{0}".format(isys))
        histo.SetDirectory(0)
        return histo

    def get_stat_uncert_list_one_cent(self, isys):
        self.stat_uncert_list_one_cent = []
        histo = self.list_plot.Get("h1statuncert_param{0}".format(isys))
        histo.SetDirectory(0)
        self.stat_uncert_list_one_cent.append(histo)
        return self.stat_uncert_list_one_cent

    def get_yield_list(self, nsys):
        self.yield_list = []
        for isys in range(0,nsys):
            histo = self.get_yield_one_cent(isys)
            print(isys)
            if histo:
                print("histo isys: ", isys, flush=True)
            self.yield_list.append(histo)
        return self.yield_list

    def get_function_list(self, isys, fit_func):
        self.function_list = []
        for ipt in range(self.start_pt, self.npt-1):
            pt1 = self.arr_pt[ipt]
            pt2 = self.arr_pt[ipt+1]
            if len(fit_func) == 1:
                fit_func_for_one_pT = fit_func[0]
            else:
                fit_func_for_one_pT = fit_func[ipt]
            histo = self.list_plot.Get("f1total_pt_{0}_{1}_syst_{2}_{3}".format(pt1, pt2, isys, fit_func_for_one_pT))
            if histo:
                self.function_list.append(histo)
            else:
                continue
        return self.function_list

    def get_subtracted_list(self, isys, scale_mix_range):
        self.subtracted_list = []
        for ipt in range(self.start_pt, self.npt-1):
            pt1 = self.arr_pt[ipt]
            pt2 = self.arr_pt[ipt+1]
            if len(scale_mix_range) == 1:
                scale_mix_range_for_one_pT = scale_mix_range[0]
            else:
                scale_mix_range_for_one_pT = scale_mix_range[ipt]

            histo = self.list_plot.Get("h1mgg_pt_{0}_{1}_syst_{2}_scale_{3}_{4}".format(pt1, pt2,isys, scale_mix_range_for_one_pT[0], scale_mix_range_for_one_pT[1]))
            if histo:
                self.subtracted_list.append(histo)
            else:
                continue

        return self.subtracted_list

    def get_same_list(self, isys):
        self.same_list = []
        for ipt in range(self.start_pt, self.npt-1):
            pt1 = self.arr_pt[ipt]
            pt2 = self.arr_pt[ipt+1]
            histo = self.list_plot.Get("h1mgg_same_pt_{0}_{1}_syst_{2}".format(pt1, pt2, isys))
            if histo:
                self.same_list.append(histo)
            else:
                continue
        return self.same_list

    def get_mixed_list(self, isys, scale_mix_range):
        self.mixed_list = []
        for ipt in range(self.start_pt, self.npt-1):
            pt1 = self.arr_pt[ipt]
            pt2 = self.arr_pt[ipt+1]
            if len(scale_mix_range) == 1:
                scale_mix_range_for_one_pT = scale_mix_range[0]
            else:
                scale_mix_range_for_one_pT = scale_mix_range[ipt]
            histo = self.list_plot.Get("h1mgg_mix_scaled_pt_{0}_{1}_syst_{2}_scale_{3}_{4}".format(pt1, pt2, isys, scale_mix_range_for_one_pT[0], scale_mix_range_for_one_pT[1]))
            if histo:
                self.mixed_list.append(histo)
            else:
                continue
        return self.mixed_list

    def get_histo_from_pkl(self):
        if os.path.exists(self.filename_pkl): #if in the current directory this file exists
            with open(self.filename_pkl, 'rb') as raw_file:
                raw = pickle.load(raw_file) #Retrieving the data (raw yield) that we saved before
        return raw

    def get_histo_from_other_file(self, hist_name, path_comp):
        #if meson == "pi0":
        #    dir_1 = self.rootfile.Get("Pi013TeV")
        #if meson == "eta":
        #    dir_1 = self.rootfile.Get("Eta13TeV")
        #hist = dir_1.Get(hist_name)
        if path_comp != "":
            dir_1 = self.rootfile.Get(path_comp)
            hist = self.dir_1.Get(hist_name)
        else:
            hist = self.rootfile.Get(hist_name)
        return hist

    def get_eff_acc_BR_from_other_file(self, hist_name, path_comp):
        if path_comp != "":
            dir_1 = self.rootfile.Get(path_comp)
            acc = self.dir_1.Get(hist_name)
            eff = self.dir_1.Get("TrueMesonEffiPt")
        else:
            acc = self.rootfile.Get(hist_name)
            eff = self.rootfile.Get("TrueMesonEffiPt")
        eff_acc_BR = acc.Clone("acc")
        acc.Sumw2()
        eff.Sumw2()
        length = acc.GetNbinsX()
        for ibin in range(2,length+1):
            bin_acc = acc.GetBinContent(ibin)
            print("ibin: ", ibin, "bin_acc: ", bin_acc)
            bin_eff = eff.GetBinContent(ibin)
            err_acc = acc.GetBinError(ibin)/bin_acc
            err_eff = eff.GetBinError(ibin)/bin_eff
            error = bin_acc*bin_eff*np.sqrt(err_acc*err_acc + err_eff*err_eff)
            print("error: ", error)
            eff_acc_BR.SetBinContent(ibin, bin_acc*bin_eff*0.988)
            eff_acc_BR.SetBinError(ibin, error)
        return eff_acc_BR

    def get_corr_yield_from_other_file(self, path_comp, isys, bin_var):
        #pt_bin_other = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.1, 1.2, 1.3, 1.4,
        #1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0,
        #9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0]
        dir_1 = self.rootfile.Get(path_comp)
        if isys == 0:
            histo1 = dir_1.Get("Hist1D_y1")
            histo2 = dir_1.Get("Hist1D_y2")
            histo_error_stat_1 = dir_1.Get("Hist1D_y1_e1")
            histo_error_stat_2 = dir_1.Get("Hist1D_y2_e1")
            histo_error_syst_1 = dir_1.Get("Hist1D_y1_e2")
            histo_error_syst_2 = dir_1.Get("Hist1D_y2_e2")
            histo_error_syst_uncorr_1 = dir_1.Get("Hist1D_y1_e3")
            histo_error_syst_uncorr_2 = dir_1.Get("Hist1D_y2_e3")
            corr1 = histo1.Clone("corr1")
            corr2 = histo2.Clone("corr2")

        else:
            histo1 = dir_1.Get("Hist1D_y{0}".format(isys+2))
            histo_error_stat_1 = dir_1.Get("Hist1D_y{0}_e1".format(isys+2))
            histo_error_syst_1 = dir_1.Get("Hist1D_y{0}_e2".format(isys+2))
            histo_error_syst_uncorr_1 = dir_1.Get("Hist1D_y{0}_e3".format(isys+2))
            corr1 = histo1.Clone("corr1")

        nbin = corr1.GetNbinsX()+1
        if isys == 0:
            for ibin in range(nbin):
                content1 = corr1.GetBinContent(ibin)
                content2 = corr2.GetBinContent(ibin)
                content1 = (content1+content2)/2
                corr1.SetBinContent(ibin, content1)
                error_stat_1 = histo_error_stat_1.GetBinContent(ibin) + histo_error_syst_uncorr_1.GetBinContent(ibin)
                error_stat_2 = histo_error_stat_2.GetBinContent(ibin) + histo_error_syst_uncorr_2.GetBinContent(ibin)
                error_syst_1 = histo_error_syst_1.GetBinContent(ibin)
                error_syst_2 = histo_error_syst_2.GetBinContent(ibin)
                error_stat = np.sqrt(error_stat_1*error_stat_1 + error_stat_2*error_stat_2)/2
                error_syst = (error_syst_1+error_syst_2)/2
                error = np.sqrt(error_stat*error_stat+error_syst*error_syst)
                corr1.SetBinError(ibin, error)
                print(histo_error_stat_1.GetBinContent(ibin),"+", histo_error_syst_uncorr_1.GetBinContent(ibin), " error_stat1: ", error_stat_1)
                print(histo_error_stat_2.GetBinContent(ibin),"+", histo_error_syst_uncorr_2.GetBinContent(ibin), " error_stat2: ", error_stat_2)
                print("error_stat: ", error_stat)
                print("error_syst1: ", error_syst_1, " error_syst_2: ", error_syst_2, "error_syst: ", error_syst)
                print("error: ", error)

        else:
            for ibin in range(nbin):
                error_stat = histo_error_stat_1.GetBinContent(ibin) + histo_error_syst_uncorr_1.GetBinContent(ibin)
                error_syst = histo_error_syst_1.GetBinContent(ibin)
                error = np.sqrt(error_stat*error_stat+error_syst*error_syst)
                corr1.SetBinError(ibin, error)

        new_hist = TH1D("p_T1", "p_T1", len(bin_var) - 1, np.array(bin_var))
        new_h_err_stat = TH1D("err_stat", "err_stat", len(bin_var) - 1, np.array(bin_var))
        new_h_err_syst = TH1D("err_syst", "err_syst", len(bin_var) - 1, np.array(bin_var))
        new_hist.Sumw2()
        xaxis = corr1.GetXaxis()
        hist_notrebin = corr1.Clone("hist_notrebin")

        for bin in range(0, len(bin_var) - 1):
            b1 = bin_var[bin]
            b2 = bin_var[bin + 1] #basically value of x-axis
            bin_b1 = corr1.GetXaxis().FindBin(b1 + 1e-6) #the bin nmber of the value of x-axis
            bin_b2 = corr1.GetXaxis().FindBin(b2 - 1e-6)

            error = c_double(0.0)
            # content = hist.IntegralAndError(bin_b1, bin_b2, error, "")
            content = 0
            error_stat1 = 0
            error_stat2 = 0
            error_syst_1 = 0
            error_syst_2 = 0
            n_syst = 0
            print("bin1: ", bin_b1, bin_b2)
            for i in range(bin_b1, bin_b2 +1):
                n_syst +=1
                content = content+ corr1.GetBinContent(i)*(xaxis.GetBinUpEdge(i)-xaxis.GetBinLowEdge(i))
                #content = content+ corr1.GetBinContent(i)
                error_stat1 += (histo_error_stat_1.GetBinContent(i) + histo_error_syst_uncorr_1.GetBinContent(i))*(histo_error_stat_1.GetBinContent(i) + histo_error_syst_uncorr_1.GetBinContent(i)) #e_pt1*e_pt1+e_pt2*e_pt2+...
                error_syst_1 += histo_error_syst_1.GetBinContent(i) #e_pt1+e_pt2+...
                if isys == 0:
                    error_stat2 += (histo_error_stat_2.GetBinContent(i) + histo_error_syst_uncorr_2.GetBinContent(i))*(histo_error_stat_2.GetBinContent(i)+ histo_error_syst_uncorr_2.GetBinContent(i))
                    error_syst_2 += histo_error_syst_2.GetBinContent(i)
            if isys == 0:
                error_stat = np.sqrt(error_stat1+error_stat2)
                error_syst = (error_syst_1+error_syst_2)/(2*n_syst)
            else:
                error_stat = np.sqrt(error_stat1) #sqrt(e_pt1*e_pt1+e_pt2*e_pt2+...)
                error_syst = (error_syst_1)/(n_syst) #(e_pt1+e_pt2+..)/n_pt
            error = np.sqrt(error_stat*error_stat+error_syst*error_syst) #sqrt(e_stat*e_stat+e_syst*e_syst)

            #error = np.sqrt(content)

            bin_center = (b1 + b2) / 2  # Take the center of the bin to fill
            bin_width = b2-b1 #basically the width of one bin
            print("in rebin function error:", error)
            print("content: ", content)
            #
            new_content = content / bin_width
            # new_hist.Fill(bin_center, new_content)  # Filling with this center of each bin
            new_hist.SetBinContent(new_hist.FindBin(bin_center), new_content)
            new_h_err_stat.SetBinContent(new_h_err_stat.FindBin(bin_center), error_stat)
            new_h_err_syst.SetBinContent(new_h_err_stat.FindBin(bin_center), error_syst)
            #new_hist.SetBinError(bin, error) # / bin_width)  # calculate error to the right bin#
            new_hist.SetBinError(new_hist.FindBin(bin_center), error)

        if new_hist is None:
            print(f"Failed to load histogram: {hist_name}")

        new_hist.Scale(0.5)
        hist_notrebin.Scale(0.5)
        #corr11.Rebin(2)

        return hist_notrebin, new_hist, new_h_err_stat, new_h_err_syst

    def get_corr_yield_from_other_file_cent20(self, path_comp, isys, bin_var):
        #pt_bin_other = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.1, 1.2, 1.3, 1.4,
        #1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0,
        #9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0]
        dir_1 = self.rootfile.Get(path_comp)
        if isys == 0:
            histo1 = dir_1.Get("Hist1D_y1")
            histo2 = dir_1.Get("Hist1D_y2")
            histo3 = dir_1.Get("Hist1D_y3")
            histo_error_stat_1 = dir_1.Get("Hist1D_y1_e1")
            histo_error_stat_2 = dir_1.Get("Hist1D_y2_e1")
            histo_error_stat_3 = dir_1.Get("Hist1D_y3_e1")
            histo_error_syst_1 = dir_1.Get("Hist1D_y1_e2")
            histo_error_syst_2 = dir_1.Get("Hist1D_y2_e2")
            histo_error_syst_3 = dir_1.Get("Hist1D_y3_e2")
            corr3 = histo3.Clone("corr3")
        if isys != 0 and isys != 4:
            histo1 = dir_1.Get("Hist1D_y{0}".format(isys+isys+2))
            histo2 = dir_1.Get("Hist1D_y{0}".format(isys+isys+3))
            histo_error_stat_1 = dir_1.Get("Hist1D_y{0}_e1".format(isys+isys+2))
            histo_error_stat_2 = dir_1.Get("Hist1D_y{0}_e1".format(isys+isys+3))
            histo_error_syst_1 = dir_1.Get("Hist1D_y{0}_e2".format(isys+isys+2))
            histo_error_syst_2 = dir_1.Get("Hist1D_y{0}_e2".format(isys+isys+3))
            print("isys_file: ", isys)

        if isys != 4:
            corr1 = histo1.Clone("corr1")
            corr2 = histo2.Clone("corr2")
        else:
            corr1 = dir_1.Get("Hist1D_y10")


        if isys !=4:
            nbin = corr1.GetNbinsX()+1
            for ibin in range(nbin):
                content1 = corr1.GetBinContent(ibin)
                content2 = corr2.GetBinContent(ibin)
                content1 = (content1+content2)/2
                corr1.SetBinContent(ibin, content1)
                if isys!=0:
                    error_stat_1 = histo_error_stat_1.GetBinContent(ibin)
                    error_stat_2 = histo_error_stat_2.GetBinContent(ibin)
                    error_syst_1 = histo_error_syst_1.GetBinContent(ibin)
                    error_syst_2 = histo_error_syst_2.GetBinContent(ibin)
                    error_stat = np.sqrt(error_stat_1*error_stat_1 + error_stat_2*error_stat_2)/2
                    error_syst = (error_syst_1+error_syst_2)/2
                    error = np.sqrt(error_stat*error_stat+error_syst*error_syst)
                    corr1.SetBinError(ibin, error)
            if isys == 0:
                nbin = corr1.GetNbinsX()+1
                for ibin in range(nbin):
                    content1 = corr1.GetBinContent(ibin)
                    content3 = corr3.GetBinContent(ibin)
                    content1 = (content1+content3)/2
                    corr1.SetBinContent(ibin, content1)
                    error_stat_1 = histo_error_stat_1.GetBinContent(ibin)
                    error_stat_2 = histo_error_stat_2.GetBinContent(ibin)
                    error_stat_3 = histo_error_stat_3.GetBinContent(ibin)
                    error_syst_1 = histo_error_syst_1.GetBinContent(ibin)
                    error_syst_2 = histo_error_syst_2.GetBinContent(ibin)
                    error_syst_3 = histo_error_syst_3.GetBinContent(ibin)
                    error_stat = np.sqrt( (error_stat_1*error_stat_1 + error_stat_2*error_stat_2)/4 + error_stat_3*error_stat_3)/2
                    error_syst = ( (error_syst_1+error_syst_2)/2 + error_syst_3)/2
                    error = np.sqrt(error_stat*error_stat+error_syst*error_syst)
                    corr1.SetBinError(ibin, error)

            else:
                for ibin in range(nbin):
                    error_stat = histo_error_stat_1.GetBinContent(ibin)
                    error_syst = histo_error_syst_1.GetBinContent(ibin)
                    error = np.sqrt(error_stat*error_stat+error_syst*error_syst)
                    corr1.SetBinError(ibin, error)

        new_hist = TH1D("p_T1", "p_T1", len(bin_var) - 1, np.array(bin_var))
        new_h_err_stat = TH1D("err_stat", "err_stat", len(bin_var) - 1, np.array(bin_var))
        new_h_err_syst = TH1D("err_syst", "err_syst", len(bin_var) - 1, np.array(bin_var))
        new_hist.Sumw2()
        xaxis = corr1.GetXaxis()
        hist_notrebin = corr1.Clone("hist_notrebin")

        for bin in range(0, len(bin_var) - 1):
            b1 = bin_var[bin]
            b2 = bin_var[bin + 1] #basically value of x-axis
            bin_b1 = corr1.GetXaxis().FindBin(b1 + 1e-6) #the bin nmber of the value of x-axis
            bin_b2 = corr1.GetXaxis().FindBin(b2 - 1e-6)

            error = c_double(0.0)
            # content = hist.IntegralAndError(bin_b1, bin_b2, error, "")
            content = 0
            error_stat1 = 0
            error_stat2 = 0
            error_syst_1 = 0
            error_syst_2 = 0
            n_syst = 0
            print("bin1: ", bin_b1, bin_b2)
            for i in range(bin_b1, bin_b2 +1):
                n_syst +=1
                content = content+ corr1.GetBinContent(i)*(xaxis.GetBinUpEdge(i)-xaxis.GetBinLowEdge(i))
                #content = content+ corr1.GetBinContent(i)
                error_stat1 += histo_error_stat_1.GetBinContent(i)*histo_error_stat_1.GetBinContent(i) #e_pt1*e_pt1+e_pt2*e_pt2+...
                error_syst_1 += histo_error_syst_1.GetBinContent(i) #e_pt1+e_pt2+...
                if isys == 0:
                    error_stat2 += histo_error_stat_2.GetBinContent(i)*histo_error_stat_2.GetBinContent(i)
                    error_syst_2 += histo_error_syst_2.GetBinContent(i)
            if isys == 0:
                error_stat = np.sqrt(error_stat1+error_stat2)
                error_syst = (error_syst_1+error_syst_2)/(2*n_syst)
            else:
                error_stat = np.sqrt(error_stat1) #sqrt(e_pt1*e_pt1+e_pt2*e_pt2+...)
                error_syst = (error_syst_1)/(n_syst) #(e_pt1+e_pt2+..)/n_pt
            error = np.sqrt(error_stat*error_stat+error_syst*error_syst) #sqrt(e_stat*e_stat+e_syst*e_syst)

            #error = np.sqrt(content)

            bin_center = (b1 + b2) / 2  # Take the center of the bin to fill
            bin_width = b2-b1 #basically the width of one bin
            print("in rebin function error:", error)
            print("content: ", content)
            #
            new_content = content / bin_width
            # new_hist.Fill(bin_center, new_content)  # Filling with this center of each bin
            new_hist.SetBinContent(new_hist.FindBin(bin_center), new_content)
            new_h_err_stat.SetBinContent(new_h_err_stat.FindBin(bin_center), error_stat)
            new_h_err_syst.SetBinContent(new_h_err_stat.FindBin(bin_center), error_syst)
            #new_hist.SetBinError(bin, error) # / bin_width)  # calculate error to the right bin#
            new_hist.SetBinError(new_hist.FindBin(bin_center), error)

        if new_hist is None:
            print(f"Failed to load histogram: {hist_name}")

        new_hist.Scale(0.5)
        hist_notrebin.Scale(0.5)
        #corr11.Rebin(2)

        return hist_notrebin, new_hist, new_h_err_stat, new_h_err_syst

    def get_corr_yield_from_other_file_cent20_another(self, path_comp, isys, bin_var):
        #pt_bin_other = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.1, 1.2, 1.3, 1.4,
        #1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0,
        #9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0]
        dir_1 = self.rootfile.Get(path_comp)
        hist_list = []
        error_stat_list = []
        error_syst_list = []

        if isys == 0:
            histo1 = dir_1.Get("Hist1D_y1")
            histo2 = dir_1.Get("Hist1D_y2")
            histo3 = dir_1.Get("Hist1D_y3")
            histo_error_stat_1 = dir_1.Get("Hist1D_y1_e1")
            histo_error_stat_2 = dir_1.Get("Hist1D_y2_e1")
            histo_error_stat_3 = dir_1.Get("Hist1D_y3_e1")
            histo_error_syst_uncorr_1 = dir_1.Get("Hist1D_y1_e3")
            error_stat_list.append(histo_error_stat_1 + histo_error_syst_uncorr_1)
            histo_error_syst_uncorr_2 = dir_1.Get("Hist1D_y2_e3")
            error_stat_list.append(histo_error_stat_2 + histo_error_syst_uncorr_2)
            histo_error_syst_uncorr_3 = dir_1.Get("Hist1D_y3_e3")
            error_stat_list.append(histo_error_stat_3 + histo_error_syst_uncorr_3)
            histo_error_syst_1 = dir_1.Get("Hist1D_y1_e2")
            error_syst_list.append(histo_error_syst_1)
            histo_error_syst_2 = dir_1.Get("Hist1D_y2_e2")
            error_syst_list.append(histo_error_syst_2)
            histo_error_syst_3 = dir_1.Get("Hist1D_y3_e2")
            error_syst_list.append(histo_error_syst_3)
            corr3 = histo3.Clone("corr3")
            hist_list.append(corr3)
        if isys != 0 and isys != 4:
            histo1 = dir_1.Get("Hist1D_y{0}".format(isys+isys+2))
            histo2 = dir_1.Get("Hist1D_y{0}".format(isys+isys+3))
            histo_error_stat_1 = dir_1.Get("Hist1D_y{0}_e1".format(isys+isys+2))
            histo_error_stat_2 = dir_1.Get("Hist1D_y{0}_e1".format(isys+isys+3))
            histo_error_syst_uncorr_1 = dir_1.Get("Hist1D_y{0}_e3".format(isys+isys+2))
            error_stat_list.append(histo_error_stat_1 + histo_error_syst_uncorr_1)
            histo_error_syst_uncorr_2 = dir_1.Get("Hist1D_y{0}_e3".format(isys+isys+3))
            error_stat_list.append(histo_error_stat_2 + histo_error_syst_uncorr_2)
            histo_error_syst_1 = dir_1.Get("Hist1D_y{0}_e2".format(isys+isys+2))
            error_syst_list.append(histo_error_syst_1)
            histo_error_syst_2 = dir_1.Get("Hist1D_y{0}_e2".format(isys+isys+3))
            error_syst_list.append(histo_error_syst_2)
            print("isys_file: ", isys)

        if isys != 4:
            corr1 = histo1.Clone("corr1")
            corr2 = histo2.Clone("corr2")
            hist_list.append(corr1)
            hist_list.append(corr2)
        else:
            corr1 = dir_1.Get("Hist1D_y10")
            hist_list.append(corr1)

        new_hist_list = [0]*len(hist_list)
        new_h_err_stat_list = [0]*len(error_stat_list)
        new_h_err_syst_list = [0]*len(error_syst_list)
        for i in range(len(hist_list)):
            new_hist_list[i] = TH1D("p_T1_{0}".format(i), "p_T1_{0}".format(i), len(bin_var) - 1, np.array(bin_var))
            new_hist_list[i].Sumw2()
            if error_syst_list and error_stat_list:
                new_h_err_stat_list[i] = TH1D("err_stat_{0}".format(i), "err_stat_{0}".format(i), len(bin_var) - 1, np.array(bin_var))
                new_h_err_syst_list[i] = TH1D("err_syst_{0}".format(i), "err_syst_{0}".format(i), len(bin_var) - 1, np.array(bin_var))
        xaxis = corr1.GetXaxis()
        hist_notrebin = corr1.Clone("hist_notrebin")

        for ihist in range(len(new_hist_list)):
            for bin in range(0, len(bin_var) - 1):
                b1 = bin_var[bin]
                b2 = bin_var[bin + 1] #basically value of x-axis
                bin_b1 = hist_list[ihist].GetXaxis().FindBin(b1 + 1e-6) #the bin nmber of the value of x-axis
                bin_b2 = hist_list[ihist].GetXaxis().FindBin(b2 - 1e-6)

                error = c_double(0.0)
                # content = hist.IntegralAndError(bin_b1, bin_b2, error, "")
                content = 0
                error_stat_1 = 0
                error_stat_2 = 0
                error_syst_1 = 0
                error_syst_2 = 0
                n_syst = 0
                print("bin1: ", bin_b1, bin_b2)
                for i in range(bin_b1, bin_b2 +1):
                    n_syst +=1
                    content = content+ hist_list[ihist].GetBinContent(i)*(xaxis.GetBinUpEdge(i)-xaxis.GetBinLowEdge(i))
                    #content = content+ corr1.GetBinContent(i)
                    error_stat_1 += error_stat_list[ihist].GetBinContent(i)*error_stat_list[ihist].GetBinContent(i) #e_pt1*e_pt1+e_pt2*e_pt2+...
                    error_syst_1 += error_syst_list[ihist].GetBinContent(i) #e_pt1+e_pt2+...
                    error_stat = np.sqrt(error_stat_1) #sqrt(e_pt1*e_pt1+e_pt2*e_pt2+...)
                    error_syst = (error_syst_1)/(n_syst) #(e_pt1+e_pt2+..)/n_pt
                error = np.sqrt(error_stat*error_stat+error_syst*error_syst) #sqrt(e_stat*e_stat+e_syst*e_syst)

                bin_center = (b1 + b2) / 2  # Take the center of the bin to fill
                bin_width = b2-b1 #basically the width of one bin
                print("in rebin function error:", error)
                print("content: ", content)
                #
                new_content = content / bin_width
                # new_hist.Fill(bin_center, new_content)  # Filling with this center of each bin
                new_hist_list[ihist].SetBinContent(new_hist_list[ihist].FindBin(bin_center), new_content)
                new_h_err_stat_list[ihist].SetBinContent(new_h_err_stat_list[ihist].FindBin(bin_center), error_stat)
                new_h_err_syst_list[ihist].SetBinContent(new_h_err_syst_list[ihist].FindBin(bin_center), error_syst)
                #new_hist.SetBinError(bin, error) # / bin_width)  # calculate error to the right bin#
                new_hist_list[ihist].SetBinError(new_hist_list[ihist].FindBin(bin_center), error)

        #if new_hist is None:
        #    print(f"Failed to load histogram: {hist_name}")

        new_cent_20_hist = new_hist_list[0].Clone("new_cent_20_hist")
        new_cent_20_err_stat = new_h_err_stat_list[0].Clone("new_cent_20_err_stat")
        new_cent_20_err_syst = new_h_err_syst_list[0].Clone("new_cent_20_err_syst")
        nbin = new_hist_list[0].GetNbinsX()+1
        for ibin in range(nbin):
            content = 0
            error_stat = 0
            error_syst = 0
            for ihist in range(len(new_hist_list)):
                if len(new_hist_list) < 3:
                    content += new_hist_list[ihist].GetBinContent(ibin)
                    error_stat += new_h_err_stat_list[ihist].GetBinContent(ibin)*new_h_err_stat_list[ihist].GetBinContent(ibin)
                    error_syst += new_h_err_syst_list[ihist].GetBinContent(ibin)
                else:
                    if ihist < len(new_hist_list)-1:
                        content += new_hist_list[ihist].GetBinContent(ibin)/2
                        error_stat += new_h_err_stat_list[ihist].GetBinContent(ibin)*new_h_err_stat_list[ihist].GetBinContent(ibin)/4
                        error_syst += new_h_err_syst_list[ihist].GetBinContent(ibin)/2
                    else:
                        content += new_hist_list[ihist].GetBinContent(ibin)
                        error_stat += new_h_err_stat_list[ihist].GetBinContent(ibin)*new_h_err_stat_list[ihist].GetBinContent(ibin)
                        error_syst += new_h_err_syst_list[ihist].GetBinContent(ibin)
            if len(new_hist_list) == 1:
                new_cent_20_hist.SetBinContent(ibin, content)
                error_stat = np.sqrt(error_stat)
            else:
                error_stat = np.sqrt(error_stat)/2
                error_syst = error_syst/2
                new_cent_20_hist.SetBinContent(ibin, content/2)
            error = np.sqrt(error_stat*error_stat+error_syst*error_syst)
            new_cent_20_err_syst.SetBinContent(ibin, error_syst)
            new_cent_20_err_stat.SetBinContent(ibin, error_stat)
            new_cent_20_hist.SetBinError(ibin, error)

        #corr11.Rebin(2)
        new_cent_20_hist.Scale(0.5)
        hist_notrebin.Scale(0.5)

        return hist_notrebin, new_cent_20_hist, new_cent_20_err_stat, new_cent_20_err_syst