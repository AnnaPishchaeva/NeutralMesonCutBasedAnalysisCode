# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3
# This code was written by Alica Enderich (Febuary 2024)
# This code was modified by Julia Schlägel (July 2024)


import numpy as np
import datetime
import ROOT
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine
from ROOT import gStyle, gROOT, gSystem
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TH1
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle, kGray, kOpenCircle
import re
import numpy as np
import datetime
import math
import ctypes
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
gStyle.SetPadTickX(1);
gStyle.SetPadTickY(1);
gStyle.SetErrorX(0)
gStyle.SetEndErrorSize(5)
gStyle.SetLineScalePS(1)
from Histo_formatting import *

class Plot_inv_mass:
    """_summary_
    """

    def __init__(self):

        print("default constructor is called");
    def __init__(
        self,
        meson: str ,
        ssname, plotting = "none"
        ) -> str:
        """_summary_

        Args:
            meson (str): _description_
            filename (int): _description_
            dirname (_type_): _description_

        Returns:
            str: _description_
        """
        print("#################################################")
        print("Plot_inv_mass.py")
        print("#################################################")
        self.meson = meson;
        self.ssname = ssname
        self.arr_pt = np.array([0,1,2,3,4,5], dtype=float);
        self.number_columns = 1
        self.number_rows = 1
        if self.meson == "pi0":
            self.plotting_range = [0., 0.3];
        if self.meson == "eta":
            self.plotting_range = [0.3, 0.8]
        self.plotting = plotting

    def set_rows_columns(self, arr_pt):
        number = np.sqrt(len(arr_pt))
        right_number = math.ceil(number)
        self.number_rows = right_number
        self.number_columns = right_number

    def draw_line(self, drawing_option, pad, place, legend, yield_option, yMin, yMax, mean, npt, iPt):
        pad.cd(place)

        if drawing_option == "inv_mass_subtracted":
            line1 = set_line(self.fit_range_for_one_pT[0], yMin, self.fit_range_for_one_pT[0], 1.6*yMax, kGray+1,2,1)
            line1.Draw("");

            line2 = set_line(self.fit_range_for_one_pT[1], yMin, self.fit_range_for_one_pT[1], 1.6*yMax, kGray+1,2,1)
            line2.Draw("");
            if (iPt == npt-1 or npt == 1):
                legend.AddEntry(line1, "fitting range", "l")

        if drawing_option == "inv_mass_same_mix":
            line1 = set_line(self.scale_mix_range_for_one_pT[0], yMin, self.scale_mix_range_for_one_pT[0], 1.6*yMax, kGray+2,1,1)
            line1.Draw("");

            line2 = set_line(self.scale_mix_range_for_one_pT[1], yMin, self.scale_mix_range_for_one_pT[1], 1.6*yMax, kGray+2,1,1)
            line2.Draw("");
            if (iPt == npt-1 or npt == 1):
                legend.AddEntry(line1, "integration range for scaling", "l")

        if yield_option =="yield":
            yield_line1 = set_line(mean-self.yield_range_for_one_pT[0], yMin, mean-self.yield_range_for_one_pT[0], 1.6*yMax, kGray+2,1,1)
            yield_line1.Draw("");
            yield_line2 = set_line(mean+self.yield_range_for_one_pT[1], yMin, mean+self.yield_range_for_one_pT[1], 1.6*yMax, kGray+2,1,1)
            yield_line2.Draw("");
            if (iPt == npt-1 or npt == 1):
                legend.AddEntry(yield_line1, "integration range of raw yield", "l")

    def draw_mean(self, pad, place, legend, npt, iPt, histogram_list2, yMin, fit_min):
        pad.cd(place)
        mean = histogram_list2[iPt].GetMaximumX();
        name = ""
        if mean == fit_min:
            if self.meson == "eta":
                mean = 0.535
            if self.meson == "pi0":
                mean = 0.135
        mean_line = set_line(mean, yMin, mean, 0.6*histogram_list2[iPt].GetMaximum(), kGreen+2,1,1)
        mean_line.Draw("");
        if (iPt == npt-1 or npt == 1):
            if self.meson == "eta":
                legend.AddEntry(mean_line, "mean #eta mass","l")
            if self.meson == "pi0":
                legend.AddEntry(mean_line, "mean #pi^{0} mass","l")
        return mean

    def draw_title(self, drawing_option, canvas, title, info_decay, is_one_pad = False, iPt = 0):
        canvas.cd()
        if drawing_option == "inv_mass_subtracted" and not is_one_pad:
            if info_decay.get_meson() == "pi0":
                title.AddText("Invariant mass of #pi^{0} with asymmetric Gaussian")
            if info_decay.get_meson() == "eta":
                title.AddText("Invariant mass of #eta with asymmetric Gaussian")
        if drawing_option == "inv_mass_same_mix" and not is_one_pad:
            title.AddText("Same and (scaled) mixed event");
        if drawing_option == "parameters":
            if info_decay.get_meson() == "pi0":
                title.AddText("Parameters for asymmetric Gaussian fit on invariant mass of #pi^{0}")
            if info_decay.get_meson() == "eta":
                title.AddText("Parameters for asymmetric Gaussian fit on invariant mass of #eta");
        if drawing_option == "parameters_all_cent":
            if info_decay.get_meson() == "pi0":
                title.AddText("Parameters for asymmetric Gaussian fit on invariant mass of #pi^{0}")
            if info_decay.get_meson() == "eta":
                title.AddText("Parameters for asymmetric Gaussian fit on invariant mass of #eta");
        if is_one_pad:
                startPt = info_decay.arr_pt[iPt+info_decay.get_start_pt()];
                endPt = info_decay.arr_pt[iPt+info_decay.get_start_pt()+1];
                title.AddText("{:.2f} GeV/#it{{c}} < #it{{p}}_{{T}}  < {:.2f} GeV/#it{{c}}".format(startPt, endPt));
        title.Draw();
        ROOT.SetOwnership(title,False);

    def set_range_param(self, i_par):
        if self.meson == "pi0":
            yMin_array = [1e-8, 120, 0, 0.1*1e-3,0, 0]
            yMax_array = [5e-2, 150, 0, 10*1e-3, 0, 0]
        if self.meson == "eta":
            yMin_array = [0, 500, 6*1e-3, 10*1e-4, 4*1e-3, 0, 0, 0]
            yMax_array = [0, 600 ,22*1e-3, 200*1e-4, 20*1e-3, 0, 0, 0]

        return yMin_array[i_par], yMax_array[i_par]

    def draw_parameters(self, pad, histo_list, place, decay_channel, typ, i_par):
        if histo_list[i_par].GetEntries() == 0.0:
            return

        pad.cd(place);
        pad.cd(place).SetTopMargin(0.15);
        pad.cd(place).SetBottomMargin(0.15);
        pad.cd(place).SetRightMargin(0.15);
        pad.cd(place).SetLeftMargin(0.25);

        title_par = histo_list[i_par].GetTitle()
        SetTitle(title_par)

        yMin_, yMax_ = 10000, -10000
        Nbins = histo_list[i_par].GetNbinsX()+1
        for ibin in range(Nbins):
            if yMax_ < histo_list[i_par].GetBinContent(ibin) and histo_list[i_par].GetBinContent(ibin) != 0:
                yMax_ = histo_list[i_par].GetBinContent(ibin)
            if yMin_ > histo_list[i_par].GetBinContent(ibin) and histo_list[i_par].GetBinContent(ibin) > 0:
                yMin_ = histo_list[i_par].GetBinContent(ibin)
        if i_par == 1: #a.k.a. if it's a mean
            yMax_ = yMax_*1000.*1.03
            yMin_ = yMin_*1000.*0.98
        else:
            yMax_ = 1.7*yMax_
            yMin_ = 0.4*yMin_

        frame = pad.cd(place).DrawFrame(0., yMin_, 20., yMax_);
        FrameSettings(frame, "#it{p}_{T} (GeV/#it{c})", "d#it{{N}}{}/d#it{{M}}{}".format(decay_channel, decay_channel), 0., 20.)
        if i_par ==0:
            pad.cd(place).SetLogy();
        if i_par == 1:
            histo_list[i_par].Scale(1000.)
        if self.meson == "pi0":
            yTitle_list = ["#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{#pi^{0}}}{d#it{p}_{T}} counts #upoint (GeV/#it{c})^{-1}",
                            "MeV/#it{c}^{2}","#it{#lambda} (GeV/#it{c}^{2})", "#it{#sigma} (GeV/#it{c}^{2})",
                            "#it{fwhm/2.36}(GeV/#it{c}^{2})", "#it{#chi}^{2}/#it{DoF}" ]
        if self.meson == "eta":
            yTitle_list = ["#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{#eta}}{d#it{p}_{T}} counts #upoint (GeV/#it{c})^{-1}",
                            "MeV/#it{c}^{2}","#it{#lambda} (GeV/#it{c}^{2})", "#it{#sigma} (GeV/#it{c}^{2})",
                            "#it{fwhm/2.36}(GeV/#it{c}^{2})", "#it{#chi}^{2}/#it{DoF}" ]

        frame.GetYaxis().SetTitle(yTitle_list[i_par]);
        frame.GetYaxis().SetTitleOffset(1.2)
        histo_list[i_par].GetYaxis().SetNdivisions(505, kTRUE)
        if  i_par != 5:
            if typ == "mc":
                DrawHisto(histo_list[i_par], kBlue, "E1 ,same", 1.5, kOpenCircle)
            else:
                DrawHisto(histo_list[i_par], kBlue, "E1 ,same", 1.5)
        else:
            if typ == "mc":
                DrawHisto(histo_list[i_par], kBlue, "P ,same", 1.5, kOpenCircle)
            else:
                DrawHisto(histo_list[i_par], kBlue, "P same", 1.5)
        histo_list[i_par].SetDirectory(ROOT.gDirectory)

    def draw_parameters_all_cent(self, pad, histo_list, place, decay_channel, color, typ, i_par):

        pad.cd(place);
        pad.cd(place).SetTopMargin(0.15);
        pad.cd(place).SetBottomMargin(0.15);
        pad.cd(place).SetRightMargin(0.15);
        pad.cd(place).SetLeftMargin(0.25);

        title_par = histo_list[0][i_par].GetTitle()
        SetTitle(title_par)

        yMin_, yMin_1 = 10000, 10000
        yMax_, yMax_1 = -10000, -10000
        for i in range(len(histo_list)-1):
            histo_list[i][i_par].BufferEmpty(-1)
            histo_list[i][i_par].BufferEmpty(-1)
            Nbins = histo_list[i][i_par].GetNbinsX()+1
            for ibin in range(Nbins):
                if yMax_ < histo_list[i][i_par].GetBinContent(ibin) and histo_list[i][i_par].GetBinContent(ibin) != 0:
                    yMax_ = histo_list[i][i_par].GetBinContent(ibin)
                if yMin_ > histo_list[i][i_par].GetBinContent(ibin) and histo_list[i][i_par].GetBinContent(ibin) > 0:
                    yMin_ = histo_list[i][i_par].GetBinContent(ibin)
            Nbins1 = histo_list[i+1][i_par].GetNbinsX()+1
            for ibin in range(Nbins1):
                if yMax_1 < histo_list[i+1][i_par].GetBinContent(ibin) and histo_list[i+1][i_par].GetBinContent(ibin) != 0:
                    yMax_1 = histo_list[i+1][i_par].GetBinContent(ibin)
                if yMin_1 > histo_list[i+1][i_par].GetBinContent(ibin) and histo_list[i+1][i_par].GetBinContent(ibin) > 0:
                    yMin_1 = histo_list[i+1][i_par].GetBinContent(ibin)
        if yMax_1 > yMax_:
            yMax_ = yMax_1
        if yMin_1 < yMin_:
            yMin_ = yMin_1
        if i_par == 1:
            if yMin_ < 1:
                yMax_ = yMax_*1000
                yMin_ = yMin_*1000
            yMax_ = 1.03*yMax_
            yMin_ = 0.98*yMin_
        else:
            yMax_ = 1.7*yMax_
            yMin_ = 0.4*yMin_

        if i_par == 0:
                #yMin_, yMax_ = 1e-8, 5e-2
                pad.cd(place).SetLogy()
        frame = pad.cd(place).DrawFrame(0., yMin_, 20., yMax_);
        FrameSettings(frame, "#it{p}_{T} (GeV/#it{c})", "d#it{{N}}{}/d#it{{M}}{}".format(decay_channel, decay_channel), 0., 20.)
        if self.meson == "pi0":
            yTitle_list = ["#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{#pi^{0}}}{d#it{p}_{T}} (GeV/#it{c})^{-1}",
                            "MeV/#it{c}^{2}","#it{#lambda} (GeV/#it{c}^{2})", "#it{#sigma} (GeV/#it{c}^{2})",
                            "#it{fwhm/2.36}(GeV/#it{c}^{2})", "#it{#chi}^{2}/#it{DoF}" ]
        if self.meson == "eta":
            yTitle_list = ["#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{#eta}}{d#it{p}_{T}} (GeV/#it{c})^{-1}",
                            "MeV/#it{c}^{2}","#it{#lambda} (GeV/#it{c}^{2})", "#it{#sigma} (GeV/#it{c}^{2})",
                            "#it{fwhm/2.36}(GeV/#it{c}^{2})", "#it{#chi}^{2}/#it{DoF}" ]

        print(i_par, "i_par")
        frame.GetYaxis().SetTitle(yTitle_list[i_par])
        for icut in range(len(histo_list)):
            if i_par == 1: #and self.plotting == "mc_and_data":
                histo_list[icut][i_par].Scale(1000.)
            if  i_par != 5:
                if typ == "mc" or (icut == 1 and self.plotting == "mc_and_data"):
                    DrawHistoCombined(histo_list[icut][i_par], color[icut], "E1 ,same", 1.5, kOpenCircle)
                else:
                    DrawHistoCombined(histo_list[icut][i_par], color[icut], "E1 ,same", 1.5)
            else:
                if typ == "mc" or (icut == 1 and self.plotting == "mc_and_data"):
                    DrawHistoCombined(histo_list[icut][i_par], color[icut], "P ,same", 1.5, kOpenCircle)
                else:
                    DrawHistoCombined(histo_list[icut][i_par], color[icut], "P same", 1.5)
            histo_list[icut][i_par].SetDirectory(ROOT.gDirectory)

    def draw_inv_mass_and_fit(self, drawing_option, pad, info_decay, iPt, place, histogram_list, histogram_list2, is_one_pad):
        no_entries = histogram_list[iPt].GetEntries()
        print("entries histo: ", no_entries)
        if no_entries == 0: #if the number of entries is really low
            print ("centrality: ",info_decay.get_cent(),"pt: ",iPt, "entries:", no_entries)
            return 0,0
        startPt = info_decay.arr_pt[iPt+info_decay.get_start_pt()];
        print(info_decay.arr_pt[iPt], "+3","+4")
        endPt = info_decay.arr_pt[iPt+info_decay.get_start_pt()+1];

        if not is_one_pad:
            pad.SetFillColor(0)
            print("is one pad: ", is_one_pad)
            print("place ", place)
            pad.cd(place).SetTopMargin(0.15);
            pad.cd(place).SetBottomMargin(0.15);
            pad.cd(place).SetRightMargin(0.15);
            pad.cd(place).SetLeftMargin(0.15); # -> change values later to make neater

        title_pt = "{:.2f} GeV/#it{{c}} < #it{{p}}_{{T}}  < {:.2f} GeV/#it{{c}}".format(startPt, endPt);

        yMin = 0.;
        yMax = 0.;
        if isinstance(histogram_list[iPt],TH1):
            for i in range (histogram_list[iPt].GetXaxis().FindBin(self.plotting_range[0]), histogram_list[iPt].GetXaxis().FindBin(self.plotting_range[1])):
                if histogram_list[iPt].GetBinContent(i) < yMin:
                    if drawing_option == "inv_mass_subtracted":
                        yMin = -20
                    if drawing_option == "inv_mass_same_mix":
                        yMin = histogram_list[iPt].GetBinContent(i)

                if histogram_list[iPt].GetBinContent(i) > yMax:
                    yMax = histogram_list[iPt].GetBinContent(i);

        else:
            yMax = histogram_list[iPt].GetMaximum(self.plotting_range[0], self.plotting_range[1]);
            if drawing_option == "inv_mass_subtracted":
                yMin = -40
            if drawing_option == "inv_mass_same_mix":
                yMin = histogram_list[iPt].GetMinimum(self.plotting_range[0], self.plotting_range[1]);

        frame = pad.cd(place).DrawFrame(self.plotting_range[0], yMin, self.plotting_range[1], 1.6*yMax);

        FrameSettings(frame, "#it{{M}}{} (GeV/#it{{c}}^{{2}})".format(info_decay.get_decay_channel()), "d#it{{N}}{}/d#it{{M}}{} [counts]".format(info_decay.get_decay_channel(), info_decay.get_decay_channel()),
                    self.plotting_range[0], self.plotting_range[1], is_one_pad)

        if info_decay.get_typ() == "mc":
            DrawHisto(histogram_list[iPt], kBlue, "Hsame", 1.5, kOpenCircle)
        else:
            DrawHisto(histogram_list[iPt], kBlue, "Hsame", 1.5);
        DrawHisto(histogram_list2[iPt], kRed, "same", 1.5);

        if is_one_pad:
            frame.GetYaxis().SetLabelSize(0.045);
            frame.GetYaxis().SetTitleSize(0.045);
            frame.GetXaxis().SetLabelSize(0.045);
            frame.GetXaxis().SetTitleSize(0.045);
            frame.GetYaxis().SetTitleOffset(1.4);
        else:
            SetTitle(title_pt)

        return yMin, yMax

    def draw_legend(self, drawing_option, pad_legend, legend, info_decay, histogram_list, histogram_list2, name_plot, BckNmb, is_one_pad = False):
        if not is_one_pad:
            pad_legend.cd()
        example_bin = 0
        npt = len(histogram_list)
        if drawing_option == "inv_mass_same_mix":
            legend.AddEntry(histogram_list[example_bin],"same evt. #it{{M}}{}".format(info_decay.decay_channel),"ep");
            legend.AddEntry(histogram_list2[example_bin], "mixed evt. #it{{M}}{} (scaled)".format(info_decay.decay_channel), "ep");
        if drawing_option == "inv_mass_subtracted":
            legend.AddEntry(histogram_list[example_bin],"same evt. #it{{M}}{} (Signal-BG)".format(info_decay.decay_channel),"ep");
            linesize         = histogram_list2[example_bin].GetLineWidth();
            histogram_list2[example_bin].SetLineWidth(linesize);
            if BckNmb == 0:
                if name_plot.Contains("FixedPzPiZero") == True:
                    legend.AddEntry(histogram_list2[example_bin], "{} #it{{M}}{} (p_{{z}} of #pi^{{0}} fixed)".format(info_decay.background_fit, info_decay.decay_channel), "l")
                else:
                    legend.AddEntry(histogram_list2[example_bin], "{} #it{{M}}{}".format(info_decay.background_fit, info_decay.decay_channel), "l")
            else:
                if name_plot.Contains("FixedPzPiZero") == True:
                    legend.AddEntry(histogram_list2[example_bin], "{} #it{{M}}{} group {} (p_{{z}} of #pi^{{0}} fixed)".format(info_decay.background_fit, info_decay.decay_channel, BckNmb), "l")
                else:
                    legend.AddEntry(histogram_list2[example_bin], "{} #it{{M}}{} group {}".format(info_decay.background_fit, info_decay.decay_channel, BckNmb), "l")

    def draw_legend_param(self, pad_legend, legend, line, start_pt, histo_list):
        pad_legend.cd()
        example_bin = self.number_columns+start_pt-1;
        marker_size = histo_list[example_bin].GetMarkerSize();
        histo_list[example_bin].SetMarkerSize(marker_size);
        legend.AddEntry(histo_list[example_bin],"fitparameter for #it{p}_{T} ranges","ep");
        linesize = histo_list[example_bin].GetLineWidth();
        histo_list[example_bin].SetLineWidth(linesize);
        if self.meson == "pi0":
            legend.AddEntry(line, "literature value of #it{#pi}^{0} mass", "L")
        if self.meson == "eta":
            legend.AddEntry(line, "literature value of #it{#eta} mass", "L")
        legend.Draw();

    def draw_legend_param_all_cent(self, pad_legend, legend, line, start_pt, histo_list, ssname_list, color, info_decay):
        pad_legend.cd()
        example_bin = self.number_columns+start_pt-1;
        marker_size = histo_list[0][example_bin].GetMarkerSize();
        for icut in range(len(histo_list)):
            histo_list[icut][example_bin].SetMarkerSize(marker_size);
            if self.plotting == "centralities":
                info_decay.set_centrality(ssname_list[icut])
                cent1, cent2 = info_decay.get_cent()
                legend.AddEntry(histo_list[icut][example_bin],"{0}% - {1}%".format(cent1, cent2),"ep");
            if self.plotting == "RV0s":
                info_decay.set_RV0(ssname_list[icut])
                V0_radius_min, V0_radius_max = info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max()
                legend.AddEntry(histo_list[icut][example_bin],"{0}<#it{{R}}_{{V^{{0}}}}<{1} cm".format(int(V0_radius_min), int(V0_radius_max)),"ep");
            if self.plotting == "mc_and_data":
                if icut == 0:
                    legend.AddEntry(histo_list[icut][example_bin],"Data","ep")
                if icut == 1:
                    legend.AddEntry(histo_list[icut][example_bin],"MC","ep")
        legend.Draw();


    def set_all_labels(self, txt_alice, start_txt_x, start_txt_y, txt_height, diff_txt, info_decay, is_one_pad = False):
        labels = []

        if is_one_pad:
            txt_height = txt_height*0.7
            y_diff = 3
            y_start = 7
        else:
            y_diff = 4
            y_start = 11
        alice = set_labels(txt_alice, start_txt_x, start_txt_y, txt_height, diff_txt, y_start+y_diff, 2, -2, 0.2)
        labels.append(alice)

        occup1, occup2 = info_decay.get_occupancy()
        if occup1 == "0" and occup2 == "10k":
            occup_range_latex = "0 < FT0Occ < 10^{4}"
        elif occup1 == "10k" and occup2 == "30k":
            occup_range_latex = "10^{4} < FT0Occ < 3x10^{4}"
        elif occup1 == "30k":
            occup_range_latex = "FT0Occ > 3x10^{4}"
        else:
            occup_range_latex = "FT0Occ: all"
        occup = set_labels(occup_range_latex, start_txt_x, start_txt_y, txt_height, diff_txt, y_start, 2, -2, 0.2)
        labels.append(occup)

        system = info_decay.get_system()
        if '-' in system:
            list_system = system.split("-")
            system = list_system[0] + "#font[122]{-}" + list_system[1]
        if info_decay.get_system() == "pp":
            energy_symbol = "#sqrt{s}"
        else:
            energy_symbol = "#sqrt{#it{s}_{NN}}"
        energy = energy_symbol + " = " + info_decay.get_energy()
        system_energy = "{0}, {1}".format(system, energy)
        latex_system_energy = set_labels(system_energy, start_txt_x, start_txt_y, txt_height, diff_txt, y_start-y_diff, 2, -3, 0.2)
        labels.append(latex_system_energy)

        period_typ = "{0}, {1}".format(info_decay.get_period(), info_decay.get_typ())
        latex_period_typ = set_labels(period_typ, start_txt_x, start_txt_y, txt_height, diff_txt, y_start-y_diff*2, 2, -4, 0.2)
        labels.append(latex_period_typ)

        reconstr = set_labels("#gamma rec. with PCM", start_txt_x, start_txt_y, txt_height, diff_txt, y_start-y_diff*3+1, 2, -6, 0.2)
        labels.append(reconstr)

        decay_name = set_decay_name(info_decay.get_meson(), info_decay.get_decay())
        latex_decay_name = set_labels(decay_name, start_txt_x, start_txt_y, txt_height, diff_txt, y_start-y_diff*4+0.5, 2, -6, 0.2)
        labels.append(latex_decay_name)

        if info_decay.get_system() != "pp":
            cent1, cent2 = info_decay.get_cent()
            latex_cent = set_labels("{0}% - {1}% centrality".format(cent1, cent2), start_txt_x, start_txt_y, txt_height, diff_txt, y_start-5*y_diff+0.5, 2, -6, 0.2)
            labels.append(latex_cent)
        else:
            y_start+=y_diff-0.5

        RV0_min = info_decay.get_V0_radius_min()
        RV0_max = info_decay.get_V0_radius_max()
        latex_RV0 = set_labels("{0} < #it{{R}}_{{V^{{0}}}} < {1} cm".format(int(RV0_min), int(RV0_max)), start_txt_x, start_txt_y, txt_height, diff_txt, y_start-6*y_diff, 2, -6, 0.2)
        labels.append(latex_RV0)

        nev = info_decay.get_Nev_cent()
        Nev = "Number of events: {}".format(f"{nev:.2e}")

        nevents = set_labels(Nev, start_txt_x, start_txt_y, txt_height, diff_txt, y_start-7*y_diff, 2, -6, 0.2)
        labels.append(nevents)

        return labels

    def plot_inv_mass(self, drawing_option, info_decay, histogram_list, histogram_list2, name_plot, yield_option, is_one_pad = False):
        TGaxis.SetMaxDigits(3);
        self.set_rows_columns(info_decay.arr_pt)
        if is_one_pad:
            self.number_rows, self.number_columns = 1, 1

        if is_one_pad:
            npt = 1
        else:
            #npt = len(info_decay.arr_pt)-info_decay.start_pt-1
            npt = len(histogram_list)
        print("npt", npt)
        print("####start pt_inv: ", info_decay.start_pt)     #start pt is here if the first few data points shall not be considered

        BckNmb=0
        margin_width_leg = 0.05;
        if not is_one_pad:
            start_txt_x     = 0.10;
            start_txt_y     = 0.5 #0.7 for larger binning, 0.5 for smaller binning
        n_pix        = 13;
        if self.number_columns > 7:
            margin_width_leg = 0.25;
            start_txt_x          = 0.05
            n_pix             = 12;

        if is_one_pad:
            canvas = TCanvas("c0","c0",0,0,800,800)
            pad = canvas.cd()
            pad.SetPad(0,0.5,3,3);
            pad.SetMargin(0.15,0.05,0.1,0.15)
            pad.SetTicks(1,1)
        else:
            canvas = set_canvas(5000, 3000, 1, 1)
            canvas.cd()
            pad = set_pad(-0.0, 0.0, 1., 0.85, 0, self.number_columns, self.number_rows, False)
            pad.Draw()

        if is_one_pad:
            title_size = 0.045
            title = set_txt(0.25,0.88,0.9,0.88,"NDC",title_size, 12)
        else:
            title_size = 0.035
            title = set_txt(0.,0.85,0.9,0.95,"NDC",title_size, 12)

        pad_legend = set_pad(-0.0, 0.0, 1., 1., 0, self.number_columns, self.number_rows, True)
        pad_legend.Draw()
        txt_alice, txt_height, diff_txt = set_txt_pad_legend(pad_legend, self.number_columns, self.number_rows, info_decay, n_pix)

        if is_one_pad:
            start_txt_x = 0.18
            txt_height = 0.009
            start_txt_y = 0.67
            diff_txt = txt_height*1.50
            legend = set_legend(0.66,0.67,1.16,0.82, 0.02, 1, TString(""), 42, margin_width_leg)
        else:
            legend = set_legend(start_txt_x, start_txt_y-8*diff_txt, 0.8, start_txt_y-20*diff_txt-0.02, 3*n_pix, 1, TString(""), 43, margin_width_leg)

        place = 0
        for iPt in range (0, npt):
            place+=1
            if place == self.number_columns:
                place +=1
            if is_one_pad:
                place = 0
                iPt = 3
                print("iPt: ", iPt, "len: ", len(histogram_list))
                if iPt >= len(histogram_list):
                    iPt = len(histogram_list) - 1
                    print("iPt: ", iPt, "len: ", len(histogram_list))
            print("iPt", iPt, "N", npt, flush=True)
            yMin, yMax = self.draw_inv_mass_and_fit(drawing_option, pad, info_decay, iPt, place, histogram_list, histogram_list2, is_one_pad)
            if (yMin == 0 and yMax ==0):
                continue
            mean = 0
            fit_range = info_decay.get_fit_range()
            if len(fit_range) == 1:
                self.fit_range_for_one_pT = fit_range[0]
            else:
                self.fit_range_for_one_pT = fit_range[iPt]

            if drawing_option == "inv_mass_subtracted":
                mean = self.draw_mean(pad, place, legend, npt, iPt, histogram_list2, yMin, self.fit_range_for_one_pT[0])

            scale_mix_range = info_decay.get_scale_mix_range()
            if len(scale_mix_range) == 1:
                self.scale_mix_range_for_one_pT = scale_mix_range[0]
            else:
                self.scale_mix_range_for_one_pT = scale_mix_range[iPt]

            yield_range = info_decay.get_yield_range()
            if len(yield_range) == 1:
                self.yield_range_for_one_pT = yield_range[0]
            else:
                self.yield_range_for_one_pT = yield_range[iPt]
            self.draw_line(drawing_option, pad, place, legend, yield_option, yMin, yMax, mean, npt, iPt)

        if is_one_pad:
            self.draw_title(drawing_option, canvas, title, info_decay, is_one_pad, iPt)
        else:
            self.draw_title(drawing_option, canvas, title, info_decay)
        self.draw_legend(drawing_option, pad_legend, legend, info_decay, histogram_list, histogram_list2, name_plot, BckNmb, is_one_pad)

        legend.Draw()
        labels = self.set_all_labels(txt_alice, start_txt_x, start_txt_y, txt_height, diff_txt, info_decay, is_one_pad)
        for label in labels:
            label.Draw()

        canvas.SaveAs(name_plot.Data());
        ROOT.SetOwnership(canvas, False)
        del pad_legend;
        del pad;
        del canvas;
        del legend


    def plot_parameters(self,info_decay, histogram_list, name_plot, ssname_list = [], is_all_cent = False):

        self.number_rows, self.number_columns = 3,3

        margin_width_leg = 0.15;
        start_txt_x     = 0.25;
        n_pix        = 13;
        if self.number_columns > 7:
            margin_width_leg = 0.25;
            start_txt_x          = 0.05
            n_pix             = 12;

        canvas = set_canvas(3000, 2000, 1, 1)
        canvas.cd()

        pad = set_pad(-0.0, 0.0, 1., 0.85, 0, self.number_columns, self.number_rows)
        pad.Draw()

        title_size = 0.040
        title = set_txt(0.,0.88,1.,0.98,"NDC", title_size, 22)

        pad_legend = set_pad(-0.0, 0.0, 1., 1., 0, self.number_columns, self.number_rows, True)
        pad_legend.Draw()

        txt_alice, txt_height, diff_txt = set_txt_pad_legend(pad_legend, self.number_columns, self.number_rows, info_decay, n_pix)
        if self.plotting == "RV0s":
            txt_size_legend = 35
            start_txt_y     = 0.55
            ncolums         = 1
        else:
            txt_size_legend = 2.8*n_pix
            ncolums         = 1
            if is_all_cent:
                ncolums     = 2
            start_txt_y     = 0.5
        legend = set_legend(start_txt_x, start_txt_y-8*diff_txt, 0.8, start_txt_y-20*diff_txt-0.02, txt_size_legend, ncolums, TString(""), 43, margin_width_leg)

        place = 0
        length = 6
        if is_all_cent:
            if self.plotting == "mc_and_data":
                color = [ROOT.kBlue, ROOT.kBlue]
            else:
                color = [kRed, 802, 401, kGreen+1, kGreen-1, kCyan ,kBlue, kMagenta, kMagenta+3, 922, 923, kBlack]

        for i_par in range (0, length):
            place+=1
            if place == self.number_columns:
                place +=1

            if is_all_cent:
                self.draw_parameters_all_cent(pad, histogram_list, place, info_decay.get_decay_channel(), color, info_decay.get_typ(), i_par)
            else:
                self.draw_parameters(pad, histogram_list, place, info_decay.get_decay_channel(), info_decay.get_typ(), i_par)

            pad.cd(place)
            if i_par == 1:
                if info_decay.get_meson() == "pi0":
                    line = set_line(0,134.97,20,134.97, kRed, 2, 1)
                if info_decay.get_meson() == "eta":
                    line = set_line(0,547.86,20,547.86, kRed, 2, 1)
                line.Draw("");
                ROOT.SetOwnership(line,False);

        if is_all_cent:
            self.draw_title("parameters_all_cent", canvas, title, info_decay)
            self.draw_legend_param_all_cent(pad_legend, legend, line, info_decay.get_start_pt(), histogram_list, ssname_list, color, info_decay)
        else:
            self.draw_title("parameters", canvas, title, info_decay)
            self.draw_legend_param(pad_legend, legend, line, info_decay.get_start_pt(), histogram_list)
        labels = self.set_all_labels(txt_alice, start_txt_x, start_txt_y, txt_height, diff_txt, info_decay)
        legend.Draw()
        for label in labels:
            label.Draw()

        canvas.SaveAs(name_plot.Data());
        ROOT.SetOwnership(canvas, False)

        del pad_legend;
        del pad;
        del canvas;
        del legend


