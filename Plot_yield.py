# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3
# This code was written by Alica Enderich (Febuary 2024)
# This code was modified by Julia Schlägel (July 2024)


import numpy as np
import datetime
import ROOT
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine, TColor
from ROOT import gStyle, gROOT, gSystem, gPad
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TH1
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle, kOpenCircle, kRainBow
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

class Plot_yield:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, decay, meson):
        print("#################################################")
        print("Plot_yield.py")
        print("#################################################")
        self.meson = meson
        self.decay = decay

    def set_N_coll(self):
        if self.cent1 == 0 and self.cent2 == 20:
            N_coll = 1286.1
        if self.cent1 == 20 and self.cent2 == 40:
            N_coll = 494.0
        if self.cent1 == 40 and self.cent2 == 60:
            N_coll = 155.0
        if self.cent1 == 60 and self.cent2 == 80:
            N_coll = 32.1
        if self.cent1 == 80 and self.cent2 == 100:
            N_coll = 4.4
        if self.cent1 == 0 and self.cent2 ==100:
            N_coll = 394.3
        if self.cent1 == 0 and self.cent2 == 10:
            N_coll = 1584.0
        if self.cent1 == 10 and self.cent2 == 20:
            N_coll = 988.2
        if self.cent1 == 20 and self.cent2 == 30:
            N_coll = 618.0
        if self.cent1 == 30 and self.cent2 == 40:
            N_coll = 369.9
        if self.cent1 == 40 and self.cent2 == 50:
            N_coll = 206.1
        if self.cent1 == 50 and self.cent2 == 60:
            N_coll = 103.7
        if self.cent1 == 60 and self.cent2 == 70:
            N_coll = 46.1
        if self.cent1 == 70 and self.cent2 == 80:
            N_coll = 18.1
        if self.cent1 == 80 and self.cent2 == 90:
            N_coll = 6.5
        if self.cent1 == 90 and self.cent2 == 100:
            N_coll = 2.3
        return N_coll

    def set_y_min_y_max(self, fHistoParameter, type_of_histo, scaling, ratio):
        yMax_, yMax_1 = -100, -100
        yMin_, yMin_1 = 100, 100
        if len(fHistoParameter) == 1:
            Nbins = fHistoParameter[0].GetNbinsX()+1
            for ibin in range(Nbins):
                if yMax_ < fHistoParameter[0].GetBinContent(ibin) and fHistoParameter[0].GetBinContent(ibin) > 0:
                    yMax_ = fHistoParameter[0].GetBinContent(ibin)
                if yMin_ > fHistoParameter[0].GetBinContent(ibin) and fHistoParameter[0].GetBinContent(ibin) > 0:
                    yMin_ = fHistoParameter[0].GetBinContent(ibin)
        else:
            for i in range(len(fHistoParameter)-1):
                fHistoParameter[i].BufferEmpty(-1)
                fHistoParameter[i+1].BufferEmpty(-1)
                Nbins = fHistoParameter[i].GetNbinsX()+1
                for ibin in range(Nbins):
                    if yMax_ < fHistoParameter[i].GetBinContent(ibin) and fHistoParameter[i].GetBinContent(ibin) != 0:
                        yMax_ = fHistoParameter[i].GetBinContent(ibin)
                    if yMin_ > fHistoParameter[i].GetBinContent(ibin) and fHistoParameter[i].GetBinContent(ibin) > 0:
                        yMin_ = fHistoParameter[i].GetBinContent(ibin)
                Nbins1 = fHistoParameter[i+1].GetNbinsX()+1
                for ibin in range(Nbins1):
                    if yMax_1 < fHistoParameter[i+1].GetBinContent(ibin) and fHistoParameter[i+1].GetBinContent(ibin) != 0:
                        yMax_1 = fHistoParameter[i+1].GetBinContent(ibin)
                    if yMin_1 > fHistoParameter[i+1].GetBinContent(ibin) and fHistoParameter[i+1].GetBinContent(ibin) > 0:
                        yMin_1 = fHistoParameter[i+1].GetBinContent(ibin)
                if yMax_1 > yMax_:
                    yMax_ = yMax_1
                if yMin_1 < yMin_:
                    yMin_ = yMin_1

        if ratio:
            yMax_ = 1.01*yMax_
            yMin_ = 0.01*yMin_
        else:
            yMax_ = 1.4*yMax_
            yMin_ = 0.4*yMin_
            if len(fHistoParameter)!=1:
                yMax_ = 1.2*yMax_*math.pow(scaling, len(fHistoParameter))
                yMin_ = yMin_
                if type_of_histo == "RawYield":
                    yMin_ = 0.001*yMin_
        return yMin_, yMax_

    def set_y_title(self, type_of_histo, scale_pT):
        if (type_of_histo == "RawYield" and scale_pT):
            y_title = "#frac{1}{#it{N}_{ev}} #frac{1}{#it{p}_{T}}#frac{d#it{N}^{raw}}{d#it{p}_{T}} (GeV/#it{c})^{-2}"
        elif (type_of_histo == "RawYield" and scale_pT == False):
            y_title = "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{raw}}{d#it{p}_{T}} (GeV/#it{c})^{-1}"
        elif (type_of_histo == "CorrectedYield" and scale_pT):
            if self.meson == "pi0":
                y_title = "#frac{1}{2#pi} #frac{1}{#it{p}_{T}} #frac{d^{2}#it{N}^{#pi^{0}}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}"
            if self.meson == "eta":
                y_title = "#frac{1}{2#pi} #frac{1}{#it{p}_{T}} #frac{d^{2}#it{N}^{#eta}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}"
        elif (type_of_histo == "CorrectedYield" and scale_pT == False):
            if self.meson == "pi0":
                y_title = "#frac{1}{2#pi} #frac{d^{2}#it{N}^{#pi^{0}}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}"
            if self.meson == "eta":
                y_title = "#frac{1}{2#pi} #frac{d^{2}#it{N}^{#eta}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}"
        if (type_of_histo == "CorrectedYieldRcp" or type_of_histo == "RawYieldRcp"):
            y_title = "R_{CP}"
        if (type_of_histo == "Significance"):
            y_title = "#it{S}/#sqrt{#it{S}+#it{B}}"
        if (type_of_histo == "StatUncert"):
            y_title = "stat. uncert. %"
        return y_title

    def set_txt_on_canvas(self, txt, plotting, is_mc_and_data, type_of_histo, title, length, info_decay):
        txt.AddText("{0}".format(TString(title)))
        if length == 1:
            txt.SetTextSize(0.040)
            if type_of_histo == "RawYield":
                set_txt_of_plot(info_decay, 0.57, 0.78, 0.87, 0.88, False, 0.025)
            elif type_of_histo == "CorrectedYield":
                set_txt_of_plot(info_decay, 0.57, 0.78, 0.87, 0.88, True, 0.025)
            elif type_of_histo == "Significance":
                set_txt_of_plot(info_decay, 0.57, 0.78, 0.87, 0.88, False, 0.025)
            elif type_of_histo == "StatUncert":
                set_txt_of_plot(info_decay, 0.48, 0.48, 0.88, 0.58, False, 0.023)
        else:
            if type_of_histo == "RawYield":
                if plotting == "centralities":
                    set_txt_of_plot(info_decay, 0.58, 0.78, 0.88, 0.88, False, 0.025)
                else:
                    set_txt_of_plot(info_decay, 0.65, 0.82, 0.95, 0.92, False, 0.025)
            elif type_of_histo == "CorrectedYield":
                if plotting == "centralities":
                    set_txt_of_plot(info_decay, 0.58, 0.78, 0.88, 0.88, True, 0.025)
                else:
                    set_txt_of_plot(info_decay, 0.65, 0.82, 0.95, 0.92, True, 0.025)
            elif type_of_histo == "Significance":
                set_txt_of_plot(info_decay, 0.65, 0.82, 0.95, 0.92, False, 0.025)
            elif type_of_histo == "StatUncert":
                if is_mc_and_data:
                    set_txt_of_plot(info_decay, 0.59, 0.78, 0.89, 0.88, True , 0.025, True)
                else:
                    set_txt_of_plot(info_decay, 0.65, 0.82, 0.95, 0.92, True , 0.025)
        txt.Draw();
        ROOT.SetOwnership(txt,False);

    def set_legend_on_canvas(self, length, type_of_histo, plotting, is_mc_and_data):
        if length == 1:
            if type_of_histo == "StatUncert":
                leg = TLegend(0.69,0.27,0.81,0.31)
            else:
                leg = TLegend(0.65,0.57,0.77,0.61)
        else:
            if type_of_histo == "RawYield":
                if plotting == "centralities":
                    leg = TLegend(0.18,0.12,0.43,0.32)
                if plotting == "RV0s" or plotting == "diff_cent_and_diff_RV0s":
                    leg = TLegend(0.2,0.32,0.32,0.46)
                if plotting == "occupancies" or plotting == "diff_cent_and_diff_occupancies":
                    leg = TLegend(0.2,0.33,0.3,0.47)
            elif type_of_histo == "CorrectedYield":
                if plotting == "centralities":
                    leg = TLegend(0.66,0.41,0.88,0.65)
                if plotting == "RV0s" or plotting == "diff_cent_and_diff_RV0s":
                    leg = TLegend(0.71,0.51,0.83,0.65)
                if plotting == "occupancies" or plotting == "diff_cent_and_diff_occupancies":
                    leg = TLegend(0.7,0.56,0.8,0.70)
            elif type_of_histo == "StatUncert":
                if is_mc_and_data:
                    leg = TLegend(0.66,0.49,0.78,0.63)
                else:
                    leg = TLegend(0.53,0.51,0.65,0.65)
            else:
                leg = TLegend(0.53,0.51,0.65,0.65)
        leg.SetTextSize(0.025);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextFont(42);#helvetica
        return leg

    def set_title_and_scale(self, length, type_of_histo, scale_pT, plotting):
        if(type_of_histo == "RawYield"):
            scaling = 10
            if scale_pT:
                title = "Invariant Raw Yield"
            else:
                title = "Raw Yield"
        if(type_of_histo == "CorrectedYield"):
            scaling = 10
            if scale_pT:
                title = "Invariant Corrected Yield"
            else:
                title = "Corrected Yield"
        if(type_of_histo == "Significance"):
            title = "Significance"
        if(type_of_histo == "StatUncert"):
            title = "Statistical Uncertainty"
        if length == 1 or plotting != "centralities" or (type_of_histo != "RawYield" and type_of_histo != "CorrectedYield"):
            scaling = 1
        return title, scaling

    def set_fit_fun(self, A_e, M_pi, T_e, A_init, T_init, n):
        #defining the fit function
        f1TCM = TF1("TCM_fit", "[0]* (TMath::Exp(-(TMath::Power(TMath::Power(x,2)+TMath::Power([1],2), 0.5)-[1])/([2])))+ [3] * (TMath::Power((1+(TMath::Power(x,2))/(TMath::Power([4],2)*[5])), -[5]))", 0.2, 15.,6)
        f1TCM.SetParameters(A_e, M_pi, T_e, A_init, T_init, n)
        f1TCM.SetLineColor(ROOT.kBlack)

        f1TCM.SetParLimits(0, 0.5*A_e, 1.5*A_e)
        f1TCM.SetParLimits(1, M_pi, M_pi)
        f1TCM.SetParLimits(2, 0.5 * T_e, 1.5*T_e);
        f1TCM.SetParLimits(3,0.5*A_init, 1.5*A_init)
        f1TCM.SetParLimits(4, 0.5 * T_init, 1.5 * T_init)
        f1TCM.SetParLimits(5, 0.5 * n,  1.5*n)

        return f1TCM

    def plot_histo(self, type_of_histo, fHistoParameter, namePlot, info_decay, plotting, ratio = False, scale_pT = False, is_mc_and_data = False, period_mc = 0):
        TGaxis.SetMaxDigits(3);
        ssnames = []
        titlePt = fHistoParameter[0].GetTitle(); # ACCESS TITLE
        title, scaling = self.set_title_and_scale(len(fHistoParameter), type_of_histo, scale_pT, plotting)

        #finding minimum and max of the histogram
        yMin_, yMax_ = self.set_y_min_y_max(fHistoParameter, type_of_histo, scaling, False)

        c1 = TCanvas("c0","c0",0,0,800,800);
        if ratio:
            c1.Divide(1,2,1e-3,1e-3);
            p1 = c1.cd(1);
            p1.SetPad(0,0.3, 1,1)
            p1.SetMargin(0.15,0.02,0.,0.15)
            p1.SetBottomMargin(0.005)
        else:
            p1 = c1.cd()
            p1.SetPad(0,0.1,2,2);
            p1.SetMargin(0.18,0.1,0.1,0.15);
            #p1.SetTicks(1,1);
        if type_of_histo == "RawYield" or type_of_histo == "CorrectedYield" or type_of_histo == "StatUncert":
            p1.SetLogy();

        y_title = self.set_y_title(type_of_histo, scale_pT)
        frame1 = p1.DrawFrame(0., yMin_, 20., yMax_);
        set_frame_basics(frame1, "#it{p}_{T} (GeV/#it{c})", y_title, 0., 20.)
        if len(fHistoParameter)==1:
            if type_of_histo == "RawYield":
                offset_y_title = 2.5
            else:
                offset_y_title = 2.1
            offset_x_title = 1.4
            label_y_size = 0.03
            title_y_size = 0.03
            title_x_size = 0.03
        else:
            offset_y_title = 1.51
            offset_x_title = 1.0
            label_y_size = 0.04
            title_y_size = 0.045
            title_x_size = 0.045
        set_frame(frame1, title_x_size, title_y_size, offset_x_title, offset_y_title, 0.03, label_y_size, 0.01, 0.01)
        frame1.GetYaxis().SetMaxDigits(3); #1.53
        frame1.GetYaxis().SetNdivisions(510)
        ROOT.SetOwnership(frame1,False)

        leg = self.set_legend_on_canvas(len(fHistoParameter), type_of_histo, plotting, is_mc_and_data)
        markersize = 1.2

        histo = [0]*len(fHistoParameter)
        if len(fHistoParameter) == 1:
            ssnames.append(info_decay.get_ssname())
            color = [kBlue]
        else:
            color = [kRed, 802, 401, kGreen+1, kGreen-1, kCyan ,kBlue, kMagenta, kMagenta+3, 922, 923, kBlack]
            if (plotting == "centralities" or plotting == "RV0s"):
                ssnames = info_decay.get_ssname_list()
            if plotting == "diff_cent_and_diff_RV0s":
                ssnames.append(info_decay.get_ssname())
                RV0_list = info_decay.get_RV0_list()
            if plotting == "diff_cent_and_diff_occupancies":
                ssnames.append(info_decay.get_ssname())
                occup_list = info_decay.get_occupancy_list()
            if is_mc_and_data:
                ssnames.append(info_decay.get_ssname())
                color = [kBlue, kBlue]

        for icut in range(len(fHistoParameter)):
            histo[icut] = fHistoParameter[icut].Clone("{0}_{1}".format(type_of_histo, icut))
            if plotting == "centralities":
                histo[icut].Scale(math.pow(scaling, len(fHistoParameter)-icut))
            histo[icut].SetMarkerSize(markersize)
            if info_decay.get_typ() == "mc" and type_of_histo != "CorrectedYield" or (is_mc_and_data and icut == 1):
                histo[icut].SetMarkerStyle(kOpenCircle)
                histo[icut].SetLineWidth(6)
            else:
                histo[icut].SetMarkerStyle(kFullCircle)
            histo[icut].SetLineColor(color[icut])
            histo[icut].SetMarkerColor(color[icut])
            histo[icut].DrawCopy("SAME EP1")
            if plotting != "diff_cent_and_diff_occupancies" and plotting != "diff_cent_and_diff_RV0s" and not is_mc_and_data:
                info_decay.set_centrality(ssnames[icut])
                cent1, cent2 = info_decay.get_cent()
                info_decay.set_RV0(ssnames[icut])
                V0_radius_min, V0_radius_max = info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max()
            else:
                info_decay.set_centrality(ssnames[0])
                cent1, cent2 = info_decay.get_cent()
                info_decay.set_RV0(ssnames[0])
                V0_radius_min, V0_radius_max = info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max()

            if plotting == "centralities" and not is_mc_and_data and info_decay.get_system()!="pp":
                if len(fHistoParameter) == 1:
                    leg.AddEntry(histo[icut],"{0}%-{1}%".format(cent1, cent2),"ep")
                else:
                    leg.AddEntry(histo[icut],"{0}%-{1}% x {2}^{{{3}}}".format(cent1, cent2, scaling, len(fHistoParameter)-icut),"ep")
            if (plotting == "RV0s" or plotting == "diff_cent_and_diff_RV0s")and not is_mc_and_data:
                if plotting == "diff_cent_and_diff_RV0s":
                    if icut < len(RV0_list)-1:
                        V0_radius_min = RV0_list[icut]
                        V0_radius_max = RV0_list[icut+1]
                    else:
                        V0_radius_min = RV0_list[0]
                        V0_radius_max = RV0_list[len(RV0_list)-1]
                leg.AddEntry(histo[icut],"{0} < #it{{R}}_{{V^{{0}}}} < {1} cm".format(int(V0_radius_min), int(V0_radius_max)),"EP");
            if is_mc_and_data:
                if icut == 0:
                    leg.AddEntry(histo[icut],"Data {0}".format(info_decay.get_period()),"EP");
                else:
                    leg.AddEntry(histo[icut],"MC {0}".format(period_mc),"EP");
            if (plotting == "occupancies" or plotting == "diff_cent_and_diff_occupancies") and not is_mc_and_data:
                if plotting == "diff_cent_and_diff_occupancies":
                    if icut < len(occup_list)-1:
                        occup1 = occup_list[icut]
                        occup2 = occup_list[icut+1]
                    else:
                        occup1 = occup_list[0]
                        occup2 = occup_list[len(occup_list)-1]
                    if occup1 == "0" and occup2 == "10k":
                        occup_range_latex = "0 < FT0Occ < 10^{4}"
                    elif occup1 == "10k" and occup2 == "30k":
                        occup_range_latex = "10^{4} < FT0Occ < 3x10^{4}"
                    elif occup1 == "30k":
                        occup_range_latex = "FT0Occ > 3x10^{4}"
                    else:
                        occup_range_latex = "FT0Occ: all"
                leg.AddEntry(histo[icut],"{0}".format(occup_range_latex),"ep")

            histo[icut].SetDirectory(ROOT.gDirectory)

        if ratio:
            p2 = c1.cd(2);
            p2.SetPad(0,0,1,0.3);
            p2.SetMargin(0.15,0.02,0.22,0.0);
            p2.SetTicks(1,1);
            p2.SetTopMargin(0.05)
            p2.SetBottomMargin(0.3)
            h1ratio1 = [0]*len(fHistoParameter)
            for iratio in range(len(fHistoParameter)):
                h1ratio1[iratio] = fHistoParameter[iratio].Clone("h1ratio");
                h1ratio1[iratio].Reset();
                h1ratio1[iratio].Sumw2();
                h1ratio1[iratio].Divide(fHistoParameter[iratio], fHistoParameter[len(fHistoParameter)-1], 1., 1., "G");

            #finding minimum and max of the histogram
            yMin_ratio, yMax_ratio = self.set_y_min_y_max(h1ratio1, type_of_histo, scaling, ratio)

            frame2 = p2.DrawFrame(0., yMin_ratio, 20., yMax_ratio);
            if plotting == "RV0s" or plotting == "diff_cent_and_diff_RV0s":
                FrameSettings(frame2, "#it{p}_{T} (GeV/#it{c})", "#frac{#it{R}_{V^{0}}}{#it{R}_{V^{0}} all}",0., 20.)
            elif plotting == "centralities":
                FrameSettings(frame2, "#it{p}_{T} (GeV/#it{c})", "#frac{cent}{cent all}",0., 20.)
            elif plotting == "occupancies" or "diff_cent_and_diff_occupancies":
                FrameSettings(frame2, "#it{p}_{T} (GeV/#it{c})", "#frac{#it{FT0Occ}}{#it{FT0Occ} all}",0., 20.)

            set_frame(frame2, 0.10, 0.10, 1.0 , 0.6, 0.09, 0.09, 0.01, 0.01)
            frame2.GetYaxis().CenterTitle(True);
            ROOT.SetOwnership(frame2,False)

            for iratio in range(len(fHistoParameter)):
                h1ratio1[iratio].SetMarkerSize(markersize)
                if info_decay.get_typ() == "mc" and type_of_histo == "RawYield":
                    h1ratio1[iratio].SetMarkerStyle(kOpenCircle)
                else:
                    h1ratio1[iratio].SetMarkerStyle(kFullCircle)
                h1ratio1[iratio].SetLineColor(color[iratio])
                h1ratio1[iratio].SetMarkerColor(color[iratio])
                h1ratio1[iratio].DrawCopy("E1,same");

            line1 = TLine(0,1,20.,1);
            line1.SetLineColor(kBlack);
            line1.SetLineStyle(1);
            line1.SetLineWidth(1);
            line1.Draw("");
            ROOT.SetOwnership(line1,False);

        c1.cd()
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        txt = set_txt(0.1,0.9,1.,0.95,"NDC",0.045, 22)
        self.set_txt_on_canvas(txt, plotting, is_mc_and_data, type_of_histo, title, len(fHistoParameter), info_decay)

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(namePlot.Data());
        c1.Close();

        del c1;

    def plot_Rcp(self, all_cent, type_of_histo, fHistoParameter, namePlot, Period, ssnames, energy, system, scale_pT = False, err_stat_list = [], err_syst_list = [], is_other_file = False):

        TGaxis.SetMaxDigits(3);
        titlePt = fHistoParameter[0].GetTitle(); # ACCESS TITLE
        color = [kRed, 802, 401, kGreen+1, kGreen-1, kCyan ,kBlue, kMagenta, kMagenta+3, 922, kBlack]

        c1 = TCanvas("c0","c0",0,0,800,800);
        p1 = c1.cd()
        p1.SetPad(0,0.1,1,1);
        p1.SetMargin(0.10,0.1,0.1,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy();
        p1.SetLogx()

        y_title = self.set_y_title(type_of_histo, scale_pT)
        frame1 = p1.DrawFrame(0.4, 0.05, 10., 3.5);
        set_frame_basics(frame1, "p_{T} (GeV/c)", y_title, 0.4, 10.)
        set_frame(frame1, 0.045, 0.045, 1.05, 1.0, 0.045, 0.045, 0.001, 0.01)
        frame1.GetYaxis().SetMaxDigits(3);
        ROOT.SetOwnership(frame1,False)

        line1 = TLine(0,1,10.,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(1);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        leg = TLegend(0.15,0.68,0.35,0.83);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.02);
        leg.SetTextFont(42);#helvetica
        markersize = 1.5

        h1ratio1 = [0]*len(ssnames)
        print("ssname length ", len(ssnames))
        print("histo length ", len(fHistoParameter))

        for iratio in range(len(fHistoParameter)):
            self.set_subsystem(ssnames[iratio])
            self.set_centrality()
            if (self.cent1 == 80 and self.cent2 == 100) or (self.cent1 == 80 and self.cent2 == 90):
                peripheral = iratio
                N_coll_periph = self.set_N_coll()
                print("N_peripheral: ", N_coll_periph)
                break

        for iratio in range(peripheral+1):
            if is_other_file == False:
                fHistoParameter[peripheral].Sumw2(1)
                fHistoParameter[iratio].Sumw2(1)
            h1ratio1[iratio] = fHistoParameter[iratio].Clone("h1ratio");
            #h1ratio1[iratio].Reset();
            h1ratio1[iratio].SetMinimum(-3);
            h1ratio1[iratio].SetMaximum(+3);

            h1ratio1[iratio].SetMarkerSize(markersize)
            h1ratio1[iratio].SetMarkerStyle(kFullCircle)
            h1ratio1[iratio].SetLineColor(color[iratio])
            h1ratio1[iratio].SetMarkerColor(color[iratio])

            self.set_subsystem(ssnames[iratio])
            self.set_centrality()
            N_coll_cent = self.set_N_coll()
            print("N_cent: ", N_coll_cent)

            title = "#frac{ N_{coll}^{ 80-90% } yield^{one cent} }{ N_{coll}^{one cent} yield^{80-90%} }"

            h1ratio1[iratio].Divide(N_coll_periph*fHistoParameter[iratio], N_coll_cent*fHistoParameter[peripheral], 1., 1., "G");
            print("Here: ", is_other_file)
            if is_other_file:
                print("HERE1")
                h_err_stat_icent = err_stat_list[iratio]
                h_err_stat_periph = err_stat_list[peripheral]
                h_err_syst_icent = err_syst_list[iratio]
                h_err_syst_periph = err_syst_list[peripheral]
                nbin = h_err_stat_icent.GetNbinsX()+1
                for ibin in range(1, nbin):
                    arg_sq_1 = (h_err_stat_icent.GetBinContent(ibin)/fHistoParameter[iratio].GetBinContent(ibin))*(h_err_stat_icent.GetBinContent(ibin)/fHistoParameter[iratio].GetBinContent(ibin))
                    arg_sq_2 = (h_err_stat_periph.GetBinContent(ibin)/fHistoParameter[peripheral].GetBinContent(ibin))*(h_err_stat_periph.GetBinContent(ibin)/fHistoParameter[peripheral].GetBinContent(ibin))
                    error = h1ratio1[iratio].GetBinContent(ibin) * np.sqrt(arg_sq_1 + arg_sq_2)
                    #error_syst = h_err_syst_icent.GetBinContent(ibin)/h_err_syst_periph.GetBinContent(ibin)
                    #error = np.sqrt(error_stat*error_stat+error_syst*error_syst)
                    h1ratio1[iratio].SetBinError(ibin, error)
            h1ratio1[iratio].DrawCopy("E1,same");

            self.set_subsystem(ssnames[iratio])
            self.set_centrality()
            leg.AddEntry(h1ratio1[iratio],"{0}%-{1}%".format(self.cent1, self.cent2),"ep");

        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        txt = set_txt(0.,0.9,1.,0.95,"NDC",0.045, 22)
        txt.AddText("{0}".format(TString(title)));
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        txt = set_txt(0.7,0.74,1.0,0.84,"NDC",0.02, 12)
        txt.AddText("this thesis");
        txt.AddText(Period)
        txt.AddText("{0} at #sqrt{{#it{{s}}}} = {1}".format(system, energy));
        decay_name = set_decay_name(self.meson, self.decay)
        txt.AddText(decay_name)
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(namePlot.Data());
        c1.Close();

        del c1;

    def plot_pythia(self, graph, info_decay, name_plot):

        graph.SetLineWidth(5)
        graph.SetLineColor(8)

        #plot the pythia corrected yield
        canvas = TCanvas("pyt", "pyt", 0, 0, 900, 900)
        pad = canvas.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.2, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        graph.Draw("AC") #C is with a smooth line
        graph.GetYaxis().SetRangeUser(5e-7, 1.1)

        graph.GetXaxis().SetTitleSize(0.035)
        graph.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
        graph.GetYaxis().SetTitleOffset(2.0)
        graph.GetYaxis().SetTitleSize(0.035)
        y_title = self.set_y_title("CorrectedYield", False)
        graph.GetYaxis().SetTitle(y_title)

        # Define legend settings
        leg = TLegend(0.67, 0.6, 1.0, 0.7)
        leg.SetBorderSize(0)
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.025)
        leg.AddEntry(graph, "PYTHIA", "L")
        leg.Draw("")
        ROOT.SetOwnership(leg,False)

        # Add text
        set_txt_of_plot(info_decay, 0.87, 0.87, 0.87, 0.87, True, 0.027)

        #save histogram
        canvas.Modified()
        canvas.Update()
        ROOT.SetOwnership(canvas, False)
        canvas.SaveAs(name_plot.Data())

    def plot_pyth_vs_corr(self, pythia, corr, info_decay, name_plot, scaled_pT = False):

        max_h1 = corr.GetMaximum()
        max_y = 2 * max_h1

        corr.SetMaximum(max_y)

        c1 = ROOT.TCanvas("py_vs_cor", "py_vs_cor", 800, 600)
        pad = c1.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        pythia.SetLineColor(8)
        pythia.SetLineWidth(5)
        pythia.SetMaximum(10e+2)
        pythia.SetMarkerColor(8)
        pythia.SetMarkerStyle(47)
        pythia.Draw("AL")
        pythia.GetXaxis().SetTitleSize(0.035)
        pythia.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
        pythia.GetYaxis().SetTitleOffset(1.5)
        pythia.GetYaxis().SetTitleSize(0.035)
        y_title = self.set_y_title("CorrectedYield", scaled_pT)
        pythia.GetYaxis().SetTitle(y_title)

        corr.SetLineColor(ROOT.kBlue)
        corr.SetMarkerColor(ROOT.kBlue)
        corr.SetMaximum(10e+2)
        corr.SetMarkerStyle(20)
        corr.SetMarkerSize(1.2)
        corr.Draw("E1P Same")

        legend = ROOT.TLegend(0.78, 0.58, 0.88, 0.65)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.027)
        legend.AddEntry(corr, "data", "PE")
        legend.AddEntry(pythia, "PYTHIA", "L")
        legend.Draw()

        set_txt_of_plot(info_decay, 0.871, 0.87, 0.871, 0.87, True, 0.027)

        c1.Update()
        c1.SaveAs(name_plot.Data())


    def plot_fit_func(self, corr_hist_scaled, energy, info_decay, which_data, name_plot):
        #This is the fit with the Two-Component Model to the corrected yield

        c1 = ROOT.TCanvas("c_fit", "c_fit", 900,900)
        pad = c1.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.2, 0.1, 0.1, 0.1)
        pad.SetLogy()

        h_corr = corr_hist_scaled.Clone("corr_fit")
        h_corr.Scale(59.4*1e-3) #getting the differential invariant yield
        h_corr.Scale(1e12) #scaling to pico barn

        if which_data == "this_data":
            fit_fun = self.set_fit_fun(244*1e9, 0.135, 0.157, 27*1e9, 0.6, 2.96)
        else:
            fit_fun = self.set_fit_fun(536*1e9, 0.135, 0.142, 30*1e9, 0.63, 2.96)
        fit_fun.SetLineColor(ROOT.kOrange)

        h_corr.Fit(fit_fun,"TCM","",0.2, 15);
        y = h_corr.GetYaxis()
        y.SetTitle("E #frac{d^{3}#it{#sigma}}{d#it{p}^{3}} (pb GeV^{-2}#it{c}^{3})")
        y.SetTitleSize(0.035)
        y.SetTitleFont(42)
        y.SetTitleOffset(2.0)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)

        legend = ROOT.TLegend(0.52, 0.70, 0.72, 0.60)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.025)
        legend.AddEntry(h_corr, "corrected yield from data", "PE")
        legend.AddEntry(fit_fun, "Fit with two component model", "L")
        legend.Draw()

        #set_txt_system(energy, system, self.decay, cent1, cent2)
        set_txt_of_plot(info_decay, 0.871, 0.87, 0.871, 0.87, True, 0.025)
        c1.Update()
        c1.SaveAs(name_plot.Data())

        params = fit_fun.GetParameters()
        params_err = fit_fun.GetParErrors()
        for i in range(6):
            print(f"Parameter {i}: {params[i]}", "±", params_err[i])