# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3
# This code was written by Alica Enderich (Febuary 2024)
# This code was modified by Julia Schlägel (July 2024)

import numpy as np
import ROOT
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle
import re
import numpy as np
import math
import ctypes
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

from Histo_formatting import FrameSettings, DrawHisto, set_axis, set_txt, set_txt_system, set_decay_name, set_histo, set_frame, set_line, set_txt_of_plot

class Plot_efficiency:
    def __init__(self):
        print("#################################################")
        print("Plot_efficiency.py")
        print("#################################################")
        print("default constructor is called");

    def __del__(self):
        print("default destructor is called");

    def set_y_min_y_max(self, histo_list):
        yMax_, yMax_1 = -100, -100
        yMin_, yMin_1 = 100, 100
        for i in range(len(histo_list)-1):
            histo_list[i].BufferEmpty(-1)
            histo_list[i+1].BufferEmpty(-1)
            Nbins = histo_list[i].GetNbinsX()+1
            for ibin in range(Nbins):
                if yMax_ < histo_list[i].GetBinContent(ibin) and histo_list[i].GetBinContent(ibin) != 0:
                    yMax_ = histo_list[i].GetBinContent(ibin)
                if yMin_ > histo_list[i].GetBinContent(ibin) and histo_list[i].GetBinContent(ibin) > 0:
                    yMin_ = histo_list[i].GetBinContent(ibin)
            Nbins1 = histo_list[i+1].GetNbinsX()+1
            for ibin in range(Nbins1):
                if yMax_1 < histo_list[i+1].GetBinContent(ibin) and histo_list[i+1].GetBinContent(ibin) != 0:
                    yMax_1 = histo_list[i+1].GetBinContent(ibin)
                if yMin_1 > histo_list[i+1].GetBinContent(ibin) and histo_list[i+1].GetBinContent(ibin) > 0:
                    yMin_1 = histo_list[i+1].GetBinContent(ibin)
            if yMax_1 > yMax_:
                yMax_ = yMax_1
            if yMin_1 < yMin_:
                yMin_ = yMin_1
        yMax_ = 1.2*yMax_
        yMin_ = 0.4*yMin_
        return yMin_, yMax_

    def plot_eff_acc(self, h_eff_acc, name_plot, name, info_decay, teff):
        #get canvas and pad to plot the acceptance in
        can = TCanvas("eff_acc", "eff_acc", 0, 0, 900, 900)
        pad = can.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        #Define marker settings
        h_eff_acc.SetFillColor(kBlue)
        h_eff_acc.SetMarkerStyle(20)
        h_eff_acc.SetMarkerColor(kBlue)
        h_eff_acc.SetMarkerSize(1.5)
        h_eff_acc.SetLineColor(kBlue)

        #Draw points in the histogramm
        h_eff_acc.Draw()
        if teff:
            h_eff_acc.SetTitle("TEfficiency; #it{p}_{T} (GeV/#it{c}) ; Teff x A x BR")
            can.Update()
        else:
            # Define y-axis settings
            y = h_eff_acc.GetYaxis()
            if "BR" in name:
                set_axis(y, "#varepsilon x A", 0.032, 42, 1.7, 0.032, 42)
            else:
                set_axis(y, "#varepsilon x A x BR", 0.032, 42, 1.7, 0.032, 42)

            #Define x-axis settings
            x = h_eff_acc.GetXaxis()
            set_axis(x, "#it{p}_{T} (GeV/#it{c})", 0.032, 42, 1.4, 0.032, 42)

        # Add text
        set_txt_of_plot(info_decay, 0.87, 0.87, 0.87, 0.87, False, 0.025)

        # Add legend
        leg = TLegend(0.66,0.51,0.78,0.65)
        leg.SetTextSize(0.025);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextFont(42)
        if info_decay.get_plotting_sys() == "centralities":
            cent1, cent2 = info_decay.get_cent()
            leg.AddEntry(h_eff_acc,"{0}%-{1}%".format(cent1, cent2),"ep")
        if info_decay.get_plotting_sys() == "RV0s":
            V0_radius_min, V0_radius_max = info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max()
            leg.AddEntry(h_eff_acc,"{0} < #it{{R}}_{{V^{{0}}}} < {1} cm".format(int(V0_radius_min), int(V0_radius_max)),"EP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        #save histogram
        can.Modified()
        can.Update()
        ROOT.SetOwnership(can, False)
        can.SaveAs(name_plot.Data())

    def plot_eff_acc_all_syst(self, histo_list, name_plot, info_decay, name, plotting, ratio):
        TGaxis.SetMaxDigits(3);
        ssname_list = []
        if plotting == "centralities" or plotting == "RV0s":
            ssname_list = info_decay.get_ssname_list()
        if plotting == "diff_cent_and_diff_RV0s":
            RV0_list = info_decay.get_RV0_list()
            ssname_list.append(info_decay.get_ssname())
        if plotting == "diff_cent_and_diff_occupancies":
            occup_list = info_decay.get_occupancy_list()
            ssname_list.append(info_decay.get_ssname())

        yMin_, yMax_ = self.set_y_min_y_max(histo_list)

        c1 = TCanvas("c0","c0",0,0,800,800);
        if ratio:
            c1.Divide(1,2,1e-3,1e-3);
            p1 = c1.cd(1);
            p1.SetPad(0,0.3, 1,1)
            p1.SetMargin(0.15,0.02,0.,0.15);
            p1.SetBottomMargin(0.005)
        else:
            p1 = c1.cd()
            p1.SetPad(0,0.1,2,2);
            p1.SetMargin(0.18,0.1,0.1,0.15);
            #p1.SetTicks(1,1);
        p1.SetLogy();

        frame1 = p1.DrawFrame(0., yMin_, 20., yMax_);
        if name == "Efficiency_Acceptance":
            FrameSettings(frame1, "#it{p}_{T} (GeV/#it{c})", "#varepsilon x A x BR", 0., 20.)
        if name == "Efficiency_Acceptance_scaled_by_BR":
            FrameSettings(frame1, "#it{p}_{T} (GeV/#it{c})", "#varepsilon x A x BR scaled",0., 20.)
        if name == "T_Efficiency_Acceptance_shifted":
            FrameSettings(frame1, "#it{p}_{T} (GeV/#it{c})", "TEfficiency",0., 20.)
        frame1.GetYaxis().SetMaxDigits(3)
        set_frame(frame1, 0.045, 0.045, 1.0, 1.52, 0.045, 0.045, 0.01, 0.01)
        ROOT.SetOwnership(frame1,False)

        if plotting == "centralities":
            leg = TLegend(0.74, 0.39, 0.84, 0.58)
        if plotting == "RV0s" or plotting == "diff_cent_and_diff_RV0s":
            leg = TLegend(0.65, 0.45, 1.10, 0.63)
        if plotting == "occupancies" or plotting == "diff_cent_and_diff_occupancies":
            leg = TLegend(0.62, 0.45, 1.07, 0.63)
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.025);
        leg.SetTextFont(42);#helvetica
        color = [kRed, 802, 401, kGreen+1, kGreen-1, kCyan ,kBlue, kMagenta, kMagenta+3, 922, 923, kBlack]
        for icut in range(len(histo_list)):
            set_histo(histo_list[icut], 0.025, 0.025, 0.05, 0.05)
            DrawHisto(histo_list[icut], color[icut], "E1,same", 1.5, kFullCircle)
            if plotting != "diff_cent_and_diff_RV0s" and plotting != "diff_cent_and_diff_occupancies":
                info_decay.set_centrality(ssname_list[icut])
                cent1, cent2 = info_decay.get_cent()
                info_decay.set_RV0(ssname_list[icut])
                V0_radius_min, V0_radius_max = info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max()
            else:
                info_decay.set_centrality(ssname_list[0])
                cent1, cent2 = info_decay.get_cent()
                info_decay.set_RV0(ssname_list[0])
                V0_radius_min, V0_radius_max = info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max()
            if plotting == "centralities" and info_decay.get_system()!="pp":
                leg.AddEntry(histo_list[icut],"{0}%-{1}%".format(cent1, cent2),"ep")
            if plotting == "RV0s" or plotting == "diff_cent_and_diff_RV0s":
                if plotting == "diff_cent_and_diff_RV0s":
                    if icut < len(RV0_list)-1:
                        V0_radius_min = RV0_list[icut]
                        V0_radius_max = RV0_list[icut+1]
                    else:
                        V0_radius_min = RV0_list[0]
                        V0_radius_max = RV0_list[len(RV0_list)-1]
                leg.AddEntry(histo_list[icut],"{0} < #it{{R}}_{{V^{{0}}}} < {1} cm".format(int(V0_radius_min), int(V0_radius_max)),"EP");
            if plotting == "occupancies" or plotting == "diff_cent_and_diff_occupancies":
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
                leg.AddEntry(histo_list[icut],"{0}".format(occup_range_latex),"ep")

        if ratio:
            p2 = c1.cd(2);
            p2.SetPad(0,0,1,0.3);
            p2.SetMargin(0.15,0.02,0.22,0.0);
            p2.SetTicks(1,1)
            p2.SetTopMargin(0.05)
            p2.SetBottomMargin(0.3)
            #p2.SetLogy();

            line1 = TLine(0,1,20.,1);
            line1.SetLineColor(kBlack);
            line1.SetLineStyle(1);
            line1.SetLineWidth(1);
            line1.Draw("");
            ROOT.SetOwnership(line1,False);

            h1ratio1 = [0]*len(histo_list)
            for iratio in range(len(histo_list)):
                h1ratio1[iratio] = histo_list[iratio].Clone("h1ratio");
                h1ratio1[iratio].Reset();
                h1ratio1[iratio].Sumw2();
                h1ratio1[iratio].Divide(histo_list[iratio], histo_list[len(histo_list)-1], 1., 1., "G");

            #finding minimum and max of the histogram
            yMin_ratio, yMax_ratio = self.set_y_min_y_max(h1ratio1)

            frame2 = p2.DrawFrame(0., yMin_ratio, 20., yMax_ratio);
            if plotting == "RV0s" or plotting == "diff_cent_and_diff_RV0s":
                FrameSettings(frame2, "#it{p}_{T} (GeV/#it{c})", "#frac{#it{R}_{V^{0}}}{#it{R}_{V^{0}} all}",0., 20.)
            elif plotting == "centralities":
                FrameSettings(frame2, "#it{p}_{T} (GeV/#it{c})", "#frac{cent}{cent all}",0., 20.)
            elif plotting == "occupancies" or plotting == "diff_cent_and_diff_occupancies":
                FrameSettings(frame2, "#it{p}_{T} (GeV/#it{c})", "#frac{#it{FT0Occ}}{#it{FT0Occ all}}",0., 20.)
            frame2.GetXaxis().SetTitleSize(0.10);
            frame2.GetYaxis().SetTitleSize(0.10);
            frame2.GetXaxis().SetTitleOffset(1.0);
            frame2.GetYaxis().SetTitleOffset(0.7);
            frame2.GetXaxis().SetLabelSize(0.10);
            frame2.GetYaxis().SetLabelSize(0.10);
            frame2.GetYaxis().CenterTitle(True);
            frame2.GetXaxis().SetLabelOffset(0.01);
            frame2.GetYaxis().SetLabelOffset(0.01);
            ROOT.SetOwnership(frame2,False);
            for iratio in range(len(histo_list)):
                #print(f"iratio {iratio}")
                h1ratio1[iratio].SetMarkerSize(1.5)
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

        txt = set_txt(0.,0.9,1.,0.95,"NDC", 0.045, 22)
        if name == "Efficiency_Acceptance":
            txt.AddText("#varepsilon x A x BR as a function of #it{p}_{T}")
        if name == "Efficiency_Acceptance_scaled_by_BR":
            txt.AddText("#varepsilon x A as a function of #it{p}_{T}");
        if name == "T_Efficiency_Acceptance_shifted":
            txt.AddText("TEfficiency as a function of #it{p}_{T}")
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        if plotting == "centralities":
            set_txt_of_plot(info_decay, 0.87, 0.83, 0.87, 0.83, False, 0.025)
        else:
            set_txt_of_plot(info_decay, 0.95, 0.88, 0.95, 0.88, False, 0.025)

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(name_plot.Data());
        c1.Close();

    def comp_eff_Teff(self, h_eff_acc, h_Teff_acc, name_plot, info_decay):
        max_h_eff_acc = h_eff_acc.GetMaximum()
        max_y = 1.5 * max_h_eff_acc
        h_eff_acc.SetMaximum(max_y)

        c1 = ROOT.TCanvas("eff_vs_teff", "eff_vs_teff", 0, 0, 800, 600)
        pad = c1.cd()

        pad.SetPad(0.0, 0.0, 1, 1)
        pad.SetMargin(0.15,0.1,0.1,0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        DrawHisto(h_Teff_acc, ROOT.kMagenta, "E1P", 1.8, 22)

        h_eff_acc.SetTitle("Comparison efficiency vs Teff")
        h_eff_acc.SetXTitle("p_{T} (GeV/c)")
        h_eff_acc.SetYTitle("Efficiency x acc x BR")
        h_eff_acc.GetXaxis().SetRangeUser(0,12)
        DrawHisto(h_eff_acc, ROOT.kBlue, "E1P Same", 1.8, 21)

        legend = ROOT.TLegend(0.7, 0.65, 0.9, 0.55)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        legend.AddEntry(h_eff_acc, "eff x acc x BR", "LP")
        legend.AddEntry(h_Teff_acc, "Teff x acc x BR", "LP")
        legend.Draw()

        # Add text
        cent1, cent2 = info_decay.get_cent()
        set_txt_system(info_decay.get_energy(), info_decay.get_system(), info_decay.get_decay(), cent1, cent2)

        c1.Update()

        c1.Update()
        c1.SaveAs(name_plot.Data())

