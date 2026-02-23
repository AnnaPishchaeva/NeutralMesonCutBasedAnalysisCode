# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3
# This code was written by Alica Enderich (Febuary 2024)
# This code was modified and extended by Julia Schlägel (July 2024)
# This code was modified and extended by Anna Pishchaeva (October 2025)

import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX
import numpy as np
import yaml
import array
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine
from Histo_formatting import *

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)


class Compare:
    def __init__(self) -> None:
        pass
    def __init__(self, decay, meson, cent1, cent2):
        print("#################################################")
        print("Compare.py")
        print("#################################################")
        self.meson = meson
        self.decay = decay
        self.cent1 = cent1
        self.cent2 = cent2

    def set_y_min_y_max(self, histo_list, ratio = False):
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
        if not ratio:
            yMax_ = 1.2*yMax_
            yMin_ = 0.4*yMin_
        else:
            yMax_ = 1.01*yMax_
            yMin_ = 0.1*yMin_
        return yMin_, yMax_

    def set_title(self, name, h1, V0_min, V0_max):
        h1.SetXTitle("p_{T} (GeV/c)")
        if name == "acc":
            title = "Acceptance centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "Acceptance"
            h1.SetYTitle("Acceptance")
        if name == "eff":
            title = "Efficiency centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "Efficiency"
            h1.SetYTitle("Efficiency")
        if name == "corr_yield":
            title = "Corrected yield centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "#frac{1}{2#pi} #frac{d^{2}N}{dp_{T}dy} (GeV/#it{c})^{-1}"
            h1.SetYTitle("#frac{1}{2#pi} #frac{d^{2}N}{dp_{T}dy} (GeV/#it{c})^{-1}")
        if name == "inv_corr_yield":
            title = "Invariant corrected yield centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "#frac{1}{2#pip_{T}} #frac{d^{2}N}{dp_{T}dy} [GeV^{-1} c^{3}]"
            h1.SetYTitle("#frac{1}{2#pip_{T}} #frac{d^{2}N}{dp_{T}dy} [GeV^{-1} c^{3}]")
        if name == "eff_acc":
            title = "#varepsilon x A x BR centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "#varepsilon x A x BR"
            h1.SetYTitle("#varepsilon x A x BR")
        if name == "eff_acc_scaled":
            title = "#varepsilon x A centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "#varepsilon x A"
            h1.SetYTitle("#varepsilon x A")
        if name == "eff_acc_all_system":
            title = "#varepsilon x A x BR"
            y_title = "#varepsilon x A x BR"
            h1.SetYTitle("#varepsilon x A x BR")
        if name == "inv_cross_section":
            title = "Invariant cross section centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "#frac{1}{2#pip_{T}} #frac{d^{2}#sigma}{dp_{T}dy} [GeV^{-2} c^{3}]"
            h1.SetYTitle("#frac{1}{2#pip_{T}} #frac{d^{2}#sigma}{dp_{T}dy} [GeV^{-2} c^{3}]")
        if name == "corr_ratio":
            title = "#eta / #pi^{{0}} ratio centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "#eta / #pi^{0}"
            h1.SetYTitle("#eta / #pi^{0}")
        if name == "mean":
            title = "mean for centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "mean (GeV/c^{2})"
            h1.SetYTitle("mean (GeV/c^{2})")
        if name == "lambda":
            title = "#lambda for centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "#lambda (GeV/c^{2})"
            h1.SetYTitle("#lambda (GeV/c^{2})")
        if name == "sigma":
            title = "#sigma for centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "#sigma (GeV/c^{2})"
            h1.SetYTitle("#sigma (GeV/c^{2})")
        if name == "eta_pi0_ratio":
            title = "#eta\#{pi}^{{0}} for centrality {0}%-{1}%".format(self.cent1, self.cent2)
            y_title = "#eta\#{pi}^{{0}}"
            h1.SetYTitle("#eta\#{pi}^{{0}}")
        if name == "corr_diff_pkl":
            title = "Corrected yield of different systems {0}<R<{1} cm".format(V0_min, V0_max)
            y_title = "#frac{1}{2#pi} #frac{d^{2}N}{dp_{T}dy} (GeV/#it{c})^{-1}"
            h1.SetYTitle("#frac{1}{2#pi} #frac{d^{2}N}{dp_{T}dy} (GeV/#it{c})^{-1}")
        if name == "eff_acc_diff_pkl":
            title = "#varepsilon x A x BR of different systems {0}<R<{1} cm".format(V0_min, V0_max)
            y_title = "#varepsilon x A x BR"
            h1.SetYTitle("#varepsilon x A x BR")
        if name == "raw_diff_data_pkl":
            title = "Raw yield of different systems {0}<R<{1} cm Data".format(V0_min, V0_max)
            y_title = "#frac{1}{#it{N}_{ev}}#frac{d#it{N}^{raw}}{d#it{p}_{T}} (GeV/#it{c})^{-1}"
            h1.SetYTitle("#frac{1}{#it{N}_{ev}}#frac{d#it{N}^{raw}}{d#it{p}_{T}} (GeV/#it{c})^{-1}")
        if name == "raw_diff_mc_pkl":
            title = "Raw yield of different systems {0}<R<{1} cm MC".format(V0_min, V0_max)
            y_title = "#frac{1}{#it{N}_{ev}}#frac{d#it{N}^{raw}}{d#it{p}_{T}} (GeV/#it{c})^{-1}"
            h1.SetYTitle("#frac{1}{#it{N}_{ev}}#frac{d#it{N}^{raw}}{d#it{p}_{T}} (GeV/#it{c})^{-1}")
        if name == "raw_yields_diff_periods":
            title = "Raw yields of different settings {0} - {1}%".format(self.cent1, self.cent2)
            y_title = "#frac{1}{#it{N}_{ev}}#frac{d#it{N}^{raw}}{d#it{p}_{T}} (GeV/#it{c})^{-1}"
            h1.SetYTitle("#frac{1}{#it{N}_{ev}}#frac{d#it{N}^{raw}}{d#it{p}_{T}} (GeV/#it{c})^{-1}")
        return title, y_title

    def set_legend_name(comparison, legend):
        if comparison == "pcm_dalitz":
            if self.meson == "pi0":
                legend.AddEntry(h1, "#pi^{0} -> e^{+} e^{-} #gamma", "LP")
                legend.AddEntry(h2, "#pi^{0} -> #gamma #gamma", "LP")
            if self.meson == "eta":
                legend.AddEntry(h1, "#eta -> e^{+} e^{-} #gamma", "LP")
                legend.AddEntry(h2, "#eta -> #gamma #gamma", "LP")
        if comparison == "pi0_eta":
            if self.decay == "dalitz":
                legend1.AddEntry(h1, "#eta -> e^{+} e^{-} #gamma", "LP")
                legend1.AddEntry(h2, "#pi^{0} -> e^{+} e^{-} #gamma", "LP")
            if self.decay == "pcm":
                legend1.AddEntry(h1, "#eta -> #gamma #gamma", "LP")
                legend1.AddEntry(h2, "#pi^{0} -> #gamma #gamma", "LP")

    def get_graph_y_range(graph):
        n_points = graph.GetN()
        y_values = [graph.GetY()[i] for i in range(n_points)]
        y_min = min(y_values)
        y_max = max(y_values)
        return y_min, y_max

    def plot_histo_one_vs_another_w_ratio(self, h1, h2, name_plot, energy, name, comparison, V0_min, V0_max, system, period, period_mc = 0):
        TGaxis.SetMaxDigits(3);

        max_h1 = h1.GetMaximum()
        max_h2 = h2.GetMaximum()
        max_y = 1.0 * max(max_h1, max_h2)  # 2 times higher than the biggest point
        #h1.SetMaximum(max_y)
        if name == "corr_yield":
            #h1.GetYaxis().SetRangeUser(10e-7, 10e2)
            y_min, y_max = 1e-7, 1e1
        if name == "eff_acc":
            y_min, y_max = 5e-8, 5e-4
        #if name == "mean" and self.meson == "pi0":
        #    y_min, y_max = 120, 150
        #if name == "lambda" and self.meson == "pi0":
        #    y_min, y_max = 0, 3*1e-2
        #if name == "sigma" and self.meson == "pi0":
        #    y_min, y_max = 0, 1e-2
        if name == "mean" or name == "lambda" or name == "sigma":
            y_max1, y_max = -100, -100
            y_min1, y_min = 100, 100
            Nbins = h1.GetNbinsX()+1
            for ibin in range(Nbins):
                if y_max < h1.GetBinContent(ibin) and h1.GetBinContent(ibin) > 0:
                    y_max = h1.GetBinContent(ibin)
                if y_min > h1.GetBinContent(ibin) and h1.GetBinContent(ibin) > 0:
                    y_min = h1.GetBinContent(ibin)
            Nbins1 = h2.GetNbinsX()+1
            for ibin in range(Nbins1):
                if y_max1 < h2.GetBinContent(ibin) and h2.GetBinContent(ibin) > 0:
                    y_max1 = h2.GetBinContent(ibin)
                if y_min1 > h2.GetBinContent(ibin) and h2.GetBinContent(ibin) > 0:
                    y_min1 = h2.GetBinContent(ibin)
            if y_max1 > y_max:
                y_max = y_max1
            if y_min1 < y_min:
                y_min = y_min1
            y_max = 1.01*y_max
            y_min = 0.98*y_min
            print("y_min ", y_min, "name: ", name)

        x_min, x_max = 0, 21
        # Make a Canvas
        c1 = ROOT.TCanvas("c1", "c1", 0, 0, 800, 800)
        c1.Divide(1, 2, 1e-3,1e-3)

        p1 = c1.cd(1)
        p1.SetPad(0.0, 0.3, 1, 1)
        p1.SetMargin(0.15,0.02,0.,0.15)
        #p1.SetTicks(1,1)
        if name!= "mean" and name!="sigma" and name!="lambda":
            p1.SetLogy()

        title, y_title = self.set_title(name, h1, V0_min, V0_max)
        frame1 = p1.DrawFrame(x_min, y_min, x_max, y_max);
        set_frame_basics(frame1, "p_{T} (GeV/c)", y_title, x_min, x_max)
        set_frame(frame1, 0.045, 0.045, 1.0, 1.52, 0.045, 0.045, 0.01, 0.01)
        frame1.GetYaxis().SetMaxDigits(3);
        #frame1.GetYaxis().SetNdivisions(510)
        ROOT.SetOwnership(frame1,False)
        h1.GetYaxis().SetTitleSize(0.05)
        h1.GetYaxis().SetLabelSize(0.05)
        if name == "mean":
            h1.SetMinimum(y_min)
            h1.SetMaximum(y_max)
        h1.GetYaxis().SetTitleOffset(1.5)

        if name == "mean" or name == "sigma" or name == "lambda":
            DrawHisto(h1, ROOT.kBlue, "E1P", 1.5, kFullCircle)
            DrawHisto(h2, ROOT.kBlue, "E1P same", 1.5, kOpenCircle)
        else:
            DrawHisto(h1, ROOT.kOrange, "E1P", 1.5, 20)
            DrawHisto(h2, ROOT.kViolet, "E1P same", 1.5, 34)

        if name == "mean":
            if self.meson == "pi0":
                line = set_line(0.1,0.13497,20,0.13497, kBlack, 2, 1)
            if self.meson == "eta":
                line = set_line(0.1,0.54786,20,0.54786, kBlack, 2, 1)
            line.SetLineStyle(2)
            line.Draw("");
            ROOT.SetOwnership(line,False);

        p2 = c1.cd(2)
        p2.SetPad(0,0,1,0.3)
        p2.SetMargin(0.15,0.02,0.22,0.0)
        p2.SetTicks(2,2)
        #p2.SetLogy()
        if comparison == "pi0_eta":
            ratio_title = "#frac{#eta}{#pi^{0}}"
        if comparison == "pcm_dalitz":
            ratio_title ="#frac{dalitz}{pcm}"
        if comparison == "other_data":
            ratio_title = "#frac{5.02 TeV}{5.36 TeV}"
        if comparison == "mc_data":
            ratio_title = "#frac{MC}{Data}"
        if name == "lambda" or name == "sigma":
            y_min_frame2 = 0.1
            y_max_frame2 = 2.0
        elif name == "mean":
            y_min_frame2 = 0.8
            y_max_frame2 = 1.2
        else:
            y_min_frame2 = 0.01
            y_max_frame2 = 5.5
        frame2 = p2.DrawFrame(x_min, y_min_frame2, x_max, y_max_frame2); #0.,0.8,15.,2.0
        set_frame_basics(frame2, "p_{T} (GeV/c)", "{0}".format(ratio_title), x_min, x_max)
        set_frame(frame2, 0.1, 0.1, 1.0, 0.7, 0.1, 0.1, 0.01, 0.01)
        frame2.GetYaxis().CenterTitle(True);
        ROOT.SetOwnership(frame2,False);

        line1 = TLine(0, 1, 21.,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(1);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        h1.Sumw2()
        h2.Sumw2()
        ratio = h2.Clone("Ratio")
        ratio.Divide(h1)
        min_ratio = ratio.GetMinimum()
        ratio.SetLineColor(ROOT.kBlue)
        ratio.SetMarkerColor(ROOT.kBlue)
        ratio.SetMarkerStyle(kFullCircle)
        ratio.Draw("E1P, same")

        c1.cd()

        txt = set_txt(0.05,0.9,1.05,0.95,"NDC",0.045, 22)
        txt.AddText(title)
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        # Add legend
        legend = ROOT.TLegend(0.75, 0.8, 0.85, 0.7)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.025)
        if comparison == "pi0_eta":
            legend.AddEntry(h1, "#pi^{0}", "P")
            legend.AddEntry(h2, "#eta", "P")
        if comparison == "pcm_dalitz":
            legend.AddEntry(h2, "e^{+} e^{-} #gamma", "P")
            legend.AddEntry(h1, "#gamma #gamma", "P")
        if comparison == "mc_data":
            legend.AddEntry(h2, "MC {0}".format(period_mc), "P")
            legend.AddEntry(h1, "Data {0}".format(period), "P")
        legend.Draw()

        # Add text
        txt = set_txt(0.85, 0.85, 0.95, 0.85,"NDC", 0.025, 33)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = set_txt(0.85,0.825,0.95,0.825,"NDC", 0.025, 33)
        txt3.AddText("{0}, #sqrt{{s}} = {1}".format(system, energy))
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        c1.Update()
        c1.SaveAs(name_plot.Data())


    def plot_histo_one_vs_another(self, h1, h2, name_plot, system, energy, name, comparison, V0_min, V0_max):
        #get canvas and pad to plot the acceptance in
        can = TCanvas("c1", "c1", 0, 0, 900, 900)
        pad = can.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        title, y_title = self.set_title(name, h1, V0_min, V0_max)

        #Draw points in the histogramm
        DrawHisto(h1, ROOT.kOrange, "E1P", 1.8, 20)
        DrawHisto(h2, ROOT.kViolet, "E1P same", 1.8, 34)
        x = h1.GetXaxis()
        set_axis(x, "p_{T} [GeV/c]", 0.035, 42, 1.3, 0.035, 42)
        y = h1.GetYaxis()
        set_axis(y, "{0}".format(y_title), 0.035, 42, 1.5, 0.035, 42)
        can.Update()

        # Define legend settings
        legend = ROOT.TLegend(0.75, 0.8, 0.85, 0.7)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        if comparison == "pi0_eta":
            legend.AddEntry(h1, "#pi^{0}", "LP")
            legend.AddEntry(h2, "#eta", "LP")
        if comparison == "pcm_dalitz":
            legend.AddEntry(h2, "e^{+} e^{-} #gamma", "LP")
            legend.AddEntry(h1, "#gamma #gamma", "LP")
        legend.Draw()

        # Add text
        txt = set_txt(0.05,0.9,1.05,0.95,"NDC",0.045, 22)
        txt.AddText(title)
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt = set_txt(0.75, 0.85, 0.85, 0.85,"NDC", 0.025, 33)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = set_txt(0.75,0.825,0.85,0.825,"NDC", 0.025, 33)
        txt3.AddText("{0}, #sqrt{{s}} = {1}".format(system, energy))
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        #save histogram
        can.Modified()
        can.Update()
        ROOT.SetOwnership(can, False)
        can.SaveAs(name_plot.Data())

    def plot_ratio(self, ratio, name_plot, system, energy, name, comparison, V0_min, V0_max):
        #get canvas and pad to plot the acceptance in
        can = TCanvas("c1", "c1", 0, 0, 900, 900)
        pad = can.cd()
        pad.SetPad(0.0, 0.01, 1, 1)
        pad.SetMargin(0.15, 0.1, 0.1, 0.1)
        pad.SetTicks(1,1)
        pad.SetLogy()

        title, y_title = self.set_title(name, ratio, V0_min, V0_max)
        if comparison == "eta_pi0":
            y_title = "#frac{#eta}{#pi^{0}}"
        if comparison == "pcm_dalitz":
            y_title = "#frac{dalitz}{pcm}"

        #Draw points in the histogramm
        ratio.SetMaximum(5)
        DrawHisto(ratio, ROOT.kRed, "E1P", 1.8, 20)
        x = ratio.GetXaxis()
        set_axis(x, "p_{T} [GeV/c]", 0.035, 42, 1.3, 0.035, 42)
        y = ratio.GetYaxis()
        set_axis(y, "{0}".format(y_title), 0.035, 42, 1.5, 0.035, 42)
        can.Update()

        # Define legend settings
        legend = ROOT.TLegend(0.75, 0.8, 0.85, 0.7)
        legend.SetBorderSize(0)
        legend.SetFillColor(kWhite)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        if comparison == "pi0_eta":
            legend.AddEntry(ratio, "#frac{#eta}{#pi^{0}}", "LP")
        if comparison == "pcm_dalitz":
            legend.AddEntry(ratio, "#frac{dalitz}{pcm}", "LP")
        if comparison == "other_data":
            legend.AddEntry(ratio, "#frac{5.02 TeV}{5.36 TeV}", "LP")
            line1 = TLine(0,1,20.,1);
            line1.SetLineColor(kBlack);
            line1.SetLineStyle(1);
            line1.SetLineWidth(1);
            line1.Draw("");
            ROOT.SetOwnership(line1,False);
        legend.Draw()

        # Add text
        txt = set_txt(0.05,0.9,1.05,0.95,"NDC",0.045, 22)
        txt.AddText(title)
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt = set_txt(0.75, 0.85, 0.85, 0.85,"NDC", 0.025, 33)
        txt.AddText("this thesis")
        txt.Draw()
        ROOT.SetOwnership(txt,False)

        txt3 = set_txt(0.75,0.825,0.85,0.825,"NDC", 0.025, 33)
        txt3.AddText("{0}, #sqrt{{s}} = {1}".format(system, energy))
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

        #save histogram
        can.Modified()
        can.Update()
        ROOT.SetOwnership(can, False)
        can.SaveAs(name_plot.Data())

    def comp_with_other_data(self, histo_list, name, name_plot, energy_other_file,system_other_file, info_decay, ratio = False):
        print("in", flush = True)
        TGaxis.SetMaxDigits(3);
        V0_min, V0_max = info_decay.get_V0_radius_min(), info_decay.get_V0_radius_max()
        title, y_title = self.set_title(name, histo_list[0], V0_min, V0_max)
        yMin_, yMax_ = self.set_y_min_y_max(histo_list)
        histo_list[0].SetMaximum(yMax_)
        histo_list[0].SetMinimum(yMin_)
        print("almost", flush = True)
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
        FrameSettings(frame1, "#it{p}_{T} (GeV/#it{c})", "{}".format(y_title), 0., 20.)
        set_frame(frame1, 0.045, 0.045, 1.0, 1.52, 0.045, 0.045, 0.01, 0.01)
        ROOT.SetOwnership(frame1,False)

        if name == "corr_diff_pkl" or name == "eff_acc_diff_pkl" or name == "raw_diff_data_pkl" or name == "raw_diff_mc_pkl":
            legend_names = ["Pb--Pb cent {0}%-{1}% #sqrt{{s_{{NN}}}} = 5.36 TeV".format(self.cent1, self.cent2), "O--O cent 0%-100% #sqrt{{s}} = 5.36 TeV", "pO cent 0%-100%"]
        else:
            system = info_decay.get_system()
            if '-' in system:
                list_system = system.split("-")
                system = list_system[0] + "#font[122]{-}" + list_system[1]
            if '-' in system_other_file:
                list_system_other_file = system_other_file.split("-")
                system_other_file = list_system_other_file[0] + "#font[122]{-}" + list_system_other_file[1]
            legend_names = ["#sqrt{{s_{{NN}}}} = {0}, {1}".format(info_decay.get_energy(), system), "#sqrt{{s_{{NN}}}} = {0}, {1}".format(energy_other_file, system_other_file)]
        if name =="raw_yields_diff_periods":
            #legend_names = ["EM_LHC23", "EM_LHC24ar", "EM_LHC25an"]
            legend_names = ["EM_LHC23", "EM_LHC23_new_settings"]
        if name == "eff_acc_all_system":
            legend_names = ["OO, #sqrt{s} = 5.36 TeV", "pO, #sqrt{s} = 9.61 TeV", "pp, #sqrt{s} = 5.36 TeV", "{0}, #sqrt{{s_{{NN}}}} = {1}, 0-20%".format(system, info_decay.get_energy()), "{0}, #sqrt{{s_{{NN}}}} = {1}, 0-20%".format(energy_other_file, system_other_file), "{0}, #sqrt{{s_{{NN}}}} = {1}, 20-40%".format(system, info_decay.get_energy()), "{0}, #sqrt{{s_{{NN}}}} = {1}, 20-40%".format(energy_other_file, system_other_file), "{0}, #sqrt{{s_{{NN}}}} = {1}, 40-60%".format(system, info_decay.get_energy()), "{0}, #sqrt{{s_{{NN}}}} = {1}, 40-60%".format(energy_other_file, system_other_file), "{0}, #sqrt{{s_{{NN}}}} = {1}, 60-80%".format(system, info_decay.get_energy()),"{0}, #sqrt{{s_{{NN}}}} = {1}, 60-80%".format(energy_other_file, system_other_file)]
        legend1 = ROOT.TLegend(0.61, 0.59, 0.71, 0.49)
        legend1.SetBorderSize(0)
        legend1.SetFillColor(kWhite)
        legend1.SetFillStyle(0)
        legend1.SetTextSize(0.025)

        #color = [kRed-4, kGreen-3,kBlue-6, kMagenta, kMagenta+3, 922, 923, kBlack]
        color = [kRed, 802, 401, kGreen+1, kGreen-1, kCyan ,kBlue, kMagenta, kMagenta+3, 922, 923, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack]
        count = 0

        print("HERE!", flush = True)
        print(len(color), flush = True)
        print(len(legend_names), flush = True)
        print(len(histo_list), flush = True)
        for i in range (len(histo_list)):
            print("start", i, flush = True)
            print(histo_list[i], flush = True)
            histo_list[i].SetYTitle(y_title)
            histo_list[i].SetLineColor(color[i])
            histo_list[i].SetMarkerColor(color[i])
            histo_list[i].SetMarkerStyle(kFullCircle)
            histo_list[i].SetMarkerSize(1.3)
            histo_list[i].Draw("E1P Same")
            legend1.AddEntry(histo_list[i], legend_names[i], "P")
            print("finish", i, flush = True)
        legend1.Draw()
        ROOT.SetOwnership(legend1,False);

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
                h1ratio1[iratio].Sumw2();
                h1ratio1[iratio].Divide(histo_list[iratio], histo_list[0], 1., 1., "G")

            #finding minimum and max of the histogram
            yMin_ratio, yMax_ratio = self.set_y_min_y_max(h1ratio1, True)

            frame2 = p2.DrawFrame(0., yMin_ratio, 20., yMax_ratio);
            FrameSettings(frame2, "#it{p}_{T} (GeV/#it{c})", "#frac{{{0}}}{{{1}}}".format(energy_other_file, info_decay.get_energy()),0., 20.)
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
        txt = set_txt(0.1,0.9,1.,0.95,"NDC",0.045, 22)
        txt.AddText("{}".format(TString(title)))
        txt.Draw()
        ROOT.SetOwnership(txt,False);
        set_txt_of_plot(info_decay, 0.87, 0.83, 0.87, 0.83, False, 0.025, len(histo_list), False, True)

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(name_plot.Data());
        c1.Close();




