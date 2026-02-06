# This code was written by Alica Enderich (Febuary, 2024)

import ROOT
from ROOT import TH1D, TH2D, TH3D, TGraph, TGraphErrors, TLatex, TPaveText, TLegend, TCanvas, TPad, TString, TLine
from ROOT import gROOT, gSystem, gStyle, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kFullCircle
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle, kTRUE

def FrameSettings(pad1, XTitle, YTitle, xMin, xMax, is_one_pad = False):

    if is_one_pad:
        pad1.GetYaxis().SetTitleOffset(1.1)
        pad1.GetYaxis().SetLabelSize(0.05);
        pad1.GetYaxis().SetTitleSize(0.05)
        pad1.GetXaxis().SetLabelSize(0.05);
        pad1.GetXaxis().SetTitleSize(0.05);
    else:
        pad1.GetYaxis().SetTitleOffset(1.0);
        pad1.GetYaxis().SetLabelSize(0.06);
        pad1.GetYaxis().SetTitleSize(0.06)
        pad1.GetXaxis().SetLabelSize(0.06);
        pad1.GetXaxis().SetTitleSize(0.06);
    pad1.GetYaxis().SetTitle(YTitle);
    pad1.GetYaxis().SetDecimals();

    pad1.GetXaxis().SetRangeUser(xMin, xMax);
    pad1.GetXaxis().SetTitleOffset(1.0);
    pad1.GetXaxis().SetTitle(XTitle);
    pad1.GetXaxis().SetNdivisions(507, True);
    pad1.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)

def FrameSettingsCombined(pad1, XTitle, YTitle, xMin, xMax):

    pad1.GetYaxis().SetLabelSize(0.06);
    pad1.GetYaxis().SetTitleSize(0.06);
    pad1.GetYaxis().SetTitleOffset(1.0);
    pad1.GetYaxis().SetTitle(YTitle);
    pad1.GetYaxis().SetDecimals();

    pad1.GetXaxis().SetRangeUser(xMin, xMax);
    pad1.GetXaxis().SetLabelSize(0.06);
    pad1.GetXaxis().SetTitleSize(0.06);
    pad1.GetXaxis().SetTitleOffset(1.0);
    pad1.GetXaxis().SetTitle(XTitle);
    pad1.GetXaxis().SetNdivisions(507, True);
    pad1.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)

def DrawHisto(histo1, markerColor, drawsettings, markerSize = 1.8, markerstyle = kFullCircle):
    histo1.SetMarkerStyle(markerstyle);
    histo1.SetMarkerColor(markerColor);
    histo1.SetMarkerSize(markerSize);
    histo1.SetLineColor(markerColor);
    histo1.SetLineWidth(1);
    histo1.SetFillColor(markerColor);
    histo1.SetFillStyle(0);
    histo1.Draw(drawsettings); #DrawCopy

def DrawHistoCombined(histo1, markerColor, drawsettings, markerSize = 2.2, markerstyle = kFullCircle):
    histo1.SetMarkerStyle(markerstyle);
    histo1.SetMarkerColor(markerColor);
    histo1.SetMarkerSize(markerSize);
    histo1.SetLineColorAlpha(markerColor, 0.6);
    histo1.SetLineWidth(1);
    histo1.SetFillColor(markerColor);
    histo1.SetFillStyle(0);
    histo1.DrawCopy(drawsettings);

def SetTitle(titlePt):
    TitlePlot = TPaveText(0.34, 0.92, 0.64, 0.98, "NDC")
    TitlePlot.AddText("{}".format(titlePt))
    TitlePlot.SetTextColor(1);
    TitlePlot.SetTextSize(0.08);
    TitlePlot.SetFillStyle(0)
    TitlePlot.SetBorderSize(0)
    TitlePlot.Draw();
    ROOT.SetOwnership(TitlePlot, False);

def set_canvas(w, h, tick1, tick2):
    canvas = TCanvas("c1", "", w, h)
    canvas.SetTicks(tick1,tick2);
    canvas.SetTickx();
    canvas.SetTicky();
    canvas.SetLogy(0);
    canvas.SetMargin(0, 0, 0.15, 0);
    canvas.SetFillColor(0);
    return canvas

def set_decay_name(meson, decay):
    if meson == "pi0":
        if decay == "PCM":
            decay_name = "#pi^{0} #rightarrow #gamma#gamma"
        elif decay == "Dalitz":
            decay_name = "#pi^{0} #rightarrow e^{+}e^{-}#gamma"
        else:
            raise ValueError("Not valid decay")
    if meson == "eta":
        if decay == "PCM":
            decay_name = "#eta #rightarrow #gamma#gamma"
        elif decay == "Dalitz":
            decay_name = "#eta #rightarrow e^{+}e^{-}#gamma"
        else:
            raise ValueError("Not valid decay")
    return decay_name


def set_pad(x_low, y_low, x_up, y_up, color, number_columns = 1, number_rows = 1, ispad_legend = False):
    if ispad_legend == True:
        width_legend = 1./number_columns ;
        height_legend = 1./number_rows;
        if number_columns > 7:
            width_legend = 2./number_columns;
        x_low = 1-width_legend
        y_low = 1-height_legend-0.11
    pad = TPad("p1", "", x_low, y_low, x_up, y_up, color)
    pad.SetPad(x_low, y_low, x_up, y_up)
    pad.SetFillColor(0);
    pad.GetFrame().SetFillColor(0);
    pad.SetMargin(0, 0, 0.1, 0);
    pad.SetTickx();
    pad.SetTicky();
    if (ispad_legend == False):
        pad.Divide(number_columns, number_rows, 0.0, 0.0)
    return pad

def set_line(point1_x, point1_y, point2_x, point2_y, color, style, width):
    line = TLine(point1_x,point1_y,point2_x,point2_y);
    line.SetLineColor(color);
    line.SetLineStyle(style);
    line.SetLineWidth(width);
    ROOT.SetOwnership(line,False);
    return line

def set_txt(x1,y1,x2,y2,option,size, align):
    txt = TPaveText(x1,y1,x2,y2,option);
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(align);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(size);
    return txt

def set_txt_system(energy, system, decay, cent1, cent2):
    txt = set_txt(0.875,0.88,0.875,0.88,"NDC",0.02, 33)
    txt.AddText("this thesis")
    txt.Draw()
    ROOT.SetOwnership(txt,False)

    txt2 = set_txt(0.875,0.85,0.875,0.85,"NDC",0.02, 33)
    txt2.AddText("centrality {0}%-{1}%".format(cent1, cent2))
    txt2.Draw()
    ROOT.SetOwnership(txt2,False)

    txt3 = set_txt(0.875,0.82,0.875,0.82,"NDC",0.02, 33)
    txt3.AddText("{0}, #sqrt{{s}} = {1}".format(system, energy))
    txt3.Draw()
    ROOT.SetOwnership(txt3,False)

    txt2 = set_txt(0.875,0.79,0.875,0.79,"NDC",0.02, 33)
    txt2.AddText(decay)
    txt2.Draw()
    ROOT.SetOwnership(txt2,False)

def set_txt_of_plot(info_decay, x1, y1, x2, y2, both_mc_and_data = False, text_size = 0.02, is_stat_uncert = False):

    txt = set_txt(x1,y1,x2,y2,"NDC", text_size, 33)
    txt.AddText("this thesis")
    txt.Draw()
    ROOT.SetOwnership(txt,False)
    y_diff = 0.03

    if info_decay.get_plotting_sys() != "diff_cent_and_diff_occupancies" and info_decay.get_plotting_sys() != "occupancies":
        occup1, occup2 = info_decay.get_occupancy()
        if occup1 == "0" and occup2 == "10k":
            occup_range_latex = "0 < FT0Occ < 10^{4}"
        elif occup1 == "10k" and occup2 == "30k":
            occup_range_latex = "10^{4} < FT0Occ < 3x10^{4}"
        elif occup1 == "30k":
            occup_range_latex = "FT0Occ > 3x10^{4}"
        else:
            occup_range_latex = "FT0Occ: all"
        y1, y2 = y1-y_diff, y2-y_diff
        txt2 = set_txt(x1,y1,x2,y2,"NDC",text_size, 33)
        txt2.AddText(occup_range_latex)
        txt2.Draw()
        ROOT.SetOwnership(txt2,False)


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
    y1, y2 = y1-y_diff, y2-y_diff
    txt2 = set_txt(x1,y1,x2,y2,"NDC",text_size, 33)
    txt2.AddText(system_energy)
    txt2.Draw()
    ROOT.SetOwnership(txt2,False)

    if not both_mc_and_data:
        y1, y2 = y1-y_diff, y2-y_diff
        txt3 = set_txt(x1,y1,x2,y2,"NDC",text_size, 33)
        txt3.AddText("{0}, {1}".format(info_decay.get_period(), info_decay.get_typ()))
        txt3.Draw()
        ROOT.SetOwnership(txt3,False)

    y1, y2 = y1-y_diff, y2-y_diff
    txt2 = set_txt(x1,y1,x2,y2,"NDC",text_size, 33)
    txt2.AddText("#gamma rec. with PCM")
    txt2.Draw()
    ROOT.SetOwnership(txt2,False)

    decay_name = set_decay_name(info_decay.get_meson(), info_decay.get_decay())
    y1, y2 = y1-y_diff, y2-y_diff
    txt2 = set_txt(x1,y1,x2,y2,"NDC",text_size, 33)
    txt2.AddText(decay_name)
    txt2.Draw()
    ROOT.SetOwnership(txt2,False)

    if (info_decay.get_plotting_sys() == "centralities") or (is_stat_uncert and both_mc_and_data):
        RV0_min = info_decay.get_V0_radius_min()
        RV0_max = info_decay.get_V0_radius_max()
        y1, y2 = y1-y_diff, y2-y_diff
        txt2 = set_txt(x1,y1,x2,y2,"NDC",text_size, 33)
        txt2.AddText("{0} < #it{{R}}_{{V^{{0}}}} < {1} cm".format(int(RV0_min), int(RV0_max)))
        txt2.Draw()
        ROOT.SetOwnership(txt2,False)
    if ((info_decay.get_plotting_sys() != "centralities") or (is_stat_uncert and both_mc_and_data)) and info_decay.get_system() != "pp":
        cent1, cent2 = info_decay.get_cent()
        y1, y2 = y1-y_diff, y2-y_diff
        txt2 = set_txt(x1,y1,x2,y2,"NDC",text_size, 33)
        txt2.AddText("{0}% - {1}% centrality".format(cent1, cent2))
        txt2.Draw()
        ROOT.SetOwnership(txt2,False)

    if not both_mc_and_data:
        nev = info_decay.get_Nev_cent()
        Nev = "Number of events: {}".format(f"{nev:.2e}")
        y1, y2 = y1-y_diff, y2-y_diff
        txt2 = set_txt(x1,y1,x2,y2,"NDC",text_size, 33)
        txt2.AddText(Nev)
        txt2.Draw()
        ROOT.SetOwnership(txt2,False)

def set_txt_pad_legend(pad_legend, number_columns, number_rows, info_decay, n_pix):
    pad_legend.cd()
    text_height     = 0.08;
    txt_alice = "";
    if TString(info_decay.plotting_type).CompareTo("wip")==0:
        txt_alice       = "ALICE work in progress";
    elif TString(info_decay.plotting_type).CompareTo("thesis")==0:
        txt_alice       = "this thesis"; #only this thesis
    elif TString(info_decay.plotting_type).CompareTo("performance")==0:
        txt_alice       = "ALICE performance";
    else:
        txt_alice       = "ALICE";

    if  pad_legend.XtoPixel(pad_legend.GetX2()) < pad_legend.YtoPixel(pad_legend.GetY1()):
        txt_height          = n_pix/pad_legend.XtoPixel(pad_legend.GetX2()) ;
    else:
        txt_height          = n_pix/pad_legend.YtoPixel(pad_legend.GetY1());

    print("txt_height", txt_height)
    diff_txt = txt_height*1.50; #1.50

    return txt_alice, txt_height, diff_txt

def set_labels(txt, start_txt_x, start_txt_y, txt_height, diff_txt, index1, index2, index3, index4):

    label = TPaveText(start_txt_x, start_txt_y+index1*diff_txt, start_txt_x + index2*txt_height, start_txt_y + index3*diff_txt+index4, "NDC")
    label.AddText(txt)
    SetStyleTLatex(label, txt_height*3, 1, 1, 42, True)

    return label

def set_legend(positionX, positionY, positionXRight, positionYUp, textSize, columns, header, textFont, margin):
    legend = TLegend(positionX,positionY,positionXRight,positionYUp);
    legend.SetNColumns(columns);
    legend.SetLineColor(0);
    legend.SetLineWidth(0);
    legend.SetFillColor(0);
    legend.SetFillStyle(0);
    legend.SetLineStyle(0);
    legend.SetBorderSize(0);
    legend.SetTextFont(textFont);
    legend.SetTextSize(textSize);
    if margin != 0:
        legend.SetMargin(margin);
    if header.CompareTo("")!= 0:
        legend.SetHeader(header);
    return legend;

def set_axis(axis, label, tsize, tfont, offs, lsize, lfont):
    axis.SetTitle(label)
    axis.SetTitleSize(tsize)
    axis.SetTitleFont(tfont)
    axis.SetTitleOffset(offs)
    axis.SetLabelSize(lsize)
    axis.SetLabelFont(lfont)

def set_frame_basics(frame, x_title, y_title, x1, x2):
    frame.GetXaxis().SetTitle(x_title)
    frame.GetYaxis().SetTitle(y_title)
    frame.GetXaxis().SetRangeUser(x1, x2);
    frame.GetYaxis().SetDecimals();
    frame.GetXaxis().SetNdivisions(507, True);
    frame.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)


def set_frame(frame, st_x, st_y, ot_x, ot_y, sl_x, sl_y, ol_x, ol_y):
    frame.GetXaxis().SetTitleSize(st_x);
    frame.GetYaxis().SetTitleSize(st_y);
    frame.GetXaxis().SetTitleOffset(ot_x);
    frame.GetYaxis().SetTitleOffset(ot_y);
    frame.GetXaxis().SetLabelSize(sl_x);
    frame.GetYaxis().SetLabelSize(sl_y);
    frame.GetXaxis().SetLabelOffset(ol_x);
    frame.GetYaxis().SetLabelOffset(ol_y);

def set_histo(histo, st_x, st_y, sl_x, sl_y):
    histo.GetXaxis().SetTitleSize(st_x);
    histo.GetYaxis().SetTitleSize(st_y);
    histo.GetXaxis().SetLabelSize(sl_x);
    histo.GetYaxis().SetLabelSize(sl_y);
    histo.GetXaxis().SetNdivisions(507, True)
    histo.GetYaxis().SetDecimals()

def SetStyleTLatex(text, textSize, lineWidth, textColor = 1, textFont = 42, kNDC = True, align = 11):
    text.SetTextFont(textFont);
    text.SetTextColor(textColor);
    text.SetTextSize(textSize);
    text.SetLineWidth(lineWidth);
    text.SetTextAlign(align);
    text.SetFillStyle(0)
    text.SetBorderSize(0)