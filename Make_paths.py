# This code was written by Alica Enderich (Febuary, 2024)

import os, sys, shutil

def make_paths(folder):
    subfolders = ["inv_mass_analysis", "inv_yield_calculation", "Comparison"]
    subsubfolders1 = ["Dalitz", "PCM"]
    subsubfolders2 = ["same_mixed_inv_mass", "fitted_inv_mass", "parameters", "raw_yield"]
    subsubfolders3 = ["eff_acc", "corr_yield"]
    os.makedirs(folder, exist_ok=True)
    for subfolder in subfolders:
        os.makedirs(os.path.join(folder, subfolder), exist_ok=True)

    for subsubfolder in subsubfolders1:
        subsubfolder_path = os.path.join(folder, "inv_mass_analysis", subsubfolder)
        os.makedirs(subsubfolder_path, exist_ok=True)
        for subsubsubfolder in subsubfolders2:
            subsubsubfolder_path = os.path.join(folder, "inv_mass_analysis", subsubfolder, subsubsubfolder)
            os.makedirs(subsubsubfolder_path, exist_ok=True)
        
    for subsubfolder in subsubfolders1:
        subsubfolder_path_2 = os.path.join(folder, "inv_yield_calculation", subsubfolder)
        os.makedirs(subsubfolder_path_2, exist_ok=True)
        for subsubsubfolder in subsubfolders3:
            subsubsubfolder_path = os.path.join(folder, "inv_yield_calculation", subsubfolder, subsubsubfolder)
            os.makedirs(subsubsubfolder_path, exist_ok=True)

    os.makedirs("pkl", exist_ok=True)


 