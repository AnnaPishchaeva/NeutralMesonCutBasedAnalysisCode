How to use PhotonMomentumResolution_new_output.

################################################################################
There are two main files in the code:

i) config_run.yml (main config file) which stores in "common" all of the information about what kind of mesons we're using for running (pi0, eta); what kind of data (mc, data); what kind of decay (dalitz, pcm); what kind of system (Pb-Pb); what kind of energy; where the data files are located; where the config files are locates; where the comaprison files are located.

In "common" one can also find the True/False statements. They are important since  they indicate what type of subsystems your initial root file has.
Thus,
"diff_RV0 - True" indicates that a file has subsystems that differ by ranges in the V0 radiuses (but have the same centrality). Since for "diff_RV0" one has the same centrality (usually this centrality is something like 0-80%, a.k.a. all of the centrality range), one needs to also put "all_cent - True" when he/she has "diff_RV0 - True".
"diff_cent - True" indicates that a file has subsystems that differ by centrality (but have the same range in the radius of V0); there is no need to put "all_cent - True" only if she/he wants to include the plots of centrality ranges of the full range. Also when "diff_cent - True", then all of the subsytems have the same range in radius of V0.
However, there can be a case where one analyses root files which have subsystems of different centralities and the root files themselves differ from one another by radius of V0. Now, one wants to compare different radiuses of V0 for different centralities, then "diff_Cent - True" and when one has gone through all of the files, one also puts "compare_RV0 - True" (for more explicit explanation see below "How to use RV0_bin in the config files")

If one doesn't have any subsystems (a.k.a only one subsystem in the file with only one centrality, RV0, etc), then "diff_RV0", "diff_cent", "all_cent", "compare_RV0" must be False

To use the program, one must change the path to where the code lies in "my_dir", and the path for the input data for the analysis in "my_dir_data"

The "output_folder" is the name of the folder where the output data will be stored. No need to put the path there, just the name!

For the input_file_data(mc), you must put not only the name of the file, but also the subfolder if the data lies not directly in the "my_dir_data" but in its subfolder (a.k.a. if the data lies in my_dir_data/subfolder/data.root then one writes "subfolder/data.root")

Also it allows certain programs to run using True/False in "statements".

First, we MUST allow "inv_mass_raw_yield" to run (by placing True) with both "data" and "mc". It will extract inv masses and calculate raw yield for each particle, decay, dataset (mc, data), and centrality. "inv_mass_raw_yield" means the loop over dataset (mc, data)

After that, we MUST run "eff_acc_corr_yield"  to find corrected yields and efficiencies with acceptance(eff_acc). Since for eff_acc we need mc, and we correct the raw yield of data, we use both types of dataset simultaneously.

Then, if we want to compare decays between each other (pcm and dalitz) or mesons between each other (pi0 and eta), we run either "comp_dalitz_pcm" or "comp_eta_pi0".

The line "comparison_to_other_data" is responsible for comparison of this dataset to the dataset of different energies or systems that we acquire externally.

The line "all_cent" is responsible for whether we include the centrality 0-100% in our root, pkl, pdfs files and work with it. True - we include, False - we leave it out

ii) run_new.py is file that reads off the information (with the class Read_config.py) from the config_run.yml and runs the rogram itself.

To run the program python3 run_new.py > test.log

(test.log will store the information about the code running. So that if there is an error, one can easily identify where it happened)
################################################################################

################################################################################
Which classes are used to run the code:

1) Analysis.py is the program where the essential part of the analysis happens. It calls out the class Calculation.py which extracts inv masses, calculates yields. It also compares the data and calls out the classes that write histos in root, pdf, pkl files

2) Calculation.py is the program where the extractions and calculations of inv masses, yields happen. It calls out the classes that directly calculate everything

3) Pair_analyzer.py is the program where the inv_mass with the fir parameters are extracted and raw yield, stat. uncert, and significance are calculated

4) Corrected_yield.py, Efficiency.py,... are programs that calculate corrected yield and efficiency*acceptance

5) Write_pdf.py is the program that writes down the histos into pdf files

6) Write_in_file.py is the program that writes down the histos into root, pkl files

7) Get_histo.py is the program that extracts histos from root, pkl files

8) Plot_(inv_mass, yield, efficiency) are the programs that plots particular histograms on pdf

9) Compare.py is the program that compares datasets or mesons or decays between each other and writes them down

10) Read_config.py is the program that reads the information from the main config file "config_run.yml" and gives it to "run_new.py"
################################################################################

################################################################################
Which config files are also used to run the code:

Besides using the main config file "config_run.yml" that indicates the main points to run the code (which system one uses, measons, dataset, where the dataset lies, etc),
there are also other config files (such as "config_PbPb_5.36TeV_eta_LHC23zzh_full_pcm.yml") that contains information about the subsystems that the initial root file has, also which pT range to use, etc)

How to use RV0_bin in the config files:
If one has "diff_RV0" then what matters is the first and the last element of the array of RV0_bin. They will establish the overall range of RV0 (for instance [4,30,58,90] => the overall range is Rmin = 4, Rmax= 90)
If one has "diff_cent" but also has several root files for different RV0s then one should do the following IN EACH FUNCTIONS (for instance in "inv_mass_raw_yield"):
first, go through each root file, for that in RV0_bin establish the range of RV0 that belongs to the file (for instance, file1.root has RV0 of 30 to 58 => RV0_bin [30, 58])
second, when one has gone through all of the files, one must establish RV0_bin with all of the RV0 ranges (for instance RV0_bin [4,30,58,90] and also put "is_compare_RV0")

The overall description of the config files for particles:

"common" one establishes:
- system and energy (which must be the same as ones in the main config file);
- meson
- period
- plotting_type
- background_fit (for the drawing legend of the inv_mass plot)
- pt_bin
- RV0_bin
- start_pT_bin

"data" or "mc" (depends which type of the data you are looking at):
- name or mc_name (again, depending on the type of the data) of your subsystem
- initial_fit_parameters (initial values to fit the inv mass peak)
- fit_parameters_name (indicates which value of the initial_fit_parameters and fit ranges correspond to which parameter)
- fit_param_limit_min and fit_param_limit_max gives the ranges of the fit values in percentage
- fit_func (gives the fit function of the inv mass peak).
Very IMPORTANT: the number of fit functions inside the array of fit_func can be either 1 ["fit function name"]=> it will mean that all of the pt ranges are going to be fitted with the same function
or it can be the number of len(pt_bin)-1, a.k.a ["name1", "name2", ...,"name(pt-1)"] (so, it you have len(pt_bin)=5, then len(fit_func)=4) because fit function corresponds not to the value of pt but to the range [pt1, pt2].
Then, if the number of fit functions are len(pt_bin) then each pt_range will have its own fit function.
- fit_range (gives the range in which fir function will be fitted).
Very IMPORTANT: the same thing as with the fit_func - you can have either one fit_range. The you MUST write it as [[min, max]] and not [min, max]. Then this range is applied to all of the pt ranges.
or it can be the number of len(pt_bin) -1, a.k.a. [[min1, max1], [min2, max2], ..., [min(pt-1), max(pt-1)]]
- scale_miv_range (gives the range in which inv mass from mix events will be scaled)
Very IMPORTANT: works the same way as fit_range
- yield_range (it gives values yield_min, yield_max which will be used to get the integration range to get the war yield. Integration range will be [mean-yield_min, mean+yield_max], where mean comes from the fit)
Very IMPORTANT: works the same way as fit_range
################################################################################







