This repository contains the scripts and the configuration files used (in addition to the SMASH code: https://github.com/smash-transport/smash) to produce the results presented in:
https://arxiv.org/abs/2201.05934

## Procedure

We run SMASH (for most results of this work we used commit 3b37a429119b75a141d3c3357a81b5e327099300) using the SMASH .yaml config files in this directory multiple times in parallel on different nodes of a HPC cluster.
The script that manages the job in the cluster takes care to run some simulations and prostprocess the results with the python3 script *analyze_data_ncl.py*.
The results of the postprocessing is a python3 pickled archive. In general, several jobs for the same system configuration are run in parallel and the results are combined at the end with the script (used also during intermediate steps of the job) *combine_data.py*.


After preparing the data, one needs to run once the script:
*create_pickled_list_of_files_to_analize.py*
that produces in output a few files with  pickled python dictionaries that associate the filenames of the output files to the name of the systems and a basic description that will be used in the labels of the plots.
_Please, note that the script must be edited accordingly to the names of your output files before run it._

Then one needs to run the scripts that make the plots and the files in ascii format with the content of the plots.
One can do it in one shot with a script or launching the following commands 
(adapting to the local context the paths and the file names, here we assume that the results are in the subdirectory pickled_data_final):

python3 make_1D_plots_at_center_vs_time_production_version_1.6.0.py pickled_data_final/lh_list.pickle plots_1D_at_center_vs_time_lh
python3 make_1D_plots_at_center_vs_time_production_version_1.5.1.py pickled_data_final/centr_dep_list.pickle plots_1D_at_center_vs_time_Au_centr_dep
python3 make_1D_plots_at_center_vs_time_production_version_1.5.1.py pickled_data_final/enscan_list.pickle plots_1D_at_center_vs_time_Au_enscan

python3 make_1D_plots_percentage_XY_vs_time_production_version_1.5.1.py pickled_data_final/lh_list.pickle.pickle plots_1D_percentage_lh
python3 make_1D_plots_percentage_XY_vs_time_production_version_1.5.1.py pickled_data_final/centr_dep_list.pickle plots_1D_percentage_Au_centr_dep
python3 make_1D_plots_percentage_XY_vs_time_production_version_1.5.1.py pickled_data_final/enscan_list.pickle plots_1D_percentage_Au_enscan
(these last three commands in some cases produce the same plots for single rections)

python3 make_2D_imgs_production_version_1.2.py pickled_data_final/half1.pickle plots_2D_imgs
python3 make_2D_imgs_production_version_1.2.py pickled_data_final/half2.pickle plots_2D_imgs

python3 make_2D_plots_production_version1.3.py pickled_data_final/half1.pickle plots_2D
python3 make_2D_plots_production_version1.3.py pickled_data_final/half2.pickle plots_2D

python3 make_4D_plots_production_version1.3.py pickled_data_final/half1.pickle plots_4D
python3 make_4D_plots_production_version1.3.py pickled_data_final/half2.pickle plots_4D

python3 make_plots_4D_volumes_systems_production_version_1.4.py pickled_data_final/all.pickle tables_integral

For plots 18 and 19 in version 2 (then manually picked and renamed):
python3 make_plots_4D_volumes_systems_production_version_1.4.py pickled_data_final/volume_absolute.pickle volume_absolute 
python3 make_plots_4D_volumes_systems_production_version_1.4.py pickled_data_final/volume_ratio.pickle volume_ratio

After the scripts have finished, one can decide to rename the output files in the directories, for example with something like:
- plots_1D_at_center_vs_time_lh with `for i in *; do mv $i ${i/./_light_ions.}; done`
- plots_1D_at_center_vs_time_Au_centr_dep with `for i in *; do mv $i ${i/./_Au_centrality_dependence.}; done`
- plots_1D_at_center_vs_time_Au_enscan with `for i in *; do mv $i ${i/./_Au_snn_dependence.}; done`
in the case that one wants to put the figures, which otherwise would have the same name, in the same directory.

Most of the scripts that make the plots waste a lot of memory (this is why sometimes they are run in two parts),
but, if needed, they can be adapted so to spare resources.

In the case of the plots made with Paraview, after running SMASH using the proper config files, one needs:
- in the case of event by event simulations, one needs first to average the results with the script average_tmunu_vtk_grid.py,
  e.g. with `python3 average_tmunu_vtk_grid.py data.vtk *landau_tmn*013.vtk`, getting in output the average T^munu and, for convenience, also the pressure anisotropy X and the off-diagonality Y in vtk format
- in the case of simulation with test particles, load the proper landau_tmn file (in this paper the one at timestep 13) and then use the filter calculator with the following formulas:
  `(abs(hadron_tmn_landau11-hadron_tmn_landau22)+abs(hadron_tmn_landau11-hadron_tmn_landau33)+abs(hadron_tmn_landau22-hadron_tmn_landau33))/(hadron_tmn_landau11+hadron_tmn_landau22+hadron_tmn_landau33)` for X
  `3*(abs(hadron_tmn_landau12)+abs(hadron_tmn_landau13)+abs(hadron_tmn_landau23))/(hadron_tmn_landau11+hadron_tmn_landau22+hadron_tmn_landau33)` for Y
The plots in the paper are slices at z=0.

The plots in the appendix about Au+Au collisions at Elab=80 AGeV, b=6fm, to compare with previous UrQMD results, are created with:

python3 make_plots_Elab80b6.py pickled_data_final/results_Au_Elab_80_b6.pickle.gz plots_AuElab80b6
and
python3 make_plots_Elab80b6_comparison.py 0.3 plots_AuElab80b6_comparison

directly from the processed data, in the second case with the names of the data files written inside the script.

The plots about the clumps are produced from data created with the script *analyze_data_wcl.py* (equal to *analyze_data_ncl.py*, except for *search_clumps=True*)
The script to actually make the plot is:

python3 make_plots_clumps_production_version_1.py plots_clumps

in this case the pickled file with the results is hardcoded inside the script, so one needs to edit it to use it on a file different than
pickled_data_final/results_Au_Elab_1_23_wcl.pickle.gz

The plots which compare the time evolution of system volume fraction of the same central Au+Au system at Elab=1.23 AGeV depending on the number of events are made with:
*create_pickled_list_of_files_to_analize_num_events.py* and *make_1D_plots_percentage_XY_vs_time_production_version_1.5.1.py*

The plots which compare the time evolution of the energy density, X_ebe and Y_ebe in the center of the grid with errorbands come from simulations with
SMASH's normal (not lattice) thermodynamic output in one point (see later for the config file), processed with:
*process_evo_at_center_w_stdev.py* and *plot_avg_with_sdev_at_center.py*

The tables of the integrals are made with:
python3 make_table_integrals_production_version_1.4.1.py pickled_data_final/all.pickle table_integrals
with some post editing and numbers round-off for the paper.

## Details about the scripts

Most scripts return the syntax with the list of arguments if run without any or with the wrong number of arguments.

## Details about the config files

The name of the config files should also deliver sufficient information about the target systems.

## Correspondence between figures, config files and pythons scripts:

Figures 1-4:
energy_density_evolution_at_grid_center_Au_snn_dependence.pdf
energy_density_evolution_at_grid_center_Au_snn_dependence_scaled.pdf
X_ebe_evolution_at_grid_center_avg_single_events_Au_snn_dependence.pdf
Y_ebe_evolution_at_grid_center_avg_single_events_Au_snn_dependence.pdf
X_evolution_at_grid_center_avg_Tmunu_Au_snn_dependence.pdf
Y_evolution_at_grid_center_avg_Tmunu_Au_snn_dependence.pdf
config files:
config_Au_Elab_1_23_centr_ncl.yaml
config_Au_Elab_3_4.yaml
config_Au_Elab8.yaml
config_Au_Elab12.yaml
config_Au_Ecm77.yaml
python scripts:
analyze_data_ncl.py, combine_data.py, create_pickled_list_of_files_to_analize.py, make_1D_plots_at_center_vs_time_production_version_1.5.1.py
(make_1D_plots_at_center_vs_time_production_version_1.6.0.py for "_scaled_" plots)

Note: figure energy_density_evolution_at_grid_center_Au_snn_dependence_scaled.pdf added in version 2, manually renamed)

Figures 5-6:
energy_density_evolution_at_grid_center_Au_centrality_dependence.pdf
X_evolution_at_grid_center_avg_Tmunu_Au_centrality_dependence.pdf
Y_evolution_at_grid_center_avg_Tmunu_Au_centrality_dependence.pdf
config files:
config_Au_Elab_1_23_centr_ncl.yaml
config_Au_Elab_1_23_semic.yaml
config_Au_Elab_1_23_per.yaml
python scripts:
analyze_data_ncl.py, combine_data.py, create_pickled_list_of_files_to_analize.py, make_1D_plots_at_center_vs_time_production_version_1.5.1.py

Figures 7-8:
energy_density_evolution_at_grid_center_light_ions.pdf
X_evolution_at_grid_center_avg_Tmunu_light_ions.pdf
Y_evolution_at_grid_center_avg_Tmunu_light_ions.pdf
config files:
config_Au_Elab_1_23_centr_ncl.yaml
config_Ag_Elab_1_58.yaml
config_ArKCl_1_756.yaml
config_C_Elab2.yaml
python scripts:
analyze_data_ncl.py, combine_data.py, create_pickled_list_of_files_to_analize.py, make_1D_plots_at_center_vs_time_production_version_1.5.1.py

Figures 9-10:
energy_density_evolution_at_grid_center_w_errorbands.pdf
X_ebe_evolution_at_grid_center_w_errorbands.pdf
Y_ebe_evolution_at_grid_center_w_errorbands.pdf
config files:
config_Au_Elab123_central_point.yaml
config_Ag_Elab158_central_point.yaml
config_Au_Ecm77_central_point.yaml
python scripts:
process_evo_at_center_w_stdev.py, plot_avg_with_sdev_at_center.py

Figure 11:
Au_Elab_1_23_b_0_3_3_energy_density_t_14.00.png
Au_Elab_1_23_b_0_3_3_X_t_14.00.png
Au_Elab_1_23_b_0_3_3_Y_t_14.00.png
config files:
config_Au_Elab_1_23_centr_ncl.yaml
python scripts:
analyze_data_ncl.py, combine_data.py, create_pickled_list_of_files_to_analize.py, make_2D_imgs_production_version_1.2.py

Figure 12:
Au_Elab_1_23_b_0_3_3_edens_vs_X_ebe_2Dhist_t_14.00fm.png
Au_Elab_1_23_b_0_3_3_edens_vs_X_2Dhist_avg_t_14.00fm.png
Au_Elab_1_23_b_0_3_3_edens_vs_Y_ebe_2Dhist_t_14.00fm.png
Au_Elab_1_23_b_0_3_3_edens_vs_Y_2Dhist_avg_t_14.00fm.png
Au_Elab_1_23_b_0_3_3_Y_ebe_vs_X_ebe_2Dhist_t_14.00fm.png
Au_Elab_1_23_b_0_3_3_Y_vs_X_2Dhist_avg_t_14.00fm.png
config files:
config_Au_Elab_1_23_centr_ncl.yaml
python scripts:
analyze_data_ncl.py, combine_data.py, create_pickled_list_of_files_to_analize.py, make_2D_plots_production_version_1.3.py

Figures 13-16:
percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_system_comparison_Xlim_0.30_Ylim_0.30_edens_001_Au_enscan.pdf
percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_system_comparison_Xlim_0.30_Ylim_0.30_edens_500_Au_enscan.pdf
percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_system_comparison_Xlim_0.30_Ylim_0.30_edens_001_Au_centr_dep.pdf
percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_system_comparison_Xlim_0.30_Ylim_0.30_edens_001_lh.pdf
config files:
config_Au_Elab_1_23_centr_ncl.yaml
config_Au_Elab_3_4.yaml
config_Au_Elab8.yaml
config_Au_Elab12.yaml
config_Au_Ecm77.yaml
config_Ag_Elab_1_58.yaml
config_ArKCl_1_756.yaml
config_C_Elab2.yaml
config_Au_Elab_1_23_semic.yaml
config_Au_Elab_1_23_per.yaml
python scripts:
analyze_data_ncl.py, combine_data.py, create_pickled_list_of_files_to_analize.py, make_1D_plots_percentage_XY_vs_time_production_version_1.5.1.py

Figure 17:
percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_system_comparison_Xlim_0.30_Ylim_0.30_edens_001_events_dep.pdf
config file:
config_Au_Elab_1_23_centr_ncl.yaml
python scripts:
analyze_data_ncl.py, combine_data.py, create_pickled_list_of_files_to_analize_num_events.py, make_1D_plots_percentage_XY_vs_time_production_version_1.5.1.py

Figures 18-19:
4Dvolume_system_comparison_Xlim_0.30_Ylim_0.30_edens_001.pdf
4Dvolume_ratio_of_integrals_system_comparison_Xlim_0.30_Ylim_0.30_edens_001.pdf
config files:
config_Au_Elab_1_23_centr_ncl.yaml
config_Au_Elab_3_4.yaml
config_Au_Elab8.yaml
config_Au_Elab12.yaml
config_Au_Ecm77.yaml
config_Ag_Elab_1_58.yaml
config_ArKCl_1_756.yaml
config_C_Elab2.yaml
config_Au_Elab_1_23_semic.yaml
config_Au_Elab_1_23_per.yaml
python scripts:
analyze_data_ncl.py, combine_data.py, create_pickled_list_of_files_to_analize.py, make_plots_4D_volumes_systems_production_version_1.4.py

Note: in version 2 Figures 17 and 19 have become 18 and 19 respectively, and are produced with dedicated pickled lists of files
volume_absolute.pickle and volume_ratio.pickle.

Figures 20-22:
clumps__Xlim_0.30_Ylim_0.30.pdf
fraction_of_single_cell_clumps__Xlim_0.30_Ylim_0.30.pdf
compactness_of_clumps__Xlim_0.30_Ylim_0.30.pdf
config files:
config_Au_Elab_1_23_wclumps.yaml
python scripts:
analyze_data_wcl.py, combine_data.py, make_plots_clumps_production_version_1.py

Figure 23:
X_slice_1080_events_no_spectators.jpg
Y_slice_1080_events_no_spectators.jpg
config files:
config_Au_Elab1_23_no_spect_vtk_out.yaml (example file, might not be exactly what has been used - different random seed or grid size)
python scripts:
average_tmunu_vtk_grid.py
plots made with Paraview

Figure 24:
X_slice_1080_testparticles_no_spectators.jpg
Y_slice_1080_testparticles_no_spectators.jpg
config file:
config_Au_Elab1_23_no_spect_testparticles_vtk_out.yaml (example file, might not be exactly what has been used - different random seed or grid size)
python scripts:
average_tmunu_vtk_grid.py
plots made with Paraview

Figure 25:
slice_1080_events_with_spectators.jpg
Y_slice_1080_events_with_spectators.jpg
config file:
config_Au_Elab1_23_w_spect_vtk_out.yaml
python scripts:
average_tmunu_vtk_grid.py
plots made with Paraview

Figures 26:
X_evolution_Au_Elab80_b6.pdf
config file:
config_Au_Elab_80_b6.yaml
python scripts:
analyze_data_ncl.py, combine_data.py, make_plots_Elab80b6.py

Figures 27:
X_lt_03_XZ_plane_percentage_Au_Elab80_b6.pdf
config file:
config_Au_Elab_80_b6.yaml
python scripts:
analyze_data_ncl.py, combine_data.py, make_plots_Elab80b6_comparison.py

## Details about the software versions
compiler gcc 4.8.5
Pythia (SMASH dependency) 8.303
Anaconda distribution 21.05 with Python 3.8.8, Numpy 1.20.1, Matplotlib 3.3.4
Paraview 5.9.1
