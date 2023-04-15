import pickle
infiles={}
#0
infiles["C_Elab_2_b_0_1_2"] = ("results_C_Elab_2.pickle.gz" ,"C+C, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=2AGeV, b=0-1.2fm ")

#1
infiles["ArKCl_Elab_1_756_b_0_1_784"] = ("results_ArKCl_Elab_1_756.pickle.gz" ,"Ar+KCl, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=1.756AGeV, b=0-1.784fm ")

#2
infiles["Ag_Elab_1_58_b_0_2_44"] = ("results_Ag_Elab_1_58.pickle.gz" ,"Ag+Ag, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=1.58AGeV, b=0-2.44fm ")

#3
infiles["Au_Elab_1_23_b_0_3_3"] = ("results_Au_Elab_1_23_centr.pickle.gz","Au+Au, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=1.23AGeV, b=0-3.3fm ")

#4
infiles["Au_Elab_1_23_b_3_6"] = ("results_Au_Elab_1_23_semic.pickle.gz","Au+Au, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=1.23AGeV, b=3.3-6.6fm ")

#5
infiles["Au_Elab_1_23_b_6_10"] = ("results_Au_Elab_1_23_per.pickle.gz","Au+Au, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=1.23AGeV, b=6.6-10.4fm ")

#6
infiles["Au_Elab_3_4_b_0_3_3"] = ("results_Au_Elab_3_4.pickle.gz","Au+Au, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=3.4AGeV, b=0-3.3fm ")

#7
infiles["Au_Elab_8_b_0_3_3"] = ("results_Au_Elab_8.pickle.gz","Au+Au, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=8AGeV, b=0-3.3fm ")

#8
infiles["Au_Elab_12_b_0_3_3"] = ("results_Au_Elab_12.pickle.gz","Au+Au, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=12AGeV, b=0-3.3fm ")

#9
infiles["Au_Ecm_7_7_b_0_3_3"] = ("results_Au_Ecm_7_7.pickle.gz","Au+Au, "+r"$\sqrt{s_{\mathrm{NN}}}$"+"=7.7GeV, b=0-3.3fm ")

# all
with open("all.pickle","wb") as outf:
    pickle.dump(infiles,outf)

# light ions dictionary
keys_to_use = ["C_Elab_2_b_0_1_2","ArKCl_Elab_1_756_b_0_1_784","Ag_Elab_1_58_b_0_2_44","Au_Elab_1_23_b_0_3_3"]
selected = { key : infiles[key] for key in keys_to_use }
with open("lh_list.pickle","wb") as outf:
    pickle.dump(selected,outf)

# Au+Au energy scan dictionary
keys_to_use = ["Au_Elab_1_23_b_0_3_3","Au_Elab_3_4_b_0_3_3","Au_Elab_8_b_0_3_3","Au_Elab_12_b_0_3_3","Au_Ecm_7_7_b_0_3_3"]
selected = { key : infiles[key] for key in keys_to_use }
with open("enscan_list.pickle","wb") as outf:
    pickle.dump(selected,outf)

# Au+Au centrality dependence dictionary
keys_to_use = ["Au_Elab_1_23_b_0_3_3","Au_Elab_1_23_b_3_6","Au_Elab_1_23_b_6_10"]
selected = { key : infiles[key] for key in keys_to_use }
with open("centr_dep_list.pickle","wb") as outf:
    pickle.dump(selected,outf)

# half list dictionary 1
keys_to_use = ["C_Elab_2_b_0_1_2","ArKCl_Elab_1_756_b_0_1_784","Ag_Elab_1_58_b_0_2_44","Au_Elab_1_23_b_0_3_3","Au_Elab_1_23_b_3_6"]
selected = { key : infiles[key] for key in keys_to_use }
with open("half1.pickle","wb") as outf:
    pickle.dump(selected,outf)

# half list dictionary 2
keys_to_use = ["Au_Elab_1_23_b_6_10","Au_Elab_3_4_b_0_3_3","Au_Elab_8_b_0_3_3","Au_Elab_12_b_0_3_3","Au_Ecm_7_7_b_0_3_3"]
selected = { key : infiles[key] for key in keys_to_use }
with open("half2.pickle","wb") as outf:
    pickle.dump(selected,outf)

# fig 17
keys_to_use = ["C_Elab_2_b_0_1_2","ArKCl_Elab_1_756_b_0_1_784","Ag_Elab_1_58_b_0_2_44","Au_Elab_1_23_b_6_10","Au_Elab_1_23_b_3_6","Au_Elab_1_23_b_0_3_3","Au_Elab_3_4_b_0_3_3","Au_Elab_8_b_0_3_3","Au_Elab_12_b_0_3_3","Au_Ecm_7_7_b_0_3_3"]
selected = { key : infiles[key] for key in keys_to_use }
with open("volume_absolute.pickle","wb") as outf:
    pickle.dump(selected,outf)

# fig 18
keys_to_use = ["C_Elab_2_b_0_1_2","ArKCl_Elab_1_756_b_0_1_784","Ag_Elab_1_58_b_0_2_44","Au_Elab_3_4_b_0_3_3","Au_Elab_8_b_0_3_3","Au_Elab_12_b_0_3_3","Au_Ecm_7_7_b_0_3_3","Au_Elab_1_23_b_0_3_3","Au_Elab_1_23_b_3_6","Au_Elab_1_23_b_6_10"]
selected = { key : infiles[key] for key in keys_to_use }
with open("volume_ratio.pickle","wb") as outf:
    pickle.dump(selected,outf)
