import pickle
infiles={}
#0
infiles["Au_Elab_1_23_3e"] = ("results_Au_Elab_1_23_3e.pickle.gz","3 events ")

#1
infiles["Au_Elab_1_23_15e"] = ("results_Au_Elab_1_23_15e.pickle.gz","15 events ")

#2
infiles["Au_Elab_1_23_90e"] = ("results_Au_Elab_1_23_90e.pickle.gz","90 events ")

#3
infiles["Au_Elab_1_23_360e"] = ("results_Au_Elab_1_23_360e.pickle.gz","360 events ")

#4
infiles["Au_Elab_1_23_1080e"] = ("results_Au_Elab_1_23_centr.pickle.gz","1080 events ")


with open("num_events_dep.pickle","wb") as outf:
    pickle.dump(infiles,outf)

