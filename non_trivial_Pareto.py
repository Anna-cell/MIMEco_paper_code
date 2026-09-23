#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 17:44:15 2026

@author: lamberta
"""

#Condition 1
import pandas as pd
import cobra
import mimeco
from mimeco import analysis
import matplotlib.pyplot as plt
#cobra.Configuration().solver = "glpk"

import sys
print(sys.prefix)

import os
os.environ["GRB_LICENSE_FILE"] = "/home/lamberta/gurobi/gurobi.lic"

min_medium = pd.read_csv('mimeco_validation/media/min_medium.csv', sep=',', header=0, index_col=0)

M_berkei = cobra.io.read_sbml_model('mimeco_validation/models/iAF692.xml')
M_berkei.solver = "gurobi"
M_berkei.objective = M_berkei.reactions.get_by_id("BIOMASS_Mb_30")
print(M_berkei.optimize())
        
G_meta = cobra.io.read_sbml_model('mimeco_validation/models/iAF987.xml')
G_meta.solver = "gurobi"
#G_meta.objective = "BIOMASS_Gm_GS15_core_79p20M"
G_meta.objective = G_meta.reactions.get_by_id("BIOMASS_Gm_GS15_WT_79p20M")

        
#Weird upper bounds : 
    
for ex in G_meta.exchanges:
    if ex.upper_bound < 0:
        print(ex)
        print(ex.upper_bound)
        ex.upper_bound = 0

"""    
EX_ac_e: ac_e <-- 
-6.84
EX_fe3_e: fe3_e <-- 
-49.21
"""

"""      
min_medium.loc["ac"] = 
min_medium.loc["fe3"] = 
"""
        
G_meta.solver = "gurobi"

int_score, int_type, xy = analysis.interaction_score_and_type(M_berkei, G_meta, min_medium, undescribed_metabolites_constraint="partially_constrained", 
                                                          undescribed_met_lb = -0.1, plot=True, verbose=True, retrieve_data="all")

print(int_score, int_type)


#manual plotting:

def pareto_plot(xy, model1_id, model2_id):
    plt.title("Pareto front of "+model1_id+" - "+model2_id+" metabolic interaction")
    plt.xlabel(model1_id+"'s objective value")
    plt.ylabel(model2_id+"'s objective value")
    plt.plot(xy['x'].to_numpy(), xy['y'].to_numpy(), 'dimgray', linestyle="-")
    plt.fill_between(xy['x'].to_numpy(), xy['y'].to_numpy(), color = "lightgray")
    plt.axhline(y = 1, color = '#009e73ff', linestyle = '--', linewidth = 1)
    plt.axvline(x = 1, color = '#009e73ff', linestyle = '--', linewidth = 1)
    #plt.scatter(5.33, 1.71, color="#00a600ff")
    #plt.scatter(4.91, 1.01, color="chartreuse", zorder=5, alpha=1.0, edgecolors='#00a600ff')
    #plt.plot([5.33, 4.91],[1.71, 1.01], '#00a600ff', linestyle=":")
    plt.savefig("mimeco_validation/results/Pareto.pdf", format="pdf", bbox_inches="tight")
    plt.show()

pareto_plot(xy, "M_berkei", "G_meta")
#M_berkei_obj = M_berkei.reactions.get_by_id("EX_cpd11416_c0").id
#G_meta_obj = G_meta.reactions.get_by_id("EX_cpd11416_c0").id
#potential_exchange = mimeco.analysis.crossfed_metabolites(model1 = A_ory, model2 = L_pla, medium = YPD, undescribed_metabolites_constraint="partially_constrained", 
#                                                        undescribed_met_lb = -0.0001, solver = "gurobi", model1_biomass_id = A_ory_obj, model2_biomass_id = L_pla_obj, plot = False)
