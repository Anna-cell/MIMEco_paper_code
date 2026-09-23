#Condition 4
import pandas as pd
import cobra
import mimeco
from mimeco import analysis

import sys
print(sys.prefix)

import os
os.environ["GRB_LICENSE_FILE"] = "/home/lamberta/gurobi/gurobi.lic"


min_medium = pd.read_csv('mimeco_validation/media/min_medium.csv', sep=',', header=0, index_col=0)


NoIle = cobra.io.read_sbml_model('mimeco_validation/models/e_coli_noIle.xml')
NoIle.solver = "gurobi"
NoIle.id = "E.coli_ΔIle"
NoLys = cobra.io.read_sbml_model('mimeco_validation/models/e_coli_noLys.xml')
NoLys.solver = "gurobi"
NoLys.id = "E.coli_ΔLys"

e_coli = cobra.io.read_sbml_model('mimeco_validation/models/iJO1366.xml')
e_coli, constrained_medium_dict1 = mimeco.utils.restrain_medium(e_coli, min_medium, undescribed_metabolites_constraint = "blocked")
solo_growth_e_coli = e_coli.optimize().objective_value
print(solo_growth_e_coli) 

#Check deletion
print(NoIle.reactions.get_by_id('THRD_L').bounds) #(0, 0)
print(NoLys.reactions.get_by_id('THRD_L').bounds) #(1000, 0)
print(NoIle.reactions.get_by_id('DAPDC').bounds) #(1000, 0)
print(NoLys.reactions.get_by_id('DAPDC').bounds) #(0, 0)

min_medium_altered = min_medium.copy()
min_medium_altered.loc["lys__L"] = 0.0001
min_medium_altered.loc["ile__L"] = 0.0001

#Check growth with and without Ile in the medium for e.coli auxotrophic for Ile
with NoIle:
    NoIle, constrained_medium_dict1 = mimeco.utils.restrain_medium(NoIle, min_medium_altered, undescribed_metabolites_constraint = "blocked")
    solo_growth_NoIle = NoIle.optimize().objective_value
    print(solo_growth_NoIle) #0.00034
    NoIle.reactions.get_by_id("EX_ile__L_e").lower_bound = -1
    NoIle.reactions.get_by_id("EX_lys__L_e").lower_bound = -1
    solo_growth_NoIle = NoIle.optimize().objective_value
    print(solo_growth_NoIle) #0.724

#Check growth with and without Lys in the medium for e.coli auxotrophic for Lys
with NoLys:
    NoLys, constrained_medium_dict1 = mimeco.utils.restrain_medium(NoLys, min_medium_altered, undescribed_metabolites_constraint = "blocked")
    solo_growth_NoLys = NoLys.optimize().objective_value
    print(solo_growth_NoLys) #0.00029
    NoLys.reactions.get_by_id("EX_ile__L_e").lower_bound = -1
    NoLys.reactions.get_by_id("EX_lys__L_e").lower_bound = -1
    solo_growth_NoLys = NoLys.optimize().objective_value
    print(solo_growth_NoLys) #0.726

#Check growth of WT e.coli in minimal medium with very few Ile or Lys
e_coli = cobra.io.read_sbml_model('mimeco_validation/models/iJO1366.xml')
e_coli.solver = "gurobi"
e_coli, constrained_medium_dict1 = mimeco.utils.restrain_medium(e_coli, min_medium_altered, undescribed_metabolites_constraint = "blocked",)
solo_growth_coli = e_coli.optimize().objective_value
print(solo_growth_coli) #0.712

print("WT_e_coli")
with e_coli:
    e_coli, constrained_medium_dict1 = mimeco.utils.restrain_medium(e_coli, min_medium_altered, undescribed_metabolites_constraint = "blocked")
    solo_growth_NoIle = e_coli.optimize().objective_value
    print("mm ",solo_growth_NoIle)
    e_coli.reactions.get_by_id("EX_ile__L_e").lower_bound = -1
    e_coli.reactions.get_by_id("EX_lys__L_e").lower_bound = -1
    solo_growth_e_coli = e_coli.optimize().objective_value
    print("mm + AA ",solo_growth_e_coli)

print("WT_e_coli")
int_score, int_type, xy = analysis.interaction_score_and_type(e_coli, e_coli, min_medium_altered, 
                                                          undescribed_metabolites_constraint="blocked", 
                                                          verbose = True, plot=True, retrieve_data = "all")


print(int_score) 
print(int_type) 

# Model auxotrophs together with very few Ile and Lys and see how they interact

print("delta")
int_score, int_type, xy = analysis.interaction_score_and_type(NoLys, NoIle, min_medium_altered, 
                                                          undescribed_metabolites_constraint="blocked", 
                                                          verbose = True, plot=True, retrieve_data = "all")
import matplotlib.pyplot as plt

def pareto_plot(xy, model1_id, model2_id):
    plt.title("Pareto front of "+model1_id+" - "+model2_id+" metabolic interaction")
    plt.xlabel(model1_id+"'s objective value")
    plt.ylabel(model2_id+"'s objective value")
    plt.plot(xy['x'].to_numpy(), xy['y'].to_numpy(), 'dimgray', linestyle="-")
    plt.fill_between(xy['x'].to_numpy(), xy['y'].to_numpy(), color = "lightgray")
    plt.axhline(y = 1, color = '#009e73ff', linestyle = '--', linewidth = 1)
    plt.axvline(x = 1, color = '#009e73ff', linestyle = '--', linewidth = 1)
    plt.savefig("mimeco_validation/results/Pareto_coli.pdf", format="pdf", bbox_inches="tight")
    plt.show()

pareto_plot(xy, "NoLys", "NoIle")

print(int_score) #2283562.321880122
print(int_type) #Limited mutualism


NoIle_biomass_id = "BIOMASS_Ec_iJO1366_WT_53p95M"
NoLys_biomass_id = "BIOMASS_Ec_iJO1366_WT_53p95M"

potential_exchange = mimeco.analysis.crossfed_metabolites(model1 = NoIle, model2 = NoLys, medium = min_medium_altered, undescribed_metabolites_constraint="blocked", 
                                                        solver = "gurobi", model1_biomass_id = NoIle_biomass_id, model2_biomass_id = NoLys_biomass_id, plot = True)

print(potential_exchange)