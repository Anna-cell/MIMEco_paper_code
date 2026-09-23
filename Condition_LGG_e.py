#Replicate LGG - enterocyte interaction analysis from 2024, with MIMEco

import pandas as pd
import cobra
import mimeco
from mimeco import analysis

import sys
print(sys.prefix)

import os
os.environ["GRB_LICENSE_FILE"] = "/home/lamberta/gurobi/gurobi.lic"

#Import same enterocyte model
enterocyte = cobra.io.read_sbml_model("mimeco_validation/models/enterocyte_ori.xml")
enterocyte.objective = enterocyte.reactions.get_by_id('biomass_reactionIEC01b')
enterocyte.solver = 'gurobi'

#Import same LGG model
LGG = cobra.io.read_sbml_model("mimeco_validation/models/Lactobacillus_rhamnosus_GG_GG_ATCC_53103.xml")
LGG.solver = "gurobi"
bacteria_id = 'LGG'

#Import same nutritional constraint (Western diet)
WD = pd.read_csv("/home/lamberta/mimeco_validation/media/Western_diet_BiGG.csv", index_col = 0)

#Translate into the expected format (pandas dataframe with two columns: met_id and Influx)


int_score, int_type = mimeco.analysis.enterocyte_interaction_score_and_type(model = LGG, solver = "gurobi", medium=WD, 
                                                                            undescribed_metabolites_constraint="as_is", plot=True)


potential_crossfeeding = mimeco.analysis.enterocyte_crossfed_metabolites(model = LGG, medium = WD, undescribed_metabolites_constraint = "as_is",
                                                              solver = "gurobi", model_biomass_id = 'Growth', sample_size = 1000)

potential_crossfeeding.to_csv("mimeco_validation/results/exchanges_LGG_enterocyte.csv", index=False)

import matplotlib.pyplot as plt

def pareto_plot(xy, model1_id, model2_id):
    plt.title("Pareto front of "+model1_id+" - "+model2_id+" metabolic interaction")
    plt.xlabel(model1_id+"'s objective value")
    plt.ylabel(model2_id+"'s objective value")
    plt.plot(xy['x'].to_numpy(), xy['y'].to_numpy(), 'dimgray', linestyle="-")
    plt.fill_between(xy['x'].to_numpy(), xy['y'].to_numpy(), color = "lightgray")
    plt.axhline(y = 1, color = '#009e73ff', linestyle = '--', linewidth = 1)
    plt.axvline(x = 1, color = '#009e73ff', linestyle = '--', linewidth = 1)
    plt.scatter(5.33, 1.71, color="#00a600ff")
    plt.scatter(4.91, 1.01, color="chartreuse", zorder=5, alpha=1.0, edgecolors='#00a600ff')
    plt.plot([5.33, 4.91],[1.71, 1.01], '#00a600ff', linestyle=":")
    plt.savefig("mimeco_validation/results/Pareto.pdf", format="pdf", bbox_inches="tight")
    plt.show()

#enterocyte = cobra.io.read_sbml_model("enterocyte_BiGG.xml")