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

YPD = pd.read_csv('mimeco_validation/media/YPD_media_ModelSEED.csv', sep=',', header=0, index_col=0)


A_ory = cobra.io.read_sbml_model('mimeco_validation/models/A_ory_gapseq.xml')
A_ory.solver = "gurobi"
L_pla = cobra.io.read_sbml_model('mimeco_validation/models/L_pla_gapseq.xml')
L_pla.solver = "gurobi"

int_score, int_type, xy = analysis.interaction_score_and_type(A_ory, L_pla, YPD, undescribed_metabolites_constraint="partially_constrained", 
                                                          undescribed_met_lb = -0.0001, plot=True, verbose=True, retrieve_data="all")

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
    plt.scatter(5.33, 1.71, color="#00a600ff")
    plt.scatter(4.91, 1.01, color="chartreuse", zorder=5, alpha=1.0, edgecolors='#00a600ff')
    plt.plot([5.33, 4.91],[1.71, 1.01], '#00a600ff', linestyle=":")
    plt.savefig("mimeco_validation/results/Pareto.pdf", format="pdf", bbox_inches="tight")
    plt.show()

pareto_plot(xy, "A_ory", "L_pla")
A_ory_obj = A_ory.reactions.get_by_id("EX_cpd11416_c0").id
L_pla_obj = L_pla.reactions.get_by_id("EX_cpd11416_c0").id
potential_exchange = mimeco.analysis.crossfed_metabolites(model1 = A_ory, model2 = L_pla, medium = YPD, undescribed_metabolites_constraint="partially_constrained", 
                                                        undescribed_met_lb = -0.0001, solver = "gurobi", model1_biomass_id = A_ory_obj, model2_biomass_id = L_pla_obj, plot = False)

print(potential_exchange)

import pandas as pd
from shapely.geometry import Point, LineString

#Convert DataFrame points into a continuous LineString
line = LineString(xy.to_numpy())

#Define target point
target_point = Point(5.33, 1.71)

# 4. Project the target point onto the line to find the closest point
# line.project() finds the distance ALONG the line to the nearest point.
# line.interpolate() converts that distance back into (X, Y) coordinates.
closest_pt_on_line = line.interpolate(line.project(target_point))

print(f"The closest point on the line is: ({closest_pt_on_line.x:.4f}, {closest_pt_on_line.y:.4f})")