import cobra
from mimeco import analysis
import pandas as pd
import os
import time

import sys
print(sys.prefix)



def compute_pareto_optimality(
    model,
    rxn1_id: str,
    rxn2_id: str,
    dinc: float = 0.001,
) -> pd.DataFrame:
    """
    Python translation of COBRA Toolbox computeParetoOptimality.m.

    Takes a pre-built joint community model (e.g. from MIMEco) with
    diet already applied, and two biomass reaction IDs.

    Two-pass scan (Oberhardt et al. 2010):
      Pass 1: fix rxn1 at intervals [min1..max1], maximise rxn2
      Pass 2: fix rxn2 at intervals [min2..max2], maximise rxn1

    Returns a DataFrame with columns ['x', 'y'] (rxn1, rxn2 flux values),
    sorted by x, NaNs dropped.
    """
    import numpy as np
    import pandas as pd

    def _bound_opt(m, rxn_id, direction):
        with m:
            m.objective = rxn_id
            m.objective_direction = direction
            sol = m.optimize()
            print(sol)
            return sol.objective_value if sol.status == "optimal" else None

    max1 = _bound_opt(model, rxn1_id, "max")
    min1 = _bound_opt(model, rxn1_id, "min")
    max2 = _bound_opt(model, rxn2_id, "max")
    min2 = _bound_opt(model, rxn2_id, "min")

    print(f"  {rxn1_id} range: [{min1:.4f}, {max1:.4f}]")
    print(f"  {rxn2_id} range: [{min2:.4f}, {max2:.4f}]")

    results = []

    # Pass 1: scan rxn1, maximise rxn2
    print(f"  Pass 1: {len(np.arange(min1, max1, dinc))} steps …")
    for val in np.arange(min1, max1 + dinc * 0.5, dinc):
        with model:
            r1 = model.reactions.get_by_id(rxn1_id)
            r1.lower_bound = val-0.001*val
            r1.upper_bound = val+0.001*val
            model.objective = rxn2_id
            model.objective_direction = "max"
            sol = model.optimize()
            if sol.status == "optimal":
                results.append((sol.fluxes[rxn1_id], sol.fluxes[rxn2_id]))

    # Pass 2: scan rxn2, maximise rxn1
    print(f"  Pass 2: {len(np.arange(min2, max2, dinc))} steps …")
    for val in np.arange(min2, max2 + dinc * 0.5, dinc):
        with model:
            r2 = model.reactions.get_by_id(rxn2_id)
            r2.lower_bound = val-0.001*val
            r2.upper_bound = val+0.001*val
            model.objective = rxn1_id
            model.objective_direction = "max"
            sol = model.optimize()
            if sol.status == "optimal":
                results.append((sol.fluxes[rxn1_id], sol.fluxes[rxn2_id]))

    df = (pd.DataFrame(results, columns=["x", "y"])
            .dropna()
            .sort_values("x")
            .drop_duplicates()
            .reset_index(drop=True))
    return df

t1 = time.time()
os.environ["GRB_LICENSE_FILE"] = "/home/lamberta/gurobi/gurobi.lic"

YPD = pd.read_csv('mimeco_validation/media/YPD_media_ModelSEED.csv', sep=',', header=0, index_col=0)

A_ory = cobra.io.read_sbml_model('/home/lamberta/mimeco_validation/models/A_ory_gapseq.xml')
A_ory.solver = "gurobi"
A_ory.objective = A_ory.reactions.get_by_id("EX_cpd11416_c0").id
L_pla = cobra.io.read_sbml_model('/home/lamberta/mimeco_validation/models/L_pla_gapseq.xml')
L_pla.solver = "gurobi"
L_pla.objective = L_pla.reactions.get_by_id("EX_cpd11416_c0").id

A_ory_obj = A_ory.reactions.get_by_id("EX_cpd11416_c0").id
L_pla_obj = L_pla.reactions.get_by_id("EX_cpd11416_c0").id




cobra_ecosys = analysis.extract_pairwize_model(model1 = A_ory, model2 = L_pla, solver = "gurobi", model1_biomass_id = A_ory_obj,
                                               model2_biomass_id = L_pla_obj, medium = YPD, undescribed_metabolites_constraint = "partially_constrained",
                                               undescribed_met_lb = -0.0001)

# Let MIMEco build the joint model with diet applied
# Then just run the Pareto scan on it
xy = compute_pareto_optimality(
    model=cobra_ecosys,
    rxn1_id="EX_cpd11416_c0:AOryzifermentans",
    rxn2_id="EX_cpd11416_c0:LPlantarum",
    dinc=0.0001,
)
t2 = time.time()
runtime = t2 - t1

print(f"With {len(xy)} data points, MMTpy takes {runtime} seconds to run")
print(xy)


import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

fig, ax = plt.subplots(figsize=(6, 4.5))

ax.scatter(xy["x"], xy["y"], s=5, alpha=0.7)

ax.set_xlabel("Biomass of A. ory")
ax.set_ylabel("Biomass of L. pla")
ax.set_title("Pairwise MMT Pareto Front")

ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax.grid(True, linestyle=":", alpha=0.6)

plt.tight_layout()
plt.show()

#What is we normilize by growth when alone:
# 1. Define the normalization factors
x_factor = 6.194575600148562e-05
y_factor = 7.629860937978175e-05

# 2. Plot directly using the original 'xy' to avoid dataframe manipulation errors
plt.figure(figsize=(6, 4.5))
plt.scatter(xy["x"] / x_factor, xy["y"] / y_factor, s=5, alpha=0.7)

# 3. Add labels so you can verify the axes
plt.xlabel("Biomass of A. ory (Normalized)")
plt.ylabel("Biomass of L. pla (Normalized)")
plt.title("Normalized Pairwise MMT Pareto Front")
plt.grid(True, linestyle=":", alpha=0.6)

plt.show()

#Test growth one bateria:
cobra_ecosys.objective = cobra_ecosys.reactions.get_by_id("EX_cpd11416_c0:AOryzifermentans")
print (cobra_ecosys.optimize())


#Now on another

os.environ["GRB_LICENSE_FILE"] = "/home/lamberta/gurobi/gurobi.lic"

min_medium = pd.read_csv('mimeco_validation/media/min_medium.csv', sep=',', header=0, index_col=0)


min_medium = pd.read_csv('mimeco_validation/media/min_medium.csv', sep=',', header=0, index_col=0)

M_berkei = cobra.io.read_sbml_model('mimeco_validation/models/iAF692.xml')
M_berkei.solver = "gurobi"
print(M_berkei.optimize())

G_meta = cobra.io.read_sbml_model('mimeco_validation/models/iAF987.xml')
G_meta.solver = "gurobi"
        
#Weird upper bounds : 
    
for ex in G_meta.exchanges:
    if ex.upper_bound < 0:
        print(ex)
        print(ex.upper_bound)
        ex.upper_bound = 0

        
G_meta.solver = "gurobi"

M_berkei_obj = M_berkei.reactions.get_by_id("BIOMASS_Mb_30").id
G_meta_obj = G_meta.reactions.get_by_id("BIOMASS_Gm_GS15_WT_79p20M").id

cobra_ecosys = analysis.extract_pairwize_model(model1 = M_berkei, model2 = G_meta, solver = "gurobi", model1_biomass_id = M_berkei_obj,
                                               model2_biomass_id = G_meta_obj, medium = min_medium, undescribed_metabolites_constraint = "partially_constrained",
                                               undescribed_met_lb = -0.1)

xy = compute_pareto_optimality(
    model=cobra_ecosys,
    rxn1_id="BIOMASS_Mb_30:iAF692",
    rxn2_id="BIOMASS_Gm_GS15_WT_79p20M:iAF987",
    dinc=0.0001,
)

plt.figure(figsize=(5,5))
plt.scatter(xy["x"] / x_factor, xy["y"] / y_factor, s=5, alpha=0.7)

# 3. Add labels so you can verify the axes
plt.xlabel("Biomass of M_berkei")
plt.ylabel("Biomass of G_meta")
plt.title("Normalized Pairwise MMT Pareto Front")
plt.grid(True, linestyle=":", alpha=0.6)

plt.show()

x_factor = 0.003423871875361499
y_factor = 0.03012880120670956

plt.figure(figsize=(6, 4.5))
plt.scatter(xy["x"] / x_factor, xy["y"] / y_factor, s=5, alpha=0.7)

# 3. Add labels so you can verify the axes
plt.xlabel("Biomass of M_berkei (Normalized)")
plt.ylabel("Biomass of G_meta (Normalized)")
plt.title("Normalized Pairwise MMT Pareto Front")
plt.grid(True, linestyle=":", alpha=0.6)

plt.show()