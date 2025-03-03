#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import pandas as pd
import pulp
from math import log, exp
from scipy.optimize import curve_fit
from scipy.stats import linregress
import os
import json

# Constants
R = 8.314  # J/(mol·K)

def parse_reaction(reaction_str):
    if '<=>' in reaction_str:
        sides = reaction_str.split('<=>')
    elif '=>' in reaction_str:
        sides = reaction_str.split('=>')
    else:
        raise ValueError("Invalid reaction string format. Use either '=>' for one-directional or '<=>' for bidirectional reactions.")
    reactants, products = sides
    return parse_side(reactants), parse_side(products)

def parse_side(side):
    compounds = side.split('+')
    result = {}
    for compound in compounds:
        compound = compound.strip()
        if ' ' in compound:
            parts = compound.split(' ', 1)
            if parts[0].replace('.', '', 1).isdigit():
                coeff = float(parts[0])
                metabolite = parts[1].strip()
            else:
                coeff = 1.0
                metabolite = compound.strip()
        else:
            coeff = 1.0
            metabolite = compound.strip()
        result[metabolite] = coeff
    return result

def mdf_analysis(reaction_str, delta_G0, RT):
    reactants, products = parse_reaction(reaction_str)
    prob = pulp.LpProblem("Max_min_Driving_Force_Optimization", pulp.LpMaximize)
    metabolites = set(list(reactants.keys()) + list(products.keys()))
    conc_vars = {met: pulp.LpVariable(f"log_conc_{met}", lowBound=log(1e-6), upBound=log(1e-2)) for met in metabolites}
    B = pulp.LpVariable("B")
    prob += B
    Q_term = pulp.lpSum(conc_vars[met] * coeff for met, coeff in products.items()) -              pulp.lpSum(conc_vars[met] * coeff for met, coeff in reactants.items())
    prob += -delta_G0 - RT * Q_term >= B
    prob.solve()
    optimal_conc = {met: exp(conc_vars[met].varValue) for met in conc_vars}
    return optimal_conc, reactants
def michaelis_menten(S, Vmax, Km):
    return Vmax * S / (Km + S)

def get_kcat_values(rxn_id, kcat_df):
    kcat_values = {}
    matching_rows = kcat_df[kcat_df['rxn_id'] == rxn_id]
    if matching_rows.empty:
        return kcat_values
    start_index = matching_rows.index[0]
    for i in range(start_index + 1, len(kcat_df)):
        if pd.notna(kcat_df.loc[i, 'rxn_id']):
            break
        substrate = kcat_df.loc[i, 'Substrate Name']
        kcat = kcat_df.loc[i, 'Kcat value (1/s)']
        if pd.notna(kcat):
            kcat_values[substrate] = kcat
    return kcat_values

def calculate_lb_ub(row, lb_ub_df):
    rxn_id = row['rxn_id']
    v_optimum = row['rate_limiting_v_optimum']

    if rxn_id not in lb_ub_df.index:
        return pd.Series({'lb': None, 'ub': None})

    model_lb = lb_ub_df.loc[rxn_id, 'Model_lb']
    model_ub = lb_ub_df.loc[rxn_id, 'Model_ub']

    if v_optimum is None:
        # If v_optimum is None, return the original bounds
        return pd.Series({'lb': model_lb, 'ub': model_ub})

    if model_lb == 0 and model_ub == 1000:
        lb = 0
        ub = min(abs(model_ub), abs(v_optimum*3600))
    elif model_lb == -1000 and model_ub == 1000:
        lb = -1 * min(abs(model_ub), abs(v_optimum*3600))
        ub = min(abs(model_ub), abs(v_optimum*3600))
    else:
        lb = model_lb
        ub = model_ub

    return pd.Series({'lb': lb, 'ub': ub})

def generate_initial_estimates(v_max, num_estimates=5):
    return np.linspace(v_max * 0.5, v_max * 1.5, num_estimates)

def generate_synthetic_data(S, Vmax, Km, noise_level=0.3):
    v = michaelis_menten(S, Vmax, Km)
    noise = np.random.normal(0, np.abs(v) * noise_level, v.shape)
    v_noisy = v + noise
    return v_noisy

def fit_michaelis_menten(S, v, initial_estimates):
    best_fit = None
    best_r_squared = -np.inf
    for Vmax_estimate in initial_estimates:
        try:
            if Vmax_estimate <= 0 or np.median(S) <= 0:
                continue
            popt, _ = curve_fit(michaelis_menten, S, v, p0=[Vmax_estimate, np.median(S)], bounds=([0, 0], [np.inf, np.inf]))
            Vmax_fit, Km_fit = popt
            v_fit = michaelis_menten(S, Vmax_fit, Km_fit)
            slope, intercept, r_value, _, _ = linregress(v, v_fit)
            r_squared = r_value ** 2
            if r_squared > best_r_squared:
                best_fit = (Vmax_fit, Km_fit)
                best_r_squared = r_squared
        except RuntimeError:
            continue
    return best_fit, best_r_squared


def process_model(model_name, temperature, base_path):
    T = temperature + 273.15  # Convert Celsius to Kelvin
    RT = R * T

    # Construct file paths
    test_reactions_file = os.path.join(base_path, "test_reactions", f"test_rxns_{model_name}.xlsx")
    fva_dlkcat_file = os.path.join(base_path, "FVA and DLKcat", f"FVA_and_DLKcat_{model_name}.xlsx")
    lb_ub_file = os.path.join(base_path, "Lb_Ub_excel", f"{model_name}_Lb_and_Ub.xlsx")

    # Read input files
    main_df = pd.read_excel(test_reactions_file, sheet_name=0)
    fva_df = pd.read_excel(fva_dlkcat_file, sheet_name=0)
    kcat_df = pd.read_excel(fva_dlkcat_file, sheet_name=1)
    lb_ub_df = pd.read_excel(lb_ub_file, index_col='rxn_id')

    results_list = []

    for index, row in main_df.iterrows():
        rxn_id = row['rxn_id']
        reaction_str = row['Reaction Formula']
        delta_G0 = row['Standard_dG']
        
        print(f"Processing reaction: {rxn_id}")

        fva_row = fva_df[fva_df['rxn_id'] == rxn_id]
        if fva_row.empty:
            print(f"Warning: No FVA data found for reaction {rxn_id}")
            results_list.append({
                'rxn_id': rxn_id,
                'v_min': None,
                'v_max': None,
                'Rate-limiting substrate': None,
                'Rate-limiting kcat/Km': None,
                'Enzyme concentration (M)': None,
                'Rate-limiting Vmax': None,
                'Rate-limiting Vmin': None,
                'Metabolite Results': None,
                'Standard_dG': delta_G0,
                'lb': None,
                'ub': None,
                'rate_limiting_v_optimum': None
            })
            continue

        v_min = fva_row['v_min'].values[0] / 3600
        v_max = fva_row['v_max'].values[0] / 3600

        kcat_values = get_kcat_values(rxn_id, kcat_df)
        optimal_conc, reactants = mdf_analysis(reaction_str, delta_G0, RT)

        results = {}
        for substrate, coeff in reactants.items():
            kcat = kcat_values.get(substrate)
            if kcat is None or pd.isna(kcat):
                print(f"Warning: No valid kcat value found for {substrate} in reaction {rxn_id}")
                results[substrate] = {'Vmax': None, 'Km': None, 'kcat/Km': None, 'kcat': None, 'R-squared': None}
                continue

            try:
                kcat = float(kcat)
            except ValueError:
                print(f"Warning: Invalid kcat value '{kcat}' for {substrate} in reaction {rxn_id}")
                results[substrate] = {'Vmax': None, 'Km': None, 'kcat/Km': None, 'kcat': None, 'R-squared': None}
                continue

            S = np.logspace(np.log10(optimal_conc[substrate]/10), np.log10(optimal_conc[substrate]*10), 20)
            
            try:
                initial_estimates = generate_initial_estimates(v_max)
                v_noisy = generate_synthetic_data(S, v_max, optimal_conc[substrate]/2, noise_level=0.3)
                best_fit, r_squared = fit_michaelis_menten(S, v_noisy, initial_estimates)

                if best_fit is not None:
                    Vmax_fit, Km_fit = best_fit
                    kcat_km = kcat / Km_fit
                    results[substrate] = {
                        'Vmax': Vmax_fit,
                        'Km': Km_fit,
                        'kcat/Km': kcat_km,
                        'kcat': kcat,
                        'R-squared': r_squared
                    }
                else:
                    print(f"Warning: Curve fitting failed for {substrate} in reaction {rxn_id}")
                    results[substrate] = {'Vmax': None, 'Km': None, 'kcat/Km': None, 'kcat': None, 'R-squared': None}
            except Exception as e:
                print(f"Error processing {substrate} in reaction {rxn_id}: {str(e)}")
                results[substrate] = {'Vmax': None, 'Km': None, 'kcat/Km': None, 'kcat': None, 'R-squared': None}

        if results:
            valid_results = {k: v for k, v in results.items() if v['kcat/Km'] is not None}
            if valid_results:
                rate_limiting = min(valid_results, key=lambda x: valid_results[x]['kcat/Km'])
                rate_limiting_kcat_km = valid_results[rate_limiting]['kcat/Km']
                rate_limiting_kcat = valid_results[rate_limiting]['kcat']
                Vmax_rl = valid_results[rate_limiting]['Vmax']
                enzyme_conc = Vmax_rl / rate_limiting_kcat
                rate_limiting_v_max = Vmax_rl
                rate_limiting_v_min = v_min
                rate_limiting_v_optimum = michaelis_menten(optimal_conc[rate_limiting], Vmax_rl, valid_results[rate_limiting]['Km'])
            else:
                rate_limiting = "N/A"
                rate_limiting_kcat_km = None
                enzyme_conc = None
                rate_limiting_v_max = v_max
                rate_limiting_v_min = v_min
                rate_limiting_v_optimum = None

            new_row = {
                'rxn_id': rxn_id,
                'v_min': v_min,
                'v_max': v_max,
                'Rate-limiting substrate': rate_limiting,
                'Rate-limiting kcat/Km': rate_limiting_kcat_km,
                'Enzyme concentration (M)': enzyme_conc,
                'Rate-limiting Vmax': rate_limiting_v_max,
                'Rate-limiting Vmin': rate_limiting_v_min,
                'Metabolite Results': json.dumps(results),
                'Standard_dG': delta_G0,
                'rate_limiting_v_optimum': rate_limiting_v_optimum
            }
        else:
            print(f"Warning: No results for reaction {rxn_id}")
            new_row = {
                'rxn_id': rxn_id,
                'v_min': v_min,
                'v_max': v_max,
                'Rate-limiting substrate': None,
                'Rate-limiting kcat/Km': None,
                'Enzyme concentration (M)': None,
                'Rate-limiting Vmax': None,
                'Rate-limiting Vmin': None,
                'Metabolite Results': None,
                'Standard_dG': delta_G0,
                'rate_limiting_v_optimum': None
            }

        results_list.append(new_row)

    results_df = pd.DataFrame(results_list)
    results_df[['lb', 'ub']] = results_df.apply(lambda row: calculate_lb_ub(row, lb_ub_df), axis=1)

    # Create folder for results
    folder_name = f"results_{model_name}_{temperature}C"
    folder_path = os.path.join(base_path, folder_name)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # Save results
    results_file = os.path.join(folder_path, f'{model_name}_results_{temperature}C.xlsx')
    results_df.to_excel(results_file, index=False)
    print(f"Results saved to '{results_file}'")

    new_bounds_file = os.path.join(folder_path, f'{model_name}_new_bounds_{temperature}C.xlsx')
    new_bounds_df = results_df[['rxn_id', 'lb', 'ub']]
    new_bounds_df.to_excel(new_bounds_file, index=False)
    print(f"New bounds saved to '{new_bounds_file}'")

def main():
    base_path = "/lustre/work/ssbio/mtabibian2/bsmh/should_Final/final_kinetic_code_all_models_lung"
    temperatures = [36, 37, 38, 39, 40]
    model_types = ["Patient"]
    model_numbers = range(1, 26)  # 1 to 43

    for model_type in model_types:
        for model_number in model_numbers:
            model_name = f"Model_{model_type}_{model_number}"
            print(f"\nProcessing model: {model_name}")
            for temp in temperatures:
                print(f"Processing temperature: {temp}°C")
                process_model(model_name, temp, base_path)

    print("\nAll models and temperatures processed successfully.")

if __name__ == "__main__":
    main()


# In[ ]:




