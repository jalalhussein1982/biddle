import itertools
import csv
import math
from tqdm import tqdm

# Updated computational functions based on calibrated equations

def get_delta_b(t, B, fr, Vr, Ps, wth):
    t_reserve = wth / Vr
    if 0 < t < t_reserve:
        return (B * fr * Vr * Ps) / wth
    else:
        return 0

def compute_outcomes(params, const):
    t = params['t']
    R, B, TR, TB, d, fr, fe, Vr, Va, wa, wth, t_pre = [params[k] for k in input_vars.keys()]

    TR_index = (TR - 1900) / 10.0
    TB_index = (TB - 1900) / 10.0

    TC = (TB_index ** 2) / TR_index if TR_index != 0 else 0
    T_dash = (TR_index ** 2) / TB_index if TB_index != 0 else 0

    Ps = TR_index ** (-const['k2'] * Vr) if TR_index > 0 else 0

    H = const['k1'] * (1 - fe)

    phi1 = (const['k9'] * B * fr * Ps) / (T_dash ** const['k4']) if T_dash > 0 else 0
    phi2 = (const['k3'] * B * (1 - fr)) / wth

    r0 = R - phi2 * (wth - wa)
    b0 = (B * (1 - fr) * wa) / (wth * d)

    effective_Va = Va + const['k8']
    Ca = const['k7'] * (1 - fe) * TC * b0 * effective_Va

    delta_r = Ca * Va + 2 * phi1 * Va

    t_reserve = wth / Vr

    denominator = delta_r + H * get_delta_b(t, B, fr, Vr, Ps, wth)
    if denominator != 0:
        t_candidate = (r0 - H * b0) / denominator
    else:
        t_candidate = float('inf')

    if t_candidate >= t_reserve:
        delta_b = get_delta_b(t_candidate, B, fr, Vr, Ps, wth)
        t_star = (r0 - H * b0) / (delta_r + H * delta_b)
    else:
        b_final = H + B * fr * Ps
        b_final1 = H * b0
        t_star = (r0 - b_final1 - b_final) / delta_r if delta_r != 0 else float('inf')

    t_star = max(t_star, 0)

    G = t_star * Va
    breakthrough = 1 if G >= d else 0

    CR = min(Ca * G + const['k5'], R)
    CB = min(b0 * G + (1 - Ps) * get_delta_b(t_star, B, fr, Vr, Ps, wth) * t_star + const['k6'], B)

    return {"t_star": t_star, "G": G, "breakthrough": breakthrough, "CR": CR, "CB": CB}

# Updated input handling to support multiple explicit values and range input
def parse_range(input_str, default):
    if input_str.strip() == '':
        return [default]
    # Split by '-' and remove any extra spaces
    parts = [p.strip() for p in input_str.split('-') if p.strip() != '']
    # If only one value is provided, return it as a list
    if len(parts) == 1:
        return [float(parts[0])]
    # If exactly three parts, assume a range specification (start, end, step)
    elif len(parts) == 3:
        start, end, step = map(float, parts)
        if step <= 0:
            raise ValueError("Step must be positive")
        values = []
        while start <= end:
            values.append(start)
            start += step
        return values
    else:
        # Otherwise, treat as an explicit list of values
        return [float(p) for p in parts]

input_vars = {
    "R": 100000, "B": 60000, "TR": 1973, "TB": 1961,
    "d": 15.0, "fr": 0.7, "fe": 0.9, "Vr": 50.0,
    "Va": 3.0, "wa": 6.0, "wth": 500.0, "t": 1.0
}

constants = {"k1": 2.5, "k2": 0.01, "k3": 0.4, "k4": 0.5,
             "k5": 0, "k6": 0, "k7": 5, "k8": 0.1, "k9": 0.01}

scenario_values = {}
for var, default in input_vars.items():
    user_input = input(f"Enter {var} [{default}] (single value, range start-end-step, or explicit values separated by '-'): ").strip()
    try:
        scenario_values[var] = parse_range(user_input, default)
    except ValueError as e:
        print(f"Error processing input for {var}: {e}")
        exit(1)

all_scenarios = list(itertools.product(*scenario_values.values()))
print(f"Generated {len(all_scenarios)} scenarios.")

results = []
# Loop through each scenario, print the outcome, and store the results
for scenario in tqdm(all_scenarios, desc="Running scenarios"):
    params = dict(zip(input_vars.keys(), scenario))
    outcome = compute_outcomes(params, constants)
    result = {**params, **outcome}
    print(result)
    results.append(result)

# Define the order of columns: independent variables first, then dependent variables
fieldnames = list(input_vars.keys()) + ["t_star", "G", "breakthrough", "CR", "CB"]

with open("biddle_model_results.csv", mode='w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print("All scenario results saved to biddle_model_results.csv.")
