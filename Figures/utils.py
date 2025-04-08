import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d, UnivariateSpline
import os


# --------------------------- utils.py ---------------------------
def load_experimental_data(experiments):
    """
    Takes in a list of csv files that contains data from a single experiment and returns a megacsv file with all the data and ajdusting the particle ID
    """
    df_list = []
    particle_offset = 0
    
    for exp in experiments:
        df = pd.read_csv(exp)
        df['particle'] += particle_offset
        particle_offset += df['particle'].max() + 1
        df_list.append(df)
    
    combined_df = pd.concat(df_list, ignore_index=True)
    unique_particles = sorted(combined_df['particle'].unique())
    particle_mapping = {old_id: new_id for new_id, old_id in enumerate(unique_particles)}
    # Map the 'particle' column to new IDs based on the mapping
    combined_df['particle'] = combined_df['particle'].map(particle_mapping)
    return combined_df

# --------------------------- uncertainty.py ---------------------------
def compute_ph_sym_uncertainty(df, ):
    df['pH_sigma'] = ((df['pHhigh'] - df['pH']) * (df['pH'] - df['pHlow']))**0.5
    return df

# --------------------------- filtering.py ---------------------------
def filter_particles(df, min_radius, consider_target=True, target_pH=7.5):
    filtered_particles = df[(df['radius_ref_pH'] > min_radius) & (df['consider_drop'] == 1)]['particle'].unique()
    return df[df['particle'].isin(filtered_particles)].copy()

# --------------------------- fitting.py ---------------------------
def parabola(t, a, b, c):
    return a * t**2 + b * t + c

def fit_ph_rate(trajectory_df, t0, fit_range, initial_guess, tolerance, min_pts = 4):
    condition = trajectory_df[(trajectory_df['pH'] > fit_range[0]) & (trajectory_df['pH'] < fit_range[1])]
    if len(condition) > min_pts:
        time = condition['time(sec)']
        pH = condition['pH']
        try:
            popt, pcov = curve_fit(
                parabola,
                time - t0,
                pH,
                p0=initial_guess,
                method='trf',
                ftol=tolerance,
                xtol=tolerance,
                gtol=tolerance
            )
            return popt, np.sqrt(np.diag(pcov))  # parameters and uncertainties
        except:
            return [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]
    return [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]

def get_ref_radius_and_time(trajectory_df, reference_pH):
    f_inverse = interp1d(trajectory_df['pH'], trajectory_df['time(sec)'], kind='linear', bounds_error=False, fill_value=np.nan)
    t0 = f_inverse(reference_pH)

    if np.isnan(t0):
        r0 = np.nan
    else:
        r_func = interp1d(trajectory_df['time(sec)'], trajectory_df['radius'], kind='linear', bounds_error=False, fill_value="extrapolate")
        r0 = r_func(t0)
        
        if np.isnan(r0):      
            t0 = np.nan

    return r0, t0

def get_time_for_vacuole_formation(trajectory_df, pH_v):
    if np.isnan(pH_v):
        tv = np.nan
    else:
        mask = trajectory_df['pH']<=pH_v
        if np.any(mask):
            tv = trajectory_df['time(sec)'][mask].iloc[0]
        else:
            tv = np.nan
    return tv

def get_uncertainity_in_pHv(trajectory_df, pH_v):
    if np.isnan(pH_v):
        sigmap = np.nan
        sigmam = np.nan
    else:
        mask = trajectory_df['pH']<=pH_v
        if np.any(mask):
            sigmap = trajectory_df['pHhigh'][mask].iloc[0] - trajectory_df['pH'][mask].iloc[0]
            sigmam = trajectory_df['pH'][mask].iloc[0] - trajectory_df['pHlow'][mask].iloc[0]
        else:
            sigmap = np.nan
            sigmam = np.nan
    return sigmap, sigmam

# --------------------------- plotting.py ---------------------------
def plot_fitting(data, popt, t0, pH_vacuole, fit_range, save_path, title):
    plt.figure()
    plt.plot(data['time(sec)'] - t0, data['pH'], 'o', label='data', fillstyle='none')
    plt.fill_between(data['time(sec)'] - t0, data['pHhigh'], data['pHlow'], alpha=0.5)
    plt.plot(data['time(sec)'] - t0, parabola(data['time(sec)'] - t0, *popt), label='fit', color='tab:red')
    plt.axhline(y=pH_vacuole, alpha=0.6, color='k', linestyle='--')
    plt.axhline(y=fit_range[1], alpha=0.6, color='k', linestyle='--')
    plt.ylim(7, 10)
    plt.text(0.2, 0.2, f'rate (min^-1)= {np.abs(popt[1])*60}', transform=plt.gca().transAxes)
    plt.legend()
    plt.title(f'{title}')
    plt.xlabel('Time (sec)')
    plt.ylabel('pH')
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, bbox_inches='tight')
    plt.close()
