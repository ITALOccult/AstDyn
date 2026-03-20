import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Create plots directory
os.makedirs("benchmarks/plots", exist_ok=True)

def plot_convergence():
    if not os.path.exists("benchmarks/convergence_test.csv"):
        return
    df = pd.read_csv("benchmarks/convergence_test.csv")
    plt.figure(figsize=(10, 6))
    plt.loglog(df['epsilon'], df['relative_pos_error'], 'o-', label='AAS Convergence')
    
    # Plot slope reference (M^-4)
    eps_range = np.array([df['epsilon'].min(), df['epsilon'].max()])
    slope_ref = 1e-6 * (eps_range / 1e-4)**4
    plt.loglog(eps_range, slope_ref, '--', color='gray', label='Slope 4 (Theoretical)')
    
    plt.xlabel('Precision Parameter ($\epsilon$)')
    plt.ylabel('Relative Position Error')
    plt.title('AAS Integrator Convergence (Ceres, 1 Orbit)')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    plt.savefig("benchmarks/plots/convergence_test.png")
    plt.close()

def plot_work_precision():
    if not os.path.exists("benchmarks/work_precision.csv"):
        return
    df = pd.read_csv("benchmarks/work_precision.csv")
    plt.figure(figsize=(10, 6))
    
    for ecc in df['eccentricity'].unique():
        subset = df[df['eccentricity'] == ecc]
        for integrator in subset['integrator_name'].unique():
            data = subset[subset['integrator_name'] == integrator]
            plt.loglog(data['nfe'], data['final_error'], 'o-', label=f'{integrator} (e={ecc})')
            
    plt.xlabel('Number of Force Evaluations (NFE)')
    plt.ylabel('Final Position Error [m]')
    plt.title('Work-Precision Diagram')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    plt.savefig("benchmarks/plots/work_precision.png")
    plt.close()

def plot_energy():
    if not os.path.exists("benchmarks/energy_conservation.csv"):
        return
    df = pd.read_csv("benchmarks/energy_conservation.csv")
    plt.figure(figsize=(12, 6))
    plt.semilogy(df['time_years'], df['delta_H'], label='Total Energy Error ($\delta H$)')
    plt.semilogy(df['time_years'], df['delta_H_shadow'], label='Shadow Hamiltonian Error ($\delta \\tilde{H}$)', alpha=0.7)
    
    plt.xlabel('Time [Years]')
    plt.ylabel('Relative Energy Error')
    plt.title('Long-term Energy Conservation (AAS)')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    plt.savefig("benchmarks/plots/energy_conservation.png")
    plt.close()

def plot_stm():
    if not os.path.exists("benchmarks/stm_validation.csv"):
        return
    df = pd.read_csv("benchmarks/stm_validation.csv")
    plt.figure(figsize=(10, 6))
    plt.semilogy(df['time'], df['det_error'], label='Det(Phi) - 1')
    plt.semilogy(df['time'], df['diff_with_FD'], label='STM vs Finite Differences', alpha=0.7)
    
    plt.xlabel('Time [Days]')
    plt.ylabel('Error')
    plt.title('STM Validation (Symplecticity & Accuracy)')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    plt.savefig("benchmarks/plots/stm_validation.png")
    plt.close()

if __name__ == "__main__":
    print("Generating benchmark plots...")
    plot_convergence()
    plot_work_precision()
    plot_energy()
    plot_stm()
    print("Plots saved in benchmarks/plots/")
