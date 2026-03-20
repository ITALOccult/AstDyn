import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Configurazione stile grafico scientifico
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 16,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "lines.linewidth": 1.5,
    "lines.markersize": 6
})

def plot_convergence():
    """Genera il grafico di convergenza (M3 - Slope Test)."""
    df = pd.read_csv("convergence_test.csv")
    plt.figure(figsize=(8, 6))
    
    # Filtriamo eventuali errori nulli per il log-log
    df = df[df["relative_pos_error"] > 0]
    
    plt.loglog(df["epsilon"], df["relative_pos_error"], "o-", label="AAS Position Error")
    
    # Linea di riferimento per ordine 4 (opzionale)
    x = df["epsilon"].values
    y_ref = (x**4) * (df["relative_pos_error"].iloc[0] / (x[0]**4))
    plt.loglog(x, y_ref, "--", color="gray", alpha=0.5, label="Order 4 Reference")
    
    plt.title("Convergence Analysis (AAS)")
    plt.xlabel(r"Precision Parameter $\varepsilon$")
    plt.ylabel("Relative Position Error")
    plt.legend()
    plt.savefig("convergence.pdf")
    plt.close()

def plot_work_precision():
    """Genera il diagramma Work-Precision (M3 - Efficiency)."""
    df = pd.read_csv("work_precision.csv")
    plt.figure(figsize=(8, 6))
    
    for name, group in df.groupby(["integrator_name", "eccentricity"]):
        label = f"{name[0]} (e={name[1]})"
        plt.loglog(group["nfe"], group["final_error"], "o-", label=label)
    
    plt.title("Work-Precision Diagram")
    plt.xlabel("Number of Force Evaluations (NFE)")
    plt.ylabel("Final Position Error [km]")
    plt.legend()
    plt.savefig("efficiency_comparison.pdf")
    plt.close()

def plot_energy():
    """Genera il grafico della conservazione dell'energia (M1 & M2)."""
    df = pd.read_csv("energy_conservation.csv")
    plt.figure(figsize=(8, 6))
    
    plt.semilogy(df["time_years"], df["delta_H"], label=r"Hamiltonian Error $\delta H$")
    plt.semilogy(df["time_years"], df["delta_H_shadow"], label=r"Shadow Hamiltonian Error $\delta \tilde{H}$")
    
    plt.title("Long-term Energy Conservation")
    plt.xlabel("Time [years]")
    plt.ylabel("Relative Energy Error")
    plt.legend()
    plt.savefig("energy_stability.pdf")
    plt.close()

def plot_stm():
    """Genera il grafico della simpletticità della STM (M4)."""
    df = pd.read_csv("stm_validation.csv")
    plt.figure(figsize=(8, 6))
    
    plt.semilogy(df["time"], df["det_error"], label=r"$|\det \Phi - 1|$")
    
    plt.title("STM Symplecticity Verification")
    plt.xlabel("Time [days]")
    plt.ylabel("Determinant Error")
    plt.legend()
    plt.savefig("stm_determinant.pdf")
    plt.close()

# Esecuzione generazione grafici
if __name__ == "__main__":
    plot_convergence()
    plot_work_precision()
    plot_energy()
    plot_stm()
    print("Grafici generati con successo: convergence.pdf, efficiency_comparison.pdf, energy_stability.pdf, stm_determinant.pdf")
