import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
from scipy.optimize import minimize

# Caricamento dati da Excel
filename = r'C:\Users\Andrea Blumer\Dropbox\Attuatore Lineare [Tirocinio Alessio Salvatore]\SCRIPTS\cubic_pm.xlsx'
data = pd.read_excel(filename, header=None).values

# Supponiamo che la prima colonna sia x e la prima riga siano V
x = np.asarray(data[1:, 0])  # Convertiamo esplicitamente x in un array NumPy
V = np.array([0, -12, 12, -24, 24])  # Valori di V (esclusa l'intestazione)
F = data[1:, 1:]  # Forze (matrice con le forze per i vari V)

# Filtro di Savitzky-Golay per ridurre il rumore nei dati
window_size = 5  # Imposta la dimensione della finestra
polynomial_order = 3  # Ordine del polinomio
F_filtered = savgol_filter(F, window_size, polynomial_order, axis=0)

# Interpolazione della componente indipendente da V (f1(x)) usando spline
mask = x <= 15
smoothing_factor = np.std(F_filtered[mask, :])  # Determina un fattore di smoothing basato sulla variabilità dei dati
spline_x = UnivariateSpline(x[mask], np.mean(F_filtered[mask, :], axis=1), s=smoothing_factor)

# Funzione f2(x, V): dipendente da x e V
def f2(alpha, beta, x, V):
    return alpha * x * V * np.exp(-x * V * beta)

# Funzione f3(x): parabola indipendente da V (solo in [15, 23])
def f3(a, b, c, x):
    return a * (x - 15)**2 + b * (x - 15) + c

# Funzione completa F(x, V)
def f(params, x, V):
    alpha, beta, a, b, c = params
    f1_values = spline_x(x)
    f2_values = f2(alpha, beta, x, V)
    f3_values = f3(a, b, c, x)
    x = np.asarray(x)  # Assicuriamoci che x sia un array NumPy per compatibilità con np.where
    return np.where(x <= 15, f1_values + f2_values, spline_x(15) + f2(alpha, beta, 15, V) + f3_values)

# Funzione di errore (MSE)
def mse(params):
    F_pred = f(params, X.ravel(), V_mesh.ravel()).reshape(X.shape)
    return np.mean((F_filtered - F_pred) ** 2)

# Inizializzazione dei parametri
initial_params = [1, 1, 1, 1, 1]

# Creazione della griglia x, V
X, V_mesh = np.meshgrid(x, V, indexing='ij')

# Ottimizzazione dei parametri
result = minimize(mse, initial_params, method='Nelder-Mead', options={'disp': True})
optimized_params = result.x

# Parametri ottimizzati
alpha_opt, beta_opt, a_opt, b_opt, c_opt = optimized_params

# Calcolare le forze predette
F_pred = f(optimized_params, X, V_mesh)

# Visualizzazione del confronto
fig, axes = plt.subplots(5, 1, figsize=(8, 12), sharex=True)
for i, ax in enumerate(axes):
    F_interpolated = np.interp(x, x, F[:, i])
    ax.plot(x, F_interpolated, 'o-', label=f'Dati reali (V = {V[i]})')
    ax.plot(x, F_pred[:, i], 'x-', label=f'Forza predetta (V = {V[i]})')
    ax.set_xlabel('Posizione x')
    ax.set_ylabel('Forza F')
    ax.legend()
    ax.set_title(f'Confronto forze per V = {V[i]}')

plt.tight_layout()
plt.show()
