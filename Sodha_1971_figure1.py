"""
The paper's Eq. 13, written in the dimensionless variable xi = z / R_d0
(where R_d0 = (1/2) a^2 k0(0) is the Rayleigh-like diffraction length),
becomes:

    f'' - (alpha1/2) f'  =  exp(alpha1 * xi) / f^3   -   Phi(xi, f)

with the "focusing" term Phi different for the two physical cases:

    linear    (Eq. 14a):  Phi = gamma1 * f * exp((alpha1 - alpha3) * xi)
    nonlinear (Eq. 14b):  Phi = gamma2 * exp((2*alpha1 - alpha2) * xi) / f^3

Initial conditions: f(0) = 1, f'(0) = 0 (plane wavefront entering the medium).
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# Right-hand side of the first-order system
def rhs(xi, y, alpha1, alpha2, alpha3, gamma, mode):
    f, fp = y[0], y[1]
    diffraction = np.exp(alpha1 * xi) / f**3
    if mode == "linear":
        focusing = gamma * f * np.exp((alpha1 - alpha3) * xi)
    else:
        focusing = (gamma / f**3) * np.exp((alpha1 - alpha2) * xi)
    fpp = 0.5 * alpha1 * fp + diffraction - focusing
    return [fp, fpp]


# Curve
curves = [
    dict(n=1, mode="linear",    alpha1=-0.5, gamma=0.3, xi_end = 2.4, plot_scale = 1,
         color="blue", linestyle="-",  label="1: linear, alpha1=-0.5, gamma1=0.3"),
    dict(n=2, mode="linear",    alpha1= 0.5, gamma=2.0, xi_end = 2.4, plot_scale = 1,
         color="green", linestyle="-",  label="2: linear, alpha1=+0.5, gamma1=2.0"),
    dict(n=3, mode="nonlinear", alpha1=-0.5, gamma=0.3, xi_end = 9.6, plot_scale = 1/4,
         color="red", linestyle="--", label="3: nonlinear, alpha1=-0.5, gamma2=0.3"),
    dict(n=4, mode="nonlinear", alpha1= 0.5, gamma=2.0, xi_end = 2.4, plot_scale = 1,
         color="orange", linestyle="--", label="4: nonlinear, alpha1=+0.5, gamma2=2.0"),
]


# Stop integration at beam collapse 
def collapse_event(xi, y, *args):
    return y[0] - 0.05
collapse_event.terminal = True
collapse_event.direction = -1


# Integration and plotting 
fig, ax = plt.subplots(figsize=(8, 6))

for c in curves:
    a1 = c["alpha1"]
    a2 = 1.5 * a1
    a3 = a1

    sol = solve_ivp(fun=rhs, t_span=[0.0, c["xi_end"]],
        y0=[1.0, 0.0],
        args=(a1, a2, a3, c["gamma"], c["mode"]),
        method="RK45",
        rtol=1e-9, atol=1e-11, max_step=1e-3,
        dense_output=True, events = collapse_event
    )

    xi_ode = np.linspace(0, sol.t[-1], 2000)
    f_vals = sol.sol(xi_ode)[0]
    xi_plot = xi_ode * c["plot_scale"]

    ax.plot(xi_plot, f_vals,
            color=c["color"], linestyle=c["linestyle"],
            linewidth=1.8, label=c["label"])

    print(f"Curve {c['n']} ({c['mode']:>9}): xi_ode ends at {sol.t[-1]:.3f}, f_final = {f_vals[-1]:.3f}")


# Axis styling
ax.axhline(1.0, color="grey", linewidth=0.5, linestyle=":")
ax.set_xlabel(r"$\xi = z / R_{d0}$", fontsize=12)
ax.set_ylabel(r"beam-width parameter  $f$", fontsize=12)
ax.set_xlim(0, 2.4)
ax.set_ylim(0, 3)
ax.set_title("Sodha, Ghatak & Tripathi (1971), Fig. 1", fontsize=13)
ax.legend(fontsize=9, loc="upper left")
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("plots/figure1_reproduction.png", dpi=200, bbox_inches="tight")
plt.show()
