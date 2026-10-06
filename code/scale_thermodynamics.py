"""
Thermodynamics of scale: what is exact, numerically checked.

Part A  Entropic c-function of a 1+1D lattice Dirac fermion with mass m (staggered potential):
        c(l) = 3 l dS/dl (Casini-Huerta), from correlation-matrix (Peschel) entanglement entropies.
        Checks strong subadditivity in the form S(l+2)+S(l-2) <= 2 S(l) and the monotone decrease of c.
Part B  Data processing under coarse-graining: relative entropy D(rho_l || sigma_l) between the massive and
        massless ground states restricted to a block of length l is non-decreasing in l (Lindblad-Uhlmann).
Part C  Spectral thermodynamics of scale: with tau as inverse temperature on the Laplacian spectrum,
        S(tau) = ln Z + tau <lambda> = -D(pi_tau || rho) is non-increasing, heat capacity C = tau^2 Var >= 0.
Part D  Jarzynski and Crooks relations for a physical scale protocol: a Gaussian lattice field whose mass
        (inverse correlation length) is changed in time under overdamped Langevin dynamics.

Usage:  python scale_thermodynamics.py [--show]
Figures go to ./figures (PNG, PDF), data to ./data (CSV); override with ST_FIG_DIR, ST_DATA_DIR.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ive

FIG_DIR = os.environ.get("ST_FIG_DIR", "figures")
DATA_DIR = os.environ.get("ST_DATA_DIR", "data")
rng = np.random.default_rng(20261006)


def savefig(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=150)


def savedata(name, columns, cols, meta=()):
    os.makedirs(DATA_DIR, exist_ok=True)
    arr = np.column_stack([np.asarray(c, dtype=float) for c in cols])
    with open(os.path.join(DATA_DIR, f"{name}.csv"), "w", encoding="utf-8", newline="\n") as fh:
        for line in meta:
            fh.write(f"# {line}\n")
        fh.write(",".join(columns) + "\n")
        np.savetxt(fh, arr, delimiter=",", fmt="%.10g")


# ================================================================ Parts A and B: lattice Dirac fermion
def ground_state_correlations(N, m):
    """Half-filled periodic chain H = -sum (c_j^+ c_{j+1} + h.c.) + m sum (-1)^j n_j; returns C_ij = <c_i^+ c_j>."""
    H = np.zeros((N, N))
    i = np.arange(N)
    H[i, (i + 1) % N] = H[(i + 1) % N, i] = -1.0
    H[i, i] = m * (-1.0) ** i
    E, V = np.linalg.eigh(H)
    occ = V[:, : N // 2]
    return occ @ occ.T


def block_entropy(C, l):
    z = np.clip(np.linalg.eigvalsh(C[:l, :l]), 1e-15, 1 - 1e-15)
    return float(-np.sum(z * np.log(z) + (1 - z) * np.log(1 - z)))


def gaussian_relative_entropy(C, G):
    """D(rho||sigma) for fermionic Gaussian states with correlation matrices C (rho) and G (sigma)."""
    zc = np.clip(np.linalg.eigvalsh(C), 1e-15, 1 - 1e-15)
    w, U = np.linalg.eigh(G)
    w = np.clip(w, 1e-15, 1 - 1e-15)
    h = (U * np.log((1 - w) / w)) @ U.T                    # sigma ~ exp(-sum h_ij c_i^+ c_j)
    S_rho = -np.sum(zc * np.log(zc) + (1 - zc) * np.log(1 - zc))
    lnZ = -np.sum(np.log(1 - w))
    return float(-S_rho + np.trace(h @ C) + lnZ)


def parts_ab(N=1200, masses=(0.0, 0.05, 0.1, 0.2), lmax=200):
    ls = np.arange(4, lmax + 1, 2)
    out_c, ssa = {}, {}
    G0 = ground_state_correlations(N, 0.0)
    D_rows = {}
    for m in masses:
        C = ground_state_correlations(N, m)
        S = np.array([block_entropy(C, l) for l in range(2, lmax + 3, 2)])     # l = 2, 4, ..., lmax+2
        c = 3 * ls * (S[2:] - S[:-2]) / 4.0                                      # central difference, step 2
        out_c[m] = (ls, c, S[1:-1])
        ssa[m] = np.max(S[2:] + S[:-2] - 2 * S[1:-1])
        if m > 0:
            D_rows[m] = np.array([gaussian_relative_entropy(C[:l, :l], G0[:l, :l]) for l in ls])
    print(f"[A] lattice Dirac fermion, N={N}, half filling, blocks l = 4..{lmax} (even)")
    for m in masses:
        ls_, c, _ = out_c[m]
        dc = np.diff(c)
        print(f"    m={m:4.2f}: c(4)={c[0]:.3f}  c({ls_[-1]})={c[-1]:.3f}  max SSA violation "
              f"S(l+2)+S(l-2)-2S(l) = {ssa[m]:+.1e}  monotone decrease: {bool(np.all(dc <= 1e-6)) if m > 0 else 'n/a'}"
              f" (largest increase {dc.max():+.1e})")
    print("[B] relative entropy D(rho_l(m) || sigma_l(m=0)) is non-decreasing in l (data processing):")
    for m, D in D_rows.items():
        print(f"    m={m:4.2f}: D(l=4)={D[0]:.4f}  D(l={ls[-1]})={D[-1]:.4f}  min increment {np.diff(D).min():+.1e}")
    cols, names = [ls], ["l"]
    for m in masses:
        cols += [out_c[m][1], out_c[m][2]]; names += [f"c_m{m}", f"S_m{m}"]
    savedata("entropic_c_function", names, cols, [f"lattice Dirac fermion (staggered mass), N={N}, periodic, half filling"])
    savedata("relative_entropy_vs_block", ["l"] + [f"D_m{m}" for m in D_rows], [ls] + list(D_rows.values()),
             ["D(rho_l(m) || sigma_l(0)) for fermionic Gaussian states"])
    return out_c, D_rows, ls


# ================================================================ Part C: spectral thermodynamics of scale
def part_c(d=2):
    tau = np.logspace(-3, 3, 400)
    x = 2 * tau
    r = ive(1, x) / ive(0, x)
    lnZ = d * (np.log(ive(0, x)))                         # ln[e^{-x} I0(x)]^d,  Z = return probability
    mean = d * 2 * (1 - r)                                # <lambda>_tau
    var = 4 * d * (1 - r / x - r ** 2)                     # Var(lambda) = -d<lambda>/dtau, using r' = 1 - r/x - r^2
    S = lnZ + tau * mean                                   # = -D(pi_tau || rho)  <= 0
    C = tau ** 2 * var
    ds = 2 * tau * mean                                    # spectral dimension (Paper 4)
    resid = np.max(np.abs(C - (ds / 2 - 0.5 * np.gradient(ds, np.log(tau))))[2:-2])   # C = d_s/2 - (1/2) dd_s/dln tau
    ip = np.argmax(C)
    print(f"[C] d={d} lattice spectrum: S(tau) = -D(pi_tau||rho) runs from {S[0]:.3e} to {S[-1]:.3f}; "
          f"max dS/dtau = {np.max(np.gradient(S, tau)):+.1e} (<= 0); min heat capacity = {C.min():.3e} (>= 0)")
    print(f"    heat capacity: max C = {C[ip]:.4f} = {C[ip] / d:.4f} d at tau = {tau[ip]:.3f}; C(tau={tau[-1]:.0f}) = "
          f"{C[-1]:.4f} (d/2 = {d / 2}); max |C - [d_s/2 - (1/2) d d_s/d ln tau]| = {resid:.1e}; "
          f"max d_s = {ds.max():.4f} = {ds.max() / d:.4f} d")
    savedata("spectral_thermodynamics", ["tau", "lnZ", "mean_lambda", "var_lambda", "S", "C"],
             [tau, lnZ, mean, var, S, C], [f"{d}D square lattice Laplacian spectrum, tau = inverse temperature"])
    return tau, S, C


# ================================================================ Part D: Jarzynski / Crooks for a scale protocol
def part_d(L=16, m2_i=0.01, m2_f=1.0, durations=(0.2, 2.0, 20.0), M=100000, steps=400):
    """Forward protocol stiffens the field (correlation length 10 -> 1). For the reverse (softening) protocol the
    estimator <e^{-W}> has infinite variance: for a sudden quench of a mode, E[e^{-2W}] < inf iff w_f^2/w_i^2 > 1/2."""
    k = np.arange(L)
    q2 = 4 * np.sin(np.pi * k / L) ** 2                    # lattice Laplacian eigenvalues on the ring
    dF = 0.5 * np.sum(np.log((q2 + m2_f) / (q2 + m2_i)))  # beta = 1
    ratio_rev = np.min((q2 + m2_i) / (q2 + m2_f))
    print(f"[D] Gaussian ring field, L={L}, m^2: {m2_i} -> {m2_f} (correlation length 10 -> 1); Delta F = {dF:.4f}")
    print(f"    reverse protocol: min w_f^2/w_i^2 = {ratio_rev:.3f} < 1/2, so the reverse Jarzynski estimator has "
          f"infinite variance in the sudden limit")

    def run(m2_start, m2_end, T):
        m2 = np.linspace(m2_start, m2_end, steps + 1)
        dt = T / steps
        w2 = q2 + m2[0]
        phi = rng.normal(size=(M, L)) / np.sqrt(w2)
        W = np.zeros(M)
        for n in range(steps):
            W += 0.5 * (m2[n + 1] - m2[n]) * np.sum(phi ** 2, axis=1)   # work: switch the parameter
            w2 = q2 + m2[n + 1]
            a = np.exp(-w2 * dt)
            phi = phi * a + np.sqrt((1 - a ** 2) / w2) * rng.normal(size=(M, L))  # exact OU relaxation
        return W

    rows, hists = [], {}
    for T in durations:
        Wf = run(m2_i, m2_f, T)
        Wr = run(m2_f, m2_i, T)
        ew = np.exp(-(Wf - dF))
        jar = -np.log(np.mean(ew)) + dF
        jar_se = np.std(ew) / (np.sqrt(M) * np.mean(ew))   # delta-method standard error of the Jarzynski estimate
        # Crooks: ln[P_F(W)/P_R(-W)] = W - Delta F; fit over the overlap region
        lo, hi = max(Wf.min(), (-Wr).min()), min(Wf.max(), (-Wr).max())
        bins = np.linspace(lo, hi, 41)
        pf, _ = np.histogram(Wf, bins, density=True); pr, _ = np.histogram(-Wr, bins, density=True)
        ctr = 0.5 * (bins[1:] + bins[:-1]); ok = (pf > 0) & (pr > 0) & (np.minimum(pf, pr) * np.diff(bins) * M > 50)
        slope, icpt = np.polyfit(ctr[ok], np.log(pf[ok] / pr[ok]), 1) if ok.sum() > 3 else (np.nan, np.nan)
        # Bennett acceptance ratio (optimal use of the Crooks relation, equal sample sizes)
        from scipy.optimize import brentq
        g = lambda Cc: np.sum(1 / (1 + np.exp(np.clip(Wf - Cc, -700, 700)))) - np.sum(1 / (1 + np.exp(np.clip(Wr + Cc, -700, 700))))
        bar = brentq(g, Wf.min() - 50, Wf.max() + 50)
        rows.append([T, Wf.mean(), Wf.mean() - dF, jar, jar_se, slope, bar])
        hists[T] = Wf
        print(f"    T={T:5.1f}: <W>={Wf.mean():+.4f} (dissipated {Wf.mean() - dF:.4f} >= 0); Jarzynski estimate "
              f"{jar:+.4f} +- {jar_se:.4f}; Crooks: slope of ln[P_F(W)/P_R(-W)] vs W = {slope:.3f} (expect 1), Bennett Delta F = {bar:+.4f}")
    R = np.array(rows)
    savedata("jarzynski_crooks", ["T", "mean_W", "dissipated_W", "jarzynski_dF", "jarzynski_se", "crooks_slope", "bennett_dF"],
             [R[:, i] for i in range(7)], [f"Gaussian ring field L={L}, m^2 {m2_i}->{m2_f}, exact Delta F={dF:.6f}, M={M}"])
    return dF, R, hists


def main():
    plt.rcParams["font.family"] = "DejaVu Sans"
    out_c, D_rows, ls = parts_ab()
    tau, S, C = part_c()
    dF, R, hists = part_d()

    fig, ax = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    for m, (l_, c, _) in out_c.items():
        ax[0, 0].plot(l_, c, label=f"m = {m}" + ("  (CFT, c = 1)" if m == 0 else f"  (xi = 2/m = {2 / m:.0f})"))
    ax[0, 0].axhline(1, color="grey", lw=0.5)
    ax[0, 0].set_xlabel("block length l"); ax[0, 0].set_ylabel("c(l) = 3 l dS/dl")
    ax[0, 0].set_title("Entropic c-function: monotone decrease along the flow to the IR"); ax[0, 0].legend(fontsize=8)
    for m, D in D_rows.items():
        ax[0, 1].plot(ls, D, label=f"m = {m}")
    ax[0, 1].set_xlabel("block length l"); ax[0, 1].set_ylabel("D(rho_l(m) || sigma_l(0))")
    ax[0, 1].set_title("Data processing: coarse-graining never increases distinguishability")
    ax[0, 1].legend(fontsize=8)
    a = ax[1, 0]
    a.semilogx(tau, S, label="S(tau) = -D(pi_tau || rho)")
    a2 = a.twinx(); a2.semilogx(tau, C, color="tab:red", label="C(tau) = tau^2 Var(lambda)")
    a.set_xlabel("tau (scale as inverse temperature)"); a.set_ylabel("S(tau)"); a2.set_ylabel("C(tau)", color="tab:red")
    a.set_title("Spectral thermodynamics of scale (2D lattice)")
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)
    a = ax[1, 1]
    for T, W in hists.items():
        a.hist(W, bins=np.linspace(0, 150, 151), density=True, histtype="step", label=f"T = {T}")
    a.axvline(dF, color="k", ls="--", label=f"Delta F = {dF:.3f}")
    a.set_xlabel("work W (units of k_B T)"); a.set_ylabel("P(W)")
    a.set_title("Scale protocol xi: 10 -> 1 (Jarzynski, Crooks)")
    a.legend(fontsize=8)
    savefig(fig, "scale_thermodynamics")
    print(f"figures written to {FIG_DIR}/, data to {DATA_DIR}/")
    if "--show" in sys.argv:
        plt.show()


if __name__ == "__main__":
    main()
