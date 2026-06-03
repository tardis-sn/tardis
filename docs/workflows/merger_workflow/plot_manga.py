import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, SymLogNorm


def radial_edges(r_vals):
    mids = 0.5 * (r_vals[:-1] + r_vals[1:])
    lo   = r_vals[0]  - (mids[0]  - r_vals[0])
    hi   = r_vals[-1] + (r_vals[-1] - mids[-1])
    return np.concatenate([[lo], mids, [hi]])


def log_compress(Z, diverging, pct_lo=2, pct_hi=98):
    if diverging:
        abs_vals  = np.abs(Z[Z != 0].ravel())
        linthresh = float(np.median(abs_vals)) if abs_vals.size else 1.0
        vmax      = float(np.percentile(abs_vals, pct_hi)) if abs_vals.size else 1.0
        norm      = SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax)
        C         = Z  # SymLogNorm works on raw values
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            log_pos = np.log10(Z[Z > 0].ravel())
            vmin    = np.percentile(log_pos, pct_lo) if log_pos.size else 0.0
            vmax    = np.percentile(log_pos, pct_hi) if log_pos.size else 1.0
            C       = np.clip(np.log10(np.maximum(Z, 1e-300)), vmin, vmax)
            norm    = Normalize(vmin=vmin, vmax=vmax)
    return C, norm


def plot_manga_profiles(data, nangles, selected_bin=0,
                        highlight_color='lime', figsize=(16, 11), day=100):
    theta_edges   = np.linspace(0, np.pi, nangles + 1)
    theta_centres = 0.5 * (theta_edges[:-1] + theta_edges[1:])

    panel_cfg = [
        ('v_r', 'RdBu_r', r'$v_r$  [cm s$^{-1}$]', True),
        ('dens', 'inferno', r'Density  [g cm$^{-3}$]', False),
        ('temp', 'plasma', 'Temperature  [K]', False),
    ]

    r_vals  = data.loc[0, 'r'].values
    r_edges = radial_edges(r_vals)
    Theta, R = np.meshgrid(theta_edges, r_edges)

    arrays = {
        col: np.array([data.loc[a, col].values for a in range(nangles)])
        for col, *_ in panel_cfg
    }

    norms = {}
    for col, _, _, diverging in panel_cfg:
        _, norm = log_compress(arrays[col].ravel(), diverging)
        norms[col] = norm

    fig = plt.figure(figsize=figsize)
    polar_axes = [fig.add_subplot(2, 3, i + 1, projection='polar') for i in range(3)]
    line_axes  = [fig.add_subplot(2, 3, i + 4) for i in range(3)]

    th0 = theta_edges[selected_bin]
    th1 = theta_edges[selected_bin + 1]
    th_arc = np.linspace(th0, th1, 60)
    r_lo, r_hi = r_edges[0], r_edges[-1]

    for ax_p, ax_l, (col, cmap, ylabel, diverging) in zip(polar_axes, line_axes, panel_cfg):
        norm = norms[col]
        C, _ = log_compress(arrays[col], diverging)
        pcm  = ax_p.pcolormesh(Theta, R, C.T, cmap=cmap, norm=norm, shading='flat')

        cb = fig.colorbar(pcm, ax=ax_p, pad=0.12, fraction=0.046)
        cb.set_label(
            r'$v_r$  [cm s$^{-1}$]  (symlog)' if diverging
            else r'$\log_{10}$ ' + ylabel,
            fontsize=8,
        )

        ax_p.set_thetamin(0)
        ax_p.set_thetamax(180)
        ax_p.set_theta_zero_location('N')
        ax_p.set_theta_direction(-1)
        ax_p.set_thetagrids(
            np.degrees(theta_centres),
            labels=[str(i) for i in range(nangles)],
            fontsize=6.5,
        )
        ax_p.set_title(ylabel, pad=14, fontsize=10)
        ax_p.xaxis.grid(True, linestyle=':', linewidth=0.5, color='white', alpha=0.4)
        ax_p.yaxis.grid(True, linestyle=':', linewidth=0.5, color='white', alpha=0.4)

        for th in (th0, th1):
            ax_p.plot([th, th], [r_lo, r_hi], color=highlight_color, lw=2.0)
        for r_arc in (r_lo, r_hi):
            ax_p.plot(th_arc, np.full_like(th_arc, r_arc), color=highlight_color, lw=2.0)

        d_sel = data.loc[selected_bin]
        r_sel = d_sel['r'].values
        y_sel = d_sel[col].values

        ax_l.plot(r_sel, y_sel, color='C0', lw=1.5)
        ax_l.set_xlabel('Radius [cm]')
        ax_l.set_ylabel(ylabel)
        ax_l.set_title(f'{ylabel}  —  bin {selected_bin}', fontsize=9)
        ax_l.grid(True, alpha=0.3)

        if diverging:
            abs_nonzero = np.abs(y_sel[y_sel != 0])
            linthresh   = float(np.median(abs_nonzero)) if abs_nonzero.size else 1.0
            ax_l.set_yscale('symlog', linthresh=linthresh)
            ax_l.axhline(0, color='k', lw=0.8, ls='--')
        else:
            ax_l.set_yscale('log')

    fig.suptitle(
        f'MANGA merger model — day {day} (highlighted angle bin {selected_bin})',
        fontsize=12, y=1.01,
    )
    fig.tight_layout()
    return fig
