#!/usr/bin/env python
# coding: utf-8

# # CylinderSim -- End-to-End Notebook
# 
# Set `MODE` and `PARAMS` in cell 1, then **Run All**.
# 
# | MODE | What happens |
# |------|--------------|
# | `"run"` | Compile C++, run simulation, save CSVs |
# | `"graph"` | Load existing CSVs, produce and save all Plotly graphs |
# | `"both"` | Full pipeline: run then graph |
# 
# Monitor progress while running:
# ```bash
# tail -f Simulations/051425cyl/progress.txt
# ```
# 
# **One-time conda setup:**
# ```bash
# conda install -c conda-forge eigen plotly
# 
# ```

# ## 1 -- Parameters & Mode
# Edit here. Tag this cell `parameters` to use with papermill.

# In[ ]:


import pathlib, os

# ------------------------------------------------------------
# MODE: "run" | "graph" | "both"
# ------------------------------------------------------------
MODE = "both"

PARAMS = {
    "dt":      1e-9,
    "datadt":  1000,
    "steps":   1_000_000_000,
    "temp":    298.15,
    "rho":     997.0,
    "mu":      8.8891e-4,
    "d":       5e-8,
    "l":       5e-7,
    "dsize":   10,
    "nbins":   50,
    "mode":    "save",
    "dir":     "Simulations",
    "res":     1000,
    "name":    "041426",
}

NM = 1e9
NS = 1e9
SIM_FOLDER = pathlib.Path(PARAMS['dir']) / (PARAMS['name'] + 'cyl')
print(f'MODE      : {MODE}')
print(f'Output dir: {SIM_FOLDER}')


# ## 2 -- Write & Compile `main.cpp`
# *Skipped in `graph` mode.*

# In[ ]:


if MODE in ('run', 'both'):
    import subprocess, sys

    main_cpp = r"""
#include "CylinderSim.hpp"
#include <chrono>

int main() {{
    auto t0 = std::chrono::high_resolution_clock::now();

    CylinderSim sim(
        {dt}, {datadt}, {steps},
        {temp}, {rho}, {mu},
        {d}, {l},
        {dsize}, {nbins},
        "{mode}", "{dir}", {res}, "{name}"
    );

    sim.updateProgress("Simulation started");
    sim.updateProgress("Parameters:\n" + sim.toString());
    sim.generateData();
    sim.saveSimCSV();

    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    sim.updateProgress("Elapsed: " + std::to_string(elapsed) + " s");
    return 0;
}}
""".format(**PARAMS)

    with open('main.cpp', 'w') as f:
        f.write(main_cpp)
    print('main.cpp written')

    eigen_path = os.path.join(os.environ.get('CONDA_PREFIX', '/usr'), 'include', 'eigen3')
    print(f'Eigen path: {eigen_path}')

    result = subprocess.run(
        ['g++', '-std=c++17', '-O2', f'-I{eigen_path}', 'main.cpp', '-o', 'cylinder_sim'],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f'Compile failed:{result.stderr}')
    print('Compiled successfully.')
else:
    print(f'MODE={MODE}: skipping compile.')


# ## 3 -- Run Simulation
# *Skipped in `graph` mode.*
# 
# This cell blocks until the simulation finishes. Monitor `progress.txt` in a terminal:
# ```bash
# tail -f Simulations/051425cyl/progress.txt
# ```

# In[ ]:


if MODE in ('run', 'both'):
    import subprocess, time

    print('Running simulation...')
    print(f'Progress log: {SIM_FOLDER / "progress.txt"}')

    t0 = time.time()
    result = subprocess.run(['./cylinder_sim'], capture_output=True, text=True)
    elapsed = time.time() - t0

    if result.returncode != 0:
        raise RuntimeError(f'Simulation failed:\n{result.stderr}')

    print(f'Simulation finished in {elapsed:.1f} s')
    print()
    progress_path = SIM_FOLDER / 'progress.txt'
    if progress_path.exists():
        print('--- progress.txt ---')
        print(progress_path.read_text())
else:
    print(f'MODE={MODE}: skipping simulation run.')


# ## 4 -- Load CSV Data
# *Skipped in `run` mode.*

# In[ ]:


if MODE in ('graph', 'both'):
    import pandas as pd
    import numpy as np
    from math import log, pi

    # Read parameters.txt -- single source of truth
    params_file = SIM_FOLDER / 'parameters.txt'
    raw = {}
    for line in params_file.read_text().splitlines():
        if '=' in line and not line.startswith('#'):
            k, _, v = line.partition('=')
            raw[k.strip()] = v.strip()

    def pD(k): return float(raw[k])
    def pI(k): return int(raw[k])

    DT        = pD('dt')
    DATADT    = pI('datadt')
    STEPS_RAW = pI('total_raw_steps')
    STEPS     = pI('steps(saved)')
    RES       = pI('res')
    L_M       = pD('l')
    D_M       = pD('d')
    TEMP      = pD('temp')
    RHO       = pD('rho')
    MU        = pD('mu')
    NBINS     = PARAMS['nbins']

    pos     = pd.read_csv(SIM_FOLDER / 'positions.csv')           * NM
    vel     = pd.read_csv(SIM_FOLDER / 'velocities.csv')
    thz     = pd.read_csv(SIM_FOLDER / 'thetaz.csv')
    omg     = pd.read_csv(SIM_FOLDER / 'omegas.csv')
    omg_obj = pd.read_csv(SIM_FOLDER / 'omegas_objframe.csv')
    ori     = pd.read_csv(SIM_FOLDER / 'orientations.csv')        * NM
    sqd     = pd.read_csv(SIM_FOLDER / 'squared_displacement.csv')

    n_rows   = len(pos)
    times_ns = np.arange(n_rows) * DT * DATADT * NS
    totalt   = times_ns[-1]
    n_res    = max(1, STEPS // RES)

    # Analytical diffusion reference lines
    kb      = 1.380649e-23
    vol     = pi * D_M**2 / 4 * L_M
    m       = vol * RHO
    iplanar = RHO * pi * L_M * D_M**2 * (L_M**2/3 + D_M**2/4) / 16
    gamma   = 2 * pi * MU * L_M
    znormal = 2*gamma / ((log(L_M/D_M) + 0.84) * m)
    D_trans = (kb * TEMP / m) / znormal
    expmsd  = D_trans * 2 * DT
    zrotnorm= gamma * L_M**2 / (iplanar * (log(L_M/D_M) - 0.66))
    Dr_rot  = (kb * TEMP / iplanar) / zrotnorm
    expmsad = Dr_rot * 2 * DT

    print(f'Loaded {n_rows} rows from {SIM_FOLDER}')
    print(f'Time range: 0 -- {totalt:.1f} ns')
else:
    print(f'MODE={MODE}: skipping data load.')


# ## 5 -- Graphs
# *Skipped in `run` mode.* All figures saved as interactive HTML.

# In[ ]:


if MODE in ('graph', 'both'):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    SIM_FOLDER.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Helper: save + show
    # ------------------------------------------------------------------
    def save_show(fig, filename):
        path = SIM_FOLDER / filename
        fig.write_html(str(path))
        fig.show()
        print(f'Saved {path}')

    # ------------------------------------------------------------------
    # graphxv: position & velocity histograms
    # ------------------------------------------------------------------
    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=[
            'x position (nm)', 'y position (nm)', 'z position (nm)',
            'x velocity (m/s)', 'y velocity (m/s)', 'z velocity (m/s)'
        ]
    )
    for j, col in enumerate(['x', 'y', 'z']):
        counts, bins = np.histogram(pos[col], bins=NBINS)
        fig.add_trace(go.Bar(x=bins[:-1], y=counts, marker_color='steelblue',
                             showlegend=False), row=1, col=j+1)
        fig.add_hline(y=counts.mean(), line_dash='dash', line_color='red',
                      row=1, col=j+1)
        v = vel[col]
        vcounts, vbins = np.histogram(v, bins=NBINS)
        bw = np.diff(vbins)[0]
        gauss = (len(v)*bw / (v.std()*(2*pi)**0.5)
                 * np.exp(-0.5*((vbins[:-1]-v.mean())/v.std())**2))
        fig.add_trace(go.Bar(x=vbins[:-1], y=vcounts, marker_color='steelblue',
                             showlegend=False), row=2, col=j+1)
        fig.add_trace(go.Scatter(x=vbins[:-1], y=gauss, mode='lines',
                                 line=dict(color='red'), showlegend=False),
                      row=2, col=j+1)
    fig.update_layout(title='Position & Velocity Distributions', height=700)
    save_show(fig, 'xvgraph.html')

    # ------------------------------------------------------------------
    # graphthw: angular position & velocity histograms
    # ------------------------------------------------------------------
    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=[
            'yaw angle (rad)', 'z-component of orientation', '',
            'roll omega (rad/s)', 'pitch omega (rad/s)', 'yaw omega (rad/s)'
        ]
    )
    for j, col in enumerate(['theta', 'z']):
        counts, bins = np.histogram(thz[col], bins=NBINS)
        fig.add_trace(go.Bar(x=bins[:-1], y=counts, marker_color='steelblue',
                             showlegend=False), row=1, col=j+1)
        fig.add_hline(y=counts.mean(), line_dash='dash', line_color='red',
                      row=1, col=j+1)
    for j, col in enumerate(['roll', 'pitch', 'yaw']):
        v = omg[col]
        vcounts, vbins = np.histogram(v, bins=NBINS)
        bw = np.diff(vbins)[0]
        gauss = (len(v)*bw / (v.std()*(2*pi)**0.5)
                 * np.exp(-0.5*((vbins[:-1]-v.mean())/v.std())**2))
        fig.add_trace(go.Bar(x=vbins[:-1], y=vcounts, marker_color='steelblue',
                             showlegend=False), row=2, col=j+1)
        fig.add_trace(go.Scatter(x=vbins[:-1], y=gauss, mode='lines',
                                 line=dict(color='red'), showlegend=False),
                      row=2, col=j+1)
    fig.update_layout(title='Angular Position & Velocity Distributions', height=700)
    save_show(fig, 'thwgraph.html')

    # ------------------------------------------------------------------
    # graph_rt_msd: translational MSD + position vs time
    # ------------------------------------------------------------------
    fig = make_subplots(rows=2, cols=1, subplot_titles=[
        'Squared displacement distribution', 'Position vs time'
    ])
    msd = sqd['Translational']
    counts, bins = np.histogram(msd, bins=NBINS)
    threshold = STEPS * 5e-4
    counts_t = np.where(counts < threshold, 0, counts)
    last = next((i for i, c in enumerate(counts_t) if c == 0), len(counts_t))
    fig.add_trace(go.Bar(x=bins[:last], y=counts_t[:last],
                         marker_color='steelblue', showlegend=False), row=1, col=1)
    fig.add_vline(x=msd.mean(), line_color='blue',
                  annotation_text=f'Mean={msd.mean():.3e}', row=1, col=1)
    fig.add_vline(x=expmsd, line_color='orange', line_dash='dash',
                  annotation_text=f'2Ddt={expmsd:.3e}', row=1, col=1)
    for col, color in zip(['x','y','z'], ['blue','red','green']):
        fig.add_trace(go.Scatter(x=times_ns, y=pos[col], mode='markers',
                                 marker=dict(size=2, color=color),
                                 name=f'{col} (nm)'), row=2, col=1)
    fig.update_xaxes(title_text='Displacement^2 (m^2)', row=1, col=1)
    fig.update_xaxes(title_text='Time (ns)', row=2, col=1)
    fig.update_yaxes(title_text='Frequency', row=1, col=1)
    fig.update_yaxes(title_text='Position (nm)', row=2, col=1)
    fig.update_layout(title='Translational MSD & Position vs Time', height=800)
    save_show(fig, 'rt_msd.html')

    # ------------------------------------------------------------------
    # graph_tht_msad: rotational MSD + angle vs time
    # ------------------------------------------------------------------
    fig = make_subplots(rows=2, cols=1, subplot_titles=[
        'Squared angular displacement distribution', 'Angle vs time'
    ])
    msad = sqd['Angular']
    counts, bins = np.histogram(msad, bins=NBINS)
    counts_t = np.where(counts < threshold, 0, counts)
    last = next((i for i, c in enumerate(counts_t) if c == 0), len(counts_t))
    fig.add_trace(go.Bar(x=bins[:last], y=counts_t[:last],
                         marker_color='steelblue', showlegend=False), row=1, col=1)
    fig.add_vline(x=msad.mean(), line_color='blue',
                  annotation_text=f'Mean={msad.mean():.3e}', row=1, col=1)
    fig.add_vline(x=expmsad, line_color='orange', line_dash='dash',
                  annotation_text=f'2Drdt={expmsad:.3e}', row=1, col=1)
    fig.add_trace(go.Scatter(x=times_ns, y=thz['theta'], mode='markers',
                             marker=dict(size=2, color='blue'),
                             name='yaw (rad)'), row=2, col=1)
    fig.add_trace(go.Scatter(x=times_ns, y=thz['z'], mode='markers',
                             marker=dict(size=2, color='red'),
                             name='z-component'), row=2, col=1)
    fig.update_xaxes(title_text='Angular displacement^2 (rad^2)', row=1, col=1)
    fig.update_xaxes(title_text='Time (ns)', row=2, col=1)
    fig.update_yaxes(title_text='Frequency', row=1, col=1)
    fig.update_yaxes(title_text='Angle (rad)', row=2, col=1)
    fig.update_layout(title='Rotational MSD & Angle vs Time', height=800)
    save_show(fig, 'tht_msad.html')

    # ------------------------------------------------------------------
    # graphPositions: 3-D position scatter
    # ------------------------------------------------------------------
    fig = go.Figure(go.Scatter3d(
        x=pos['x'].values[::n_res],
        y=pos['y'].values[::n_res],
        z=pos['z'].values[::n_res],
        mode='markers',
        marker=dict(size=2, color=times_ns[::n_res],
                    colorscale='Jet', colorbar=dict(title='Time (ns)'),
                    showscale=True)
    ))
    boxsize = float(pos.abs().max().max()) * 2
    fig.update_layout(
        title='Particle Brownian Motion -- 3D Positions',
        scene=dict(
            xaxis_title='x (nm)', yaxis_title='y (nm)', zaxis_title='z (nm)',
            xaxis=dict(range=[-boxsize/2, boxsize/2]),
            yaxis=dict(range=[-boxsize/2, boxsize/2]),
            zaxis=dict(range=[-boxsize/2, boxsize/2]),
            aspectmode='cube'
        ), height=700
    )
    save_show(fig, 'positions_3d.html')

    # ------------------------------------------------------------------
    # graphRod: 3-D rod orientation cones
    # ------------------------------------------------------------------
    lx = (pos['x'] - ori['x']).values[::n_res]
    ly = (pos['y'] - ori['y']).values[::n_res]
    lz = (pos['z'] - ori['z']).values[::n_res]
    ux = (2 * ori['x']).values[::n_res]
    uy = (2 * ori['y']).values[::n_res]
    uz = (2 * ori['z']).values[::n_res]
    t  = times_ns[::n_res]
    fig = go.Figure(go.Cone(
        x=lx, y=ly, z=lz, u=ux, v=uy, w=uz,
        colorscale='Jet', cmin=float(t.min()), cmax=float(t.max()),
        colorbar=dict(title='Time (ns)'),
        sizemode='absolute', sizeref=float(L_M*NM)*0.8,
        showscale=True
    ))
    l_ax = float(max(
        abs(lx).max(), abs(ly).max(), abs(lz).max(),
        abs(lx+ux).max(), abs(ly+uy).max(), abs(lz+uz).max()
    )) * 1.1
    fig.update_layout(
        title='Particle Orientation -- 3D Rod',
        scene=dict(
            xaxis_title='x (nm)', yaxis_title='y (nm)', zaxis_title='z (nm)',
            xaxis=dict(range=[-l_ax, l_ax]),
            yaxis=dict(range=[-l_ax, l_ax]),
            zaxis=dict(range=[-l_ax, l_ax]),
            aspectmode='cube'
        ), height=700
    )
    save_show(fig, 'rod_3d.html')

    # ------------------------------------------------------------------
    # graphOrientations: 3-D orientation vector tips
    # ------------------------------------------------------------------
    l_nm = L_M * NM
    fig = go.Figure(go.Scatter3d(
        x=ori['x'].values[::n_res],
        y=ori['y'].values[::n_res],
        z=ori['z'].values[::n_res],
        mode='markers',
        marker=dict(size=2, color=times_ns[::n_res],
                    colorscale='Jet', colorbar=dict(title='Time (ns)'),
                    showscale=True)
    ))
    fig.update_layout(
        title='Orientation Vector Tips -- 3D',
        scene=dict(
            xaxis_title='x (nm)', yaxis_title='y (nm)', zaxis_title='z (nm)',
            xaxis=dict(range=[-l_nm/2, l_nm/2]),
            yaxis=dict(range=[-l_nm/2, l_nm/2]),
            zaxis=dict(range=[-l_nm/2, l_nm/2]),
            aspectmode='cube'
        ), height=700
    )
    save_show(fig, 'orientations_3d.html')

    print('All graphs saved to', SIM_FOLDER)
else:
    print(f'MODE={MODE}: skipping graphs.')

