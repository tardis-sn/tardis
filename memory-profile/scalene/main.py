from tardis import run_tardis

sim = run_tardis("tardis_example.yml",
                 virtual_packet_logging=True,
                 show_convergence_plots=True,
                 export_convergence_plots=True,
                 log_level="INFO")

# scalene --memory --reduced-profile main.py
