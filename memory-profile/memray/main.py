from tardis import run_tardis

sim = run_tardis("tardis_example.yml",
                 virtual_packet_logging=True,
                 show_convergence_plots=True,
                 export_convergence_plots=True,
                 log_level="INFO")

# memray memray run --trace-python-allocators main.py
# memray flamegraph --leaks memray-main.py.12683.bin