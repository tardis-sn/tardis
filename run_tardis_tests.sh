#!/bin/bash


LOG_FILE=/home/rishmitar/tardisDev/tardis_tests_$(date +%Y%m%d_%H%M%S).log
exec > >(tee -a "$LOG_FILE") 2>&1
 
echo "=== Test run started at $(date) ==="
echo "=== Logging to $LOG_FILE ==="


REGRESSION_DATA=/home/rishmitar/tardisDev/tardis-regression-data

pytest /home/rishmitar/tardisDev/tardis/tardis/opacities/tests/test_opacity_solver.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/opacities/tests/test_tau_sobolev.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/plasma/equilibrium/tests/test_coll_ionization_rates.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/plasma/equilibrium/tests/test_coll_ionization_strengths.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/plasma/equilibrium/tests/test_heating_cooling_rates.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/plasma/equilibrium/tests/test_ion_populations.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/plasma/equilibrium/tests/test_level_populations.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/plasma/equilibrium/tests/test_rate_matrix.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/plasma/tests/test_hdf_plasma.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/plasma/tests/test_tardis_model_density_config.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/simulation/tests/test_simulation.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/spectrum/formal_integral/tests/test_source_function.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/transport/montecarlo/packets/tests/test_rpacket_tracker.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/visualization/tests/test_plot_util.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/workflows/high_energy/tests/test_tardis_he_workflow.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

pytest /home/rishmitar/tardisDev/tardis/tardis/workflows/tests/test_workflows.py \
  --tardis-regression-data=$REGRESSION_DATA --generate-reference

echo "=== Test run completed at $(date) ==="
