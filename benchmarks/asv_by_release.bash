#!/usr/bin/env bash

RELEASE_LIST=$(git tag -l "release-202[43]*" | sort -r)

readarray -t RELEASE_TAGS <<<"${RELEASE_LIST[@]}"
RELEASE_HASHES=()
for release_tag in "${RELEASE_TAGS[@]}"; do
  echo "Tag: ${release_tag}"
  HASH_COMMIT=$(git show-ref -s "${release_tag}")
  RELEASE_HASHES+=("${HASH_COMMIT}")
done
echo "RELEASE_HASHES: ${#RELEASE_HASHES[*]}"

ASV_CONFIG_PATH="/app/asv"
cd "${ASV_CONFIG_PATH}" || exit

rm -f release_hashes.txt
touch release_hashes.txt
for release_hash in "${RELEASE_HASHES[@]}"; do
  echo "${release_hash}" >>release_hashes.txt
done

function show_timed_time {
  local time=${1}
  local milliseconds="${time: -3}"
  local seconds=$((time / 1000))
  local minutes=0
  local minutes_display=""
  local hours=0
  local hours_display=""
  local days=0
  local days_display=""

  if [[ "${seconds}" -gt 59 ]]; then
    minutes=$((seconds / 60))
    seconds=$((seconds % 60))
    minutes_display="${minutes}m "
  fi

  if [[ "${minutes}" -gt 59 ]]; then
    hours=$((minutes / 60))
    minutes=$((minutes % 60))
    minutes_display="${minutes}m "
    hours_display="${hours}h "
  fi

  if [[ "${hours}" -gt 23 ]]; then
    days=$((hours / 24))
    hours=$((hours % 24))
    hours_display="${hours}h "
    days_display="${days}d "
  fi

  echo "${days_display}${hours_display}${minutes_display}${seconds}.${milliseconds}s"
}

start=$(date +%s%N | cut -b1-13)

#time asv run --quick --dry-run --bench energy_input_* HASHFILE:release_hashes.txt
#time asv run --bench energy_input_* release-2023.01.11..master
#time asv run --bench energy_input_* HASHFILE:release_hashes.txt

# time asv run --skip-existing-commits HASHFILE:release_hashes.txt

#time asv run \
#  --bench "time_create_energy_cdf([^_A-Za-z]|$)" \
#  --bench "time_move_packet([^_A-Za-z]|$)" \
#  --bench "time_pair_creation([^_A-Za-z]|$)" \
#  --bench "time_scatter_type([^_A-Za-z]|$)" \
#  --bench "time_get_perpendicular_vector([^_A-Za-z]|$)" \
#  --bench "time_klein_nishina([^_A-Za-z]|$)" \
#  --bench "time_get_with_config_item_string_item_access([^_A-Za-z]|$)" \
#  --bench "time_calculate_p_values([^_A-Za-z]|$)" \
#  --bench "time_trapezoid_integration([^_A-Za-z]|$)" \
#  --bench "time_get_inverse_doppler_factor([^_A-Za-z]|$)" \
#  --bench "time_get_inverse_doppler_factor_full_relativity([^_A-Za-z]|$)" \
#  release-2023.01.11..master

time asv run \
   --skip-existing-commits \
  --bench "BenchmarkMontecarloMontecarloNumbaOpacities" \
  ALL

end=$(date +%s%N | cut -b1-13)
runtime=$((end - start))
display_time="$(show_timed_time ${runtime})"
echo ""
echo "Time: ${display_time}"
echo ""
