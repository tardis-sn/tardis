#!/usr/bin/env bash

set -e

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

ASV_CONFIG_PATH="/app/asv"

EXCLUDE_FILES=("__init__.py" "benchmark_base.py" "benchmark_template.py")

for x_file in *_x.py; do
  EXCLUDE_FILES+=("${x_file}")
done

BENCHMARK_FILES=()
for file in *.py; do
  IS_EXCLUDED_FILE=false
  for excluded_file in "${EXCLUDE_FILES[@]}"; do
    if [[ "${excluded_file}" == "${file}" ]]; then
      IS_EXCLUDED_FILE=true
      break
    fi
  done

  if "${IS_EXCLUDED_FILE}"; then
    continue
  fi

  BENCHMARK_FILES+=("${file}")
done

cd "${ASV_CONFIG_PATH}" || exit

declare -A TIME_RECORDS
BENCHMARK_INDEX=0
total_start=$(date +%s%N | cut -b1-13)
for benchmark_file in "${BENCHMARK_FILES[@]}"; do
  echo "Benchmark file: ${benchmark_file}"
  base_name=$(basename "${benchmark_file}")
  benchmark_name="${base_name%.py}"
  start=$(date +%s%N | cut -b1-13)
  time asv run --bench "${benchmark_name}"
  end=$(date +%s%N | cut -b1-13)
  runtime=$((end - start))
  display_time="$(show_timed_time ${runtime})"
  TIME_RECORDS[${BENCHMARK_INDEX}, 0]="${benchmark_file}"
  TIME_RECORDS[${BENCHMARK_INDEX}, 1]="${runtime}"
  TIME_RECORDS[${BENCHMARK_INDEX}, 2]="${display_time}"
  BENCHMARK_INDEX=$((BENCHMARK_INDEX + 1))
  echo "Time: ${display_time}"
  echo ""
  if [[ "energy_input_util" == "${benchmark_name}" ]]; then
    break
  fi
done
total_end=$(date +%s%N | cut -b1-13)
total_runtime=$((total_end - total_start))
display_time="$(show_timed_time ${total_runtime})"
echo "Total raw time: ${total_runtime} ms"
echo "Total time: ${display_time}"
echo ""

COLUMNS=3
ROWS=$((${#TIME_RECORDS[@]} / COLUMNS))

declare -A DUPLICATE_RECORDS
for ((i = 0; i < ROWS; i++)); do
  for ((j = 0; j < COLUMNS; j++)); do
    DUPLICATE_RECORDS[${i}, ${j}]=${TIME_RECORDS[${i}, ${j}]}
  done
done

declare -A SORTED_RECORDS
for i in $(seq 0 $((ROWS - 1))); do
  NEW_ROWS=$((${#DUPLICATE_RECORDS[@]} / COLUMNS))
  HIGHER_TIME[0]=""
  HIGHER_TIME[1]="-1"
  HIGHER_TIME[2]=""
  for j in $(seq 0 $((NEW_ROWS - 1))); do
    if [[ "${DUPLICATE_RECORDS[${j}, 1]}" -gt "${HIGHER_TIME[1]}" ]]; then
      HIGHER_TIME[0]="${DUPLICATE_RECORDS[${j}, 0]}"
      HIGHER_TIME[1]="${DUPLICATE_RECORDS[${j}, 1]}"
      HIGHER_TIME[2]="${DUPLICATE_RECORDS[${j}, 2]}"
      HIGHER_TIME[3]="${j}"
    fi
  done
  SORTED_RECORDS[${i}, 0]="${HIGHER_TIME[0]}"
  SORTED_RECORDS[${i}, 1]="${HIGHER_TIME[1]}"
  SORTED_RECORDS[${i}, 2]="${HIGHER_TIME[2]}"
  DUPLICATE_RECORDS[${HIGHER_TIME[3]}, 0]=""
  DUPLICATE_RECORDS[${HIGHER_TIME[3]}, 1]="-1"
  DUPLICATE_RECORDS[${HIGHER_TIME[3]}, 2]=""
done

echo "ORIGINAL RECORDS"
echo "================"
echo "Position: Benchmark file | Duration time | Raw time in ms"
for ((i = 0; i < ROWS; i++)); do
  echo "#$((i + 1)): ${TIME_RECORDS[${i}, 0]} | ${TIME_RECORDS[${i}, 2]} | ${TIME_RECORDS[${i}, 1]}"
done
echo "Total time: ${display_time}"
echo ""

echo "SLOWER RECORDS"
echo "=============="
echo "Position: Benchmark file | Duration time | Raw time in ms"
for ((i = 0; i < ROWS; i++)); do
  echo "#$((i + 1)): ${SORTED_RECORDS[${i}, 0]} | ${SORTED_RECORDS[${i}, 2]} | ${SORTED_RECORDS[${i}, 1]}"
done
echo "Total time: ${display_time}"
echo ""

echo "FASTER RECORDS"
echo "=============="
echo "Position: Benchmark file | Duration time | Raw time in ms"
for ((i = ROWS - 1; i >= 0; i--)); do
  echo "#$((ROWS - i)): ${SORTED_RECORDS[${i}, 0]} | ${SORTED_RECORDS[${i}, 2]} | ${SORTED_RECORDS[${i}, 1]}"
done
echo "Total time: ${display_time}"
echo ""
