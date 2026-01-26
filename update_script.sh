# Generate env.yml from conda-lock.yml for ASV benchmarks
# Run this whenever conda-lock.yml is updated

set -e

echo "Generating env.yml from conda-lock.yml"
conda-lock render --kind env --filename-template "env" -p linux-64

echo "env.yml generated successfully"
wc -l env.yml
