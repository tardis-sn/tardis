

set -ex



pip check
python -c "import pandas; pandas.test(extra_args=['-m not clipboard and not single_cpu and not slow and not network and not db', '-k', 'not (_not_a_real_test)', '--no-strict-data-files'])"
exit 0
