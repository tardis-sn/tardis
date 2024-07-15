

set -ex



pytest -v tests/ -k "not test_nested_typevar_in_signature"
exit 0
