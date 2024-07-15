

set -ex



cd src && cd tests && pytest -vv --cov fastjsonschema -k "not (compile_to_code_ipv6_regex or compile_to_code_custom_format_with_refs)" --cov-report term-missing:skip-covered --no-cov-on-fail --cov-fail-under 84
exit 0
