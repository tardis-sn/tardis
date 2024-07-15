

set -ex



pip check
jupyter lab --version
jlpm --version
jupyter labextension list
jupyter lab licenses
jupyter lab path
jupyter server extension list
jupyter server extension list 1>server_extensions 2>&1
grep -iE "jupyterlab.*4\.2\.2.*OK" server_extensions
jupyter lab build --dev-build=False --minimize=False
jupyter lab clean
exit 0
