

set -ex



pip check
jupyter server extension list
jupyter server extension list 1>server_extensions 2>&1
grep -iE "jupyter_lsp.*OK" server_extensions
exit 0
