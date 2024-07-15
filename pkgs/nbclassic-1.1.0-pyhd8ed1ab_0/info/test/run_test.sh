

set -ex



pip check
jupyter nbclassic --help
jupyter nbclassic-extension --help
jupyter nbclassic-serverextension --help
jupyter nbclassic-bundlerextension --help
jupyter server extension list
jupyter server extension list 1>server_extensions 2>&1
cat server_extensions | grep -iE "nbclassic.*enabled"
cat server_extensions | grep -iE "nbclassic.*OK"
exit 0
