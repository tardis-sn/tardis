

set -ex



python -c "import pybtex.plugin; pybtex.plugin.find_plugin('pybtex.backends', 'docutils')()"
exit 0
