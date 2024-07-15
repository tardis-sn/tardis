

set -ex



test -f ${PREFIX}/lib/libopenblasp-r0.3.27.so
python -c "import ctypes; ctypes.cdll['${PREFIX}/lib/libopenblasp-r0.3.27.so']"
exit 0
