

set -ex



pip check
python -c "import contourpy as c; print(c.contour_generator(z=[[0, 1], [2, 3]]).lines(0.9))"
python -c "from contourpy.util import build_config; from pprint import pprint; pprint(build_config())"
python -c "import platform, sys; print(sys.version, platform.uname())"
exit 0
