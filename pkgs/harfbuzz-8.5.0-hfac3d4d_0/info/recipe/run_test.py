import gi
gi.require_version('HarfBuzz', '0.0')
from gi.repository import HarfBuzz as hb
import sys

if hb.buffer_create () is None:
    sys.exit(1)
