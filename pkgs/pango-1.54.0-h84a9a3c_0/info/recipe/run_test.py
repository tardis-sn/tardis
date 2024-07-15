import gi
gi.require_version('PangoCairo', '1.0')
from gi.repository import PangoCairo
import sys

fontmap = PangoCairo.FontMap.get_default()
fontmap.list_families()
