# Licensed under a 3-clause BSD style license - see LICENSE.rst
import logging
logger = logging.getLogger('tardis')
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(console_formatter)
logger.addHandler(console_handler)