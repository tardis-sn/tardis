import logging
from colorlog import ColoredFormatter

LOG_LEVEL = logging.DEBUG
LOGFORMAT = "  %(log_color)s%(levelname)-8s%(reset)s | %(log_color)s%(message)s%(reset)s"

logging.root.setLevel(LOG_LEVEL)
formatter = ColoredFormatter(LOGFORMAT)
stream = logging.StreamHandler()
stream.setLevel(LOG_LEVEL)
stream.setFormatter(formatter)
log = logging.getLogger('tardis')
log.setLevel(LOG_LEVEL)
log.addHandler(stream)

log.debug("A message for developers")
log.info("Message for users who might want to know this")
log.warn("A warning message to inform the user")
log.error("A serious error for User")
log.critical("Critical Malfunctioning in Project")

