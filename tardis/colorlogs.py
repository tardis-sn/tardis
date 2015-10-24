import logging
import colorlog

# a theme is just a dict of strings to represent each level
THEME = {logging.CRITICAL: " [!!!!!] ",
         logging.ERROR:    "  [!!!]  ",
         logging.WARNING:  "   [!]   ",
         logging.INFO:     "    i    ",
         logging.DEBUG:    "   ...   "}
        

class Log:
    def __init__(self, lvl=logging.DEBUG, format=None):
        self._lvl = lvl
        if not format:
            self.format = "  %(log_color)s%(styledname)-8s%(reset)s | %(log_color)s%(message)s%(reset)s"
        logging.root.setLevel(self._lvl)
        self.formatter = colorlog.ColoredFormatter(self.format)
        self.stream = logging.StreamHandler()
        self.stream.setLevel(self._lvl)
        self.stream.setFormatter(self.formatter)
        self.logger = logging.getLogger('pythonConfig')
        self.logger.setLevel(self._lvl)
        self.logger.addHandler(self.stream)
        self.theme = THEME
        self.extra = {"styledname": self.theme[self._lvl]}

    def critical(self, message, *args, **kwargs):
        self.logger.critical(message, 
                             extra={"styledname": self.theme[logging.CRITICAL]}, 
                             *args, **kwargs)
    crit = c = fatal = critical
    def error(self, message, *args, **kwargs):
        self.logger.error(message, 
                          extra={"styledname": self.theme[logging.ERROR]}, 
                          *args, **kwargs)
    err = e = error
    def warn(self, message, *args, **kwargs):
        self.logger.warn(message, 
                         extra={"styledname": self.theme[logging.WARNING]}, 
                         *args, **kwargs)
    warning = w = warn
    def info(self, message, *args, **kwargs):
        self.logger.info(message, 
                         extra={"styledname": self.theme[logging.INFO]}, 
                         *args, **kwargs)
    inf = nfo = i = info
    def debug(self, message, *args, **kwargs):
        self.logger.debug(message, 
                          extra={"styledname": self.theme[logging.DEBUG]}, 
                          *args, **kwargs)
    dbg = d = debug

    # other convenience functions to set the global logging level
    def _parse_level(lvl):
        if lvl == logging.CRITICAL or lvl in ("critical", "crit", "c", "fatal"):
            return logging.CRITICAL
        elif lvl == logging.ERROR or lvl in ("error", "err", "e"):
            return logging.ERROR
        elif lvl == logging.WARNING or lvl in ("warning", "warn", "w"):
            return logging.WARNING
        elif lvl == logging.INFO or lvl in ("info", "inf", "nfo", "i"):
            return logging.INFO
        elif lvl == logging.DEBUG or lvl in ("debug", "dbg", "d"):
            return logging.DEBUG
        else:
            raise TypeError("Unrecognized logging level: %s" % lvl)
        
    def level(lvl=None):
        '''Get or set the logging level.'''
        if not lvl:
            return self._lvl
        self._lvl = self._parse_level(self._lvl)
        logging.root.setLevel(self._lvl)

log = Log()
