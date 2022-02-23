import sys
import logging


# Custom log level
REPORT = 25


class CrashingHandler(logging.StreamHandler):
    def emit(self, record):
        """Unlike the method it overrides, this will not catch exceptions"""
        msg = self.format(record)
        stream = self.stream
        stream.write(msg)
        stream.write(self.terminator)
        self.flush()


class NiceFormatter(logging.Formatter):
    """
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).

    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno not in (logging.INFO, REPORT):
            record.msg = f"{record.levelname}: {record.msg}"
        return super().format(record)


def setup_logging(logger, log_to_stderr=True, minimal=False, quiet=False, debug=0):
    """
    Attach handler to the global logger object
    """
    # For --report=minimal, we need this custom log level because we want to
    # print nothing except the minimal report and therefore cannot use the
    # INFO level (and the ERROR level would give us an 'ERROR:' prefix).
    logging.addLevelName(REPORT, "REPORT")

    stream_handler = CrashingHandler(sys.stderr if log_to_stderr else sys.stdout)
    stream_handler.setFormatter(NiceFormatter())
    # debug overrides quiet overrides minimal
    if debug > 0:
        level = logging.DEBUG
    elif quiet:
        level = logging.ERROR
    elif minimal:
        level = REPORT
    else:
        level = logging.INFO
    stream_handler.setLevel(level)
    logger.setLevel(level)
    logger.addHandler(stream_handler)
