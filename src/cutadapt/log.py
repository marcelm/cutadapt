import sys
import logging


# Custom log level
REPORT = 25


class NiceFormatter(logging.Formatter):
    """
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).

    Based on http://stackoverflow.com/a/9218261/715090 .
    """
    def format(self, record):
        if record.levelno not in (logging.INFO, REPORT):
            record.msg = '{}: {}'.format(record.levelname, record.msg)
        return super().format(record)


def setup_logging(logger, stdout=False, minimal=False, quiet=False, debug=False):
    """
    Attach handler to the global logger object
    """
    # For --report=minimal, we need this custom log level because we want to
    # print nothing except the minimal report and therefore cannot use the
    # INFO level (and the ERROR level would give us an 'ERROR:' prefix).
    logging.addLevelName(REPORT, 'REPORT')

    # Due to backwards compatibility, logging output is sent to standard output
    # instead of standard error if the -o option is used.
    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setFormatter(NiceFormatter())
    # debug overrides quiet overrides minimal
    if debug:
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
