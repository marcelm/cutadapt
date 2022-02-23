import pytest

from cutadapt.__main__ import (
    main,
    parse_cutoffs,
    parse_lengths,
    CommandLineError,
    setup_logging,
)


def test_help():
    with pytest.raises(SystemExit) as e:
        main(["--help"])
    assert e.value.args[0] == 0


def test_parse_cutoffs():
    assert parse_cutoffs("5") == (0, 5)
    assert parse_cutoffs("6,7") == (6, 7)
    with pytest.raises(CommandLineError):
        parse_cutoffs("a,7")
    with pytest.raises(CommandLineError):
        parse_cutoffs("a")
    with pytest.raises(CommandLineError):
        parse_cutoffs("a,7")
    with pytest.raises(CommandLineError):
        parse_cutoffs("1,2,3")


def test_parse_lengths():
    assert parse_lengths("25") == (25,)
    assert parse_lengths("17:25") == (17, 25)
    assert parse_lengths("25:") == (25, None)
    assert parse_lengths(":25") == (None, 25)
    with pytest.raises(CommandLineError):
        parse_lengths("1:2:3")
    with pytest.raises(CommandLineError):
        parse_lengths("a:2")
    with pytest.raises(CommandLineError):
        parse_lengths("a")
    with pytest.raises(CommandLineError):
        parse_lengths("2:a")
    with pytest.raises(CommandLineError):
        parse_lengths(":")


def test_setup_logging():
    import logging

    logger = logging.getLogger(__name__)
    setup_logging(logger, log_to_stderr=False, quiet=False, minimal=False, debug=False)
    logger.info("Log message")
    setup_logging(logger, log_to_stderr=False, debug=1)
    setup_logging(logger, log_to_stderr=False, quiet=True)
    setup_logging(logger, log_to_stderr=False, minimal=True)
