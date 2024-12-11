import operator
from collections.abc import Callable, Iterable


def _acgt_table() -> bytes:
    """
    Provide a translation table that maps A, C, G, T characters to the lower four bits of a byte.

    Other characters (including possibly IUPAC characters) are mapped to the most significant bit (0x80).

    Lowercase versions are also translated, and U is treated the same as T.

    Args:
        None

    Returns:
        bytes: A translation table that maps A, C, G, T characters to the lower
        four bits of a byte.
    """
    d = {"A": 1, "C": 2, "G": 4, "T": 8, "U": 8}
    t = bytearray([0x80]) * 256
    for c, v in d.items():
        t[ord(c)] = v
        t[ord(c.lower())] = v
    return bytes(t)


def _iupac_table() -> bytes:
    """
    Provides a translation table for IUPAC characters.

    The table maps ASCII-encoded IUPAC nucleotide characters to bytes in which
    the four least significant bits are used to represent one nucleotide each.

    For the "N" wildcard, additionally the most significant bit is set (0x80),
    which allows it to match characters that are not A, C, G or T if _acgt_table
    was used to encode them.

    Whether two encoded characters x and y match can then be checked with the
    expression "x & y != 0".

    Args:
        None

    Returns:
        bytes: A translation table for IUPAC characters
    """
    A: int = 1
    C: int = 2
    G: int = 4
    T: int = 8
    iupac = {
        "X": 0,
        "A": A,
        "C": C,
        "G": G,
        "T": T,
        "U": T,
        "R": A | G,
        "Y": C | T,
        "S": G | C,
        "W": A | T,
        "K": G | T,
        "M": A | C,
        "B": C | G | T,
        "D": A | G | T,
        "H": A | C | T,
        "V": A | C | G,
        "N": A | C | G | T + 0x80,
    }
    t = bytearray(b"\0") * 256
    for c, v in iupac.items():
        t[ord(c)] = v
        t[ord(c.lower())] = v
    return bytes(t)


def _upper_table() -> bytes:
    """"""
    return bytes(range(256)).upper()


def all_matches_generator(ref: bytes, query: bytes, comp_op: Callable) -> Iterable[bytes]:
    """
    Generates all possible pair matches between two translation tables.

    Args:
        ref (bytes): A reference translation table (ASCII-encoded characters mapped to bytes)
        query (bytes): A query translation table (ASCII-encoded characters mapped to bytes)
        comp_op (Callable): A comparison operator that takes two bytes and returns a boolean

    Yields:
        Iterable[bytes]: A generator that yields all possible pair matches between two translation tables
    """
    for _, ref_char in enumerate(ref):
        matches = ""
        for j, query_char in enumerate(query):
            if j >= 128:  # Only ASCII characters supported.
                break
            if comp_op(ref_char, query_char):
                matches += chr(j)
        # NULL byte should not match anything
        yield matches.encode("ascii").replace(b"\00", b"")


def matches_lookup(ref_wildcards: bool, query_wildcards: bool) -> list[bytes]:
    """
    Based on wildcard settings, return a list of all possible matches between two translation tables

    Args:
        ref_wildcards (bool): Whether the reference table can contain wildcards
        query_wildcards (bool): Whether the query table can contain wildcards

    Returns:
        list[bytes]: A list of all possible matches between two translation tables
    """
    if (not ref_wildcards) and (not query_wildcards):
        ref_table = _upper_table()
        query_table = _upper_table()
        comp_op = operator.eq
    elif ref_wildcards and (not query_wildcards):
        ref_table = _iupac_table()
        query_table = _acgt_table()
        comp_op = operator.and_
    elif (not ref_wildcards) and query_wildcards:
        ref_table = _acgt_table()
        query_table = _iupac_table()
        comp_op = operator.and_
    else:
        ref_table = _iupac_table()
        query_table = _iupac_table()
        comp_op = operator.and_
    return list(all_matches_generator(ref_table, query_table, comp_op))
