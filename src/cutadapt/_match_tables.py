import operator


def _acgt_table():
    """
    Return a translation table that maps A, C, G, T characters to the lower
    four bits of a byte. Other characters (including possibly IUPAC characters)
    are mapped to the most significant bit (0x80).

    Lowercase versions are also translated, and U is treated the same as T.
    """
    d = dict(A=1, C=2, G=4, T=8, U=8)
    t = bytearray([0x80]) * 256
    for c, v in d.items():
        t[ord(c)] = v
        t[ord(c.lower())] = v
    return bytes(t)


def _iupac_table():
    """
    Return a translation table for IUPAC characters.

    The table maps ASCII-encoded IUPAC nucleotide characters to bytes in which
    the four least significant bits are used to represent one nucleotide each.

    For the "N" wildcard, additionally the most significant bit is set (0x80),
    which allows it to match characters that are not A, C, G or T if _acgt_table
    was used to encode them.

    Whether two encoded characters x and y match can then be checked with the
    expression "x & y != 0".
    """
    A = 1
    C = 2
    G = 4
    T = 8
    iupac = dict(
        X=0,
        A=A,
        C=C,
        G=G,
        T=T,
        U=T,
        R=A | G,
        Y=C | T,
        S=G | C,
        W=A | T,
        K=G | T,
        M=A | C,
        B=C | G | T,
        D=A | G | T,
        H=A | C | T,
        V=A | C | G,
        N=A | C | G | T + 0x80,
    )
    t = bytearray(b"\0") * 256
    for c, v in iupac.items():
        t[ord(c)] = v
        t[ord(c.lower())] = v
    return bytes(t)


def _upper_table():
    table = bytes(range(256)).upper()
    return table


def all_matches_generator(ref: bytes, query: bytes, comp_op):
    for i, ref_char in enumerate(ref):
        matches = ""
        for j, query_char in enumerate(query):
            if j >= 128:  # Only ASCII characters supported.
                break
            if bool(comp_op(ref_char, query_char)):
                matches += chr(j)
        # NULL byte should not match anything
        yield matches.encode("ascii").replace(b"\00", b"")


def matches_lookup(ref_wildcards, query_wildcards):
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
