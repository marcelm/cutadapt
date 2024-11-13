import re
from dataclasses import dataclass
from typing import Iterator, Type


@dataclass
class Token:
    value: str

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}("{self.value}")'


class StringToken(Token):
    pass


class BraceToken(Token):
    pass


class TokenizeError(Exception):
    pass


def tokenize_braces(s: str, left: str = "{", right: str = "}") -> Iterator[Token]:
    """
    >>> list(tokenize_braces(""))
    []
    >>> list(tokenize_braces("before {braced} after"))
    [StringToken("before "), BraceToken("braced"), StringToken(" after")]
    >>> list(tokenize_braces("ab{cd}{ef}"))
    [StringToken("ab"), BraceToken("cd"), BraceToken("ef")]
    >>> list(tokenize_braces("ab(cd)ef", left="(", right=")"))
    [StringToken("ab"), BraceToken("cd"), StringToken("ef")]
    """
    if len(left) != 1 or len(right) != 1 or left == right:
        raise ValueError("left and right must be unequal one-character strings")
    for value in re.split(
        f"({re.escape(left)}[^{re.escape(right)}]*{re.escape(right)})", s
    ):
        if value == "":
            continue
        if value.startswith(left) and value.endswith(right):
            value = value[1:-1]
            token_class: Type[Token] = BraceToken
        else:
            token_class = StringToken
        if left in value:
            raise TokenizeError(f"Unexpected '{left}' encountered")
        if right in value:
            raise TokenizeError(f"Unexpected '{right}' encountered")
        yield token_class(value)
