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


def tokenize_braces(s: str) -> Iterator[Token]:
    """
    >>> list(tokenize_braces(""))
    []
    >>> list(tokenize_braces("before {braced} after"))
    [StringToken("before "), BraceToken("braced"), StringToken(" after")]
    >>> list(tokenize_braces("ab{cd}{ef}"))
    [StringToken("ab"), BraceToken("cd"), BraceToken("ef")]
    """
    for value in re.split("({[^}]*})", s):
        if value == "":
            continue
        if value.startswith("{") and value.endswith("}"):
            value = value[1:-1]
            token_class: Type[Token] = BraceToken
        else:
            token_class = StringToken
        if "{" in value:
            raise TokenizeError("Unexpected '{' encountered")
        if "}" in value:
            raise TokenizeError("Unexpected '}' encountered")
        yield token_class(value)
