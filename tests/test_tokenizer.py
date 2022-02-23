import pytest

from cutadapt.tokenizer import tokenize_braces, StringToken, BraceToken, TokenizeError


def test_tokenize_braces():
    tokenize = tokenize_braces
    assert list(tokenize("")) == []
    assert list(tokenize("text")) == [StringToken("text")]
    assert list(tokenize("before {variable} after")) == [
        StringToken("before "),
        BraceToken("variable"),
        StringToken(" after"),
    ]


def test_tokenize_unexpected_braces():
    with pytest.raises(TokenizeError):
        list(tokenize_braces("abc {def{ghi}"))

    with pytest.raises(TokenizeError):
        list(tokenize_braces("abc {def} gh} i"))
