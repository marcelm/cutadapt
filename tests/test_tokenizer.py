from cutadapt.tokenizer import tokenize_braces, StringToken, BraceToken


def test_tokenize_braces():
    tokenize = tokenize_braces
    assert list(tokenize("")) == []
    assert list(tokenize("text")) == [StringToken("text")]
    assert list(tokenize("before {variable} after")) == [
        StringToken("before "), BraceToken("variable"), StringToken(" after")
    ]
