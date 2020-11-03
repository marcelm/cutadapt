from textwrap import dedent
import pytest

from dnaio import Sequence
from cutadapt.adapters import LinkedAdapter, BackAdapter, FrontAdapter, InvalidCharacter
from cutadapt.parser import AdapterParser, AdapterSpecification
from cutadapt.modifiers import ModificationInfo


def test_expand_braces():
    expand_braces = AdapterSpecification.expand_braces
    assert expand_braces('') == ''
    assert expand_braces('A') == 'A'
    assert expand_braces('A{0}') == ''
    assert expand_braces('A{1}') == 'A'
    assert expand_braces('A{2}') == 'AA'
    assert expand_braces('A{2}C') == 'AAC'
    assert expand_braces('ACGTN{3}TGACCC') == 'ACGTNNNTGACCC'
    assert expand_braces('ACGTN{10}TGACCC') == 'ACGTNNNNNNNNNNTGACCC'
    assert expand_braces('ACGTN{3}TGA{4}CCC') == 'ACGTNNNTGAAAACCC'
    assert expand_braces('ACGTN{0}TGA{4}CCC') == 'ACGTTGAAAACCC'


def test_expand_braces_fail():
    for expression in ['{', '}', '{}', '{5', '{1}', 'A{-7}', 'A{', 'A{1', 'N{7', 'AN{7', 'A{4{}',
            'A{4}{3}', 'A{b}', 'A{6X}', 'A{X6}', 'A}A']:
        with pytest.raises(ValueError):
            AdapterSpecification.expand_braces(expression)


def test_parse_file_notation(tmpdir):
    tmp_path = str(tmpdir.join('adapters.fasta'))
    with open(tmp_path, 'w') as f:
        f.write(dedent(""">first_name
            ADAPTER1
            >second_name
            ADAPTER2
            """))
    parser = AdapterParser(
        max_errors=0.2, min_overlap=4, read_wildcards=False,
        adapter_wildcards=False, indels=False)

    adapters = list(parser.parse('file:' + tmp_path, cmdline_type='back'))
    assert len(adapters) == 2
    assert adapters[0].name == 'first_name'
    assert adapters[0].sequence == 'ADAPTER1'
    assert adapters[1].name == 'second_name'
    assert adapters[1].sequence == 'ADAPTER2'
    for a in adapters:
        assert a.max_error_rate == 0.2
        assert a.min_overlap == 4
        assert not a.read_wildcards
        assert not a.adapter_wildcards
        assert not a.indels


def test_parse_not_linked():
    p = AdapterSpecification.parse
    assert p('A', 'front') == AdapterSpecification(None, None, 'A', {}, 'front')
    assert p('A', 'back') == AdapterSpecification(None, None, 'A', {}, 'back')
    assert p('A', 'anywhere') == AdapterSpecification(None, None, 'A', {}, 'anywhere')
    assert p('^A', 'front') == AdapterSpecification(None, 'anchored', 'A', {}, 'front')
    assert p('XXXA', 'front') == AdapterSpecification(None, 'noninternal', 'A', {}, 'front')
    assert p('A$', 'back') == AdapterSpecification(None, 'anchored', 'A', {}, 'back')
    assert p('AXXXX', 'back') == AdapterSpecification(None, 'noninternal', 'A', {}, 'back')
    assert p('a_name=ADAPT', 'front') == AdapterSpecification('a_name', None, 'ADAPT', {}, 'front')


def test_parse_invalid_cmdline_type():
    with pytest.raises(ValueError) as e:
        AdapterSpecification._parse('A', 'invalid_type')
    assert "cmdline_type must be front, back or anywhere" in e.value.args[0]


@pytest.mark.parametrize("spec,cmdline_type", [
    ("^XA", "front"),
    ("^AX", "front"),
    ("XA$", "back"),
    ("AX$", "back"),
])
def test_parse_double_placement_restrictions(spec, cmdline_type):
    with pytest.raises(ValueError) as e:
        AdapterSpecification._parse(spec, cmdline_type)
    assert "cannot use multiple placement restrictions" in e.value.args[0]


def test_parse_misplaced_placement_restrictions():
    with pytest.raises(ValueError) as e:
        AdapterSpecification._parse("A$", "front")
    assert "Allowed placement restrictions for a 5' adapter" in e.value.args[0]
    with pytest.raises(ValueError) as e:
        AdapterSpecification._parse("^A", "back")
    assert "Allowed placement restrictions for a 3' adapter" in e.value.args[0]


def test_restriction_to_class():
    with pytest.raises(ValueError) as e:
        AdapterSpecification._restriction_to_class("anywhere", "noninternal")
    assert "No placement may be specified" in e.value.args[0]


def test_parse_parameters():
    p = AdapterSpecification._parse_parameters
    assert p('e=0.1') == {'max_errors': 0.1}
    assert p('error_rate=0.1') == {'max_errors': 0.1}
    assert p('max_errors=2') == {'max_errors': 2}
    assert p('o=5') == {'min_overlap': 5}
    assert p('min_overlap=5') == {'min_overlap': 5}
    assert p('o=7; e=0.4') == {'min_overlap': 7, 'max_errors': 0.4}
    assert p('anywhere') == {'anywhere': True}
    assert p('required') == {'required': True}
    assert p('optional') == {'required': False}

    with pytest.raises(ValueError):
        p('e=hallo')
    with pytest.raises(KeyError):
        p('bla=0.1')
    with pytest.raises(ValueError):
        p('e=')
    with pytest.raises(KeyError) as e:
        p('e=0.1;e=0.1')
    assert "specified twice" in e.value.args[0]
    with pytest.raises(KeyError) as e:
        p('e=0.1;max_errors=0.1')
    assert "specified twice" in e.value.args[0]
    with pytest.raises(ValueError) as e:
        p('optional; required')
    assert "cannot be specified at the same time" in e.value.args[0]


def test_parse_with_parameters(tmp_path):
    parser = AdapterParser(
        max_errors=0.2, min_overlap=4, read_wildcards=False,
        adapter_wildcards=False, indels=False)
    a = parser._parse('ACGTACGT; e=0.15', 'front')
    assert isinstance(a, FrontAdapter)
    assert a.max_error_rate == 0.15
    assert a.min_overlap == 4

    a = parser._parse('ACGTAAAA; o=5; e=0.11', 'back')
    assert isinstance(a, BackAdapter)
    assert a.max_error_rate == 0.11
    assert a.min_overlap == 5

    for spec in ('thename=ACG;e=0.15 ... TGT;e=0.17', 'thename=ACG;e=0.15...TGT;e=0.17'):
        a = parser._parse(spec, 'back')
        assert isinstance(a, LinkedAdapter)
        assert a.front_adapter.max_error_rate == 0.15
        assert a.back_adapter.max_error_rate == 0.17

    with pytest.raises(ValueError) as e:
        parser._parse("A", "invalid-cmdline-type")
    assert "cmdline_type cannot be" in e.value.args[0]


def test_parse_with_adapter_sequence_as_a_path(tmp_path):
    parser = AdapterParser()
    with pytest.raises(InvalidCharacter):
        parser._parse("invalid.character", "back")
    # user forgot to write "file:"
    path = (tmp_path / "afile.fasta")
    path.write_text(">abc\nACGT\n")
    with pytest.raises(InvalidCharacter) as e:
        list(parser.parse(str(path), "back"))
    assert "A file exists named" in e.value.args[0]


def test_parse_multi():
    parser = AdapterParser()
    with pytest.raises(ValueError) as e:
        parser.parse_multi([("invalid-type", "A")])
    assert "adapter type must be" in e.value.args[0]


def test_normalize_ellipsis():
    ne = AdapterParser._normalize_ellipsis
    assert ne("ACGT", "", "back") == ("ACGT", "front")  # -a ACGT...
    assert ne("ACGT", "", "front") == ("ACGT", "front")  # -g ACGT...
    assert ne("", "ACGT", "back") == ("ACGT", "back")  # -a ...ACGT
    with pytest.raises(ValueError) as e:
        # -g ...ACGT
        ne("", "ACGT", "front")
    assert "Invalid adapter specification" in e.value.args[0]

    with pytest.raises(ValueError) as e:
        ne("A", "C", "back")
    assert "either" in e.value.args[0]
    with pytest.raises(ValueError) as e:
        ne("A", "", "anywhere")
    assert "No ellipsis" in e.value.args[0]


@pytest.mark.parametrize("seq,req1,req2", [
    ("ACG...TGT", False, False),
    ("ACG...TGT$", False, True),
    ("^ACG...TGT", True, False),
    ("^ACG...TGT$", True, True),
])
def test_anchoring_makes_front_linked_adapter_required(seq, req1, req2):
    # -a X...Y
    a = AdapterParser()._parse(seq, "back")
    assert isinstance(a, LinkedAdapter)
    assert a.front_required is req1
    assert a.back_required is req2


@pytest.mark.parametrize("r1,r2,req1,req2", [
    ("", "", False, False),
    ("", ";required", False, True),
    (";required", "", True, False),
    (";required", ";required", True, True),
    ("", ";optional", False, False),
    (";optional", "", False, False),
    (";optional", ";optional", False, False),
])
def test_linked_adapter_back_required_optional(r1, r2, req1, req2):
    # -a X...Y
    a = AdapterParser()._parse("ACG" + r1 + "...TGT" + r2, "back")
    assert isinstance(a, LinkedAdapter)
    assert a.front_required is req1
    assert a.back_required is req2


@pytest.mark.parametrize("r1,r2,exp1,exp2", [
    ("", "", True, True),
    ("", ";required", True, True),
    (";required", "", True, True),
    (";required", ";required", True, True),
    ("", ";optional", True, False),
    (";optional", "", False, True),
    (";optional", ";optional", False, False),
])
def test_linked_adapter_front_required_optional(r1, r2, exp1, exp2):
    # -g X...Y
    a = AdapterParser()._parse("ACG" + r1 + "...TGT" + r2, "front")
    assert isinstance(a, LinkedAdapter)
    assert a.front_required is exp1
    assert a.back_required is exp2


def test_linked_adapter_parameters():
    # issue #394
    a = AdapterParser(max_errors=0.17, indels=False)._parse("ACG...TGT")
    assert isinstance(a, LinkedAdapter)
    assert a.front_adapter.max_error_rate == 0.17
    assert a.back_adapter.max_error_rate == 0.17
    assert not a.front_adapter.indels
    assert not a.back_adapter.indels


def test_linked_adapter_name():
    # issue #414
    a = AdapterParser()._parse("the_name=^ACG...TGT")
    assert isinstance(a, LinkedAdapter)
    assert a.create_statistics().name == "the_name"


def test_anywhere_parameter_back():
    parser = AdapterParser(max_errors=0.2, min_overlap=4, read_wildcards=False,
        adapter_wildcards=False, indels=True)
    adapter = list(parser.parse('CTGAAGTGAAGTACACGGTT;anywhere', 'back'))[0]
    assert isinstance(adapter, BackAdapter)
    assert adapter._force_anywhere

    # TODO move the rest to a separate test
    read = Sequence('foo1', 'TGAAGTACACGGTTAAAAAAAAAA')
    from cutadapt.modifiers import AdapterCutter
    cutter = AdapterCutter([adapter])
    trimmed_read = cutter(read, ModificationInfo(read))
    assert trimmed_read.sequence == ''


def test_anywhere_parameter_front():
    parser = AdapterParser(max_errors=0.2, min_overlap=4, read_wildcards=False,
        adapter_wildcards=False, indels=True)
    adapter = list(parser.parse('CTGAAGTGAAGTACACGGTT;anywhere', 'front'))[0]
    assert isinstance(adapter, FrontAdapter)
    assert adapter._force_anywhere

    # TODO move the rest to a separate test
    read = Sequence('foo1', 'AAAAAAAAAACTGAAGTGAA')
    from cutadapt.modifiers import AdapterCutter
    cutter = AdapterCutter([adapter])
    trimmed_read = cutter(read, ModificationInfo(read))
    assert trimmed_read.sequence == ''
