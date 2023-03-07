"""
Parse adapter specifications
"""
import re
import logging
from pathlib import Path
from typing import Type, Optional, List, Tuple, Any, Dict, Iterable
from xopen import xopen
from dnaio.readers import FastaReader
from .adapters import (
    Adapter,
    FrontAdapter,
    NonInternalFrontAdapter,
    BackAdapter,
    NonInternalBackAdapter,
    AnywhereAdapter,
    PrefixAdapter,
    SuffixAdapter,
    LinkedAdapter,
    InvalidCharacter,
    RightmostFrontAdapter,
)

logger = logging.getLogger(__name__)


def parse_search_parameters(spec: str):
    """Parse key=value;key=value;key=value into a dict"""

    allowed_parameters = {
        # abbreviations
        "e": "max_error_rate",
        "error_rate": "max_errors",
        "max_error_rate": "max_errors",
        "o": "min_overlap",
        # allowed parameters
        "max_errors": None,
        "min_overlap": None,
        "anywhere": None,
        "required": None,
        "optional": None,  # If this is specified, 'required' will be set to False
        "indels": None,
        "noindels": None,
        "rightmost": None,
    }

    fields = spec.split(";")
    result: Dict[str, Any] = dict()
    for field in fields:
        field = field.strip()
        if not field:
            continue
        key, equals, value = field.partition("=")  # type: (str, str, Any)
        key = key.strip()
        if key not in allowed_parameters:
            raise KeyError(f"Unknown parameter '{key}'")
        if equals == "=" and value == "":
            raise ValueError(f"No value given for key '{key}'")
        # unabbreviate
        while allowed_parameters[key] is not None:
            key = allowed_parameters[key]  # type: ignore
        value = value.strip()
        if value == "":
            value = True
        else:
            try:
                value = int(value)
            except ValueError:
                value = float(value)
        if key in result:
            raise KeyError(f"Key '{key}' specified twice")
        result[key] = value
    if "optional" in result and "required" in result:
        raise ValueError(
            "'optional' and 'required' cannot be specified at the same time"
        )
    if "indels" in result and "noindels" in result:
        raise ValueError("'indels' and 'noindels' cannot be specified at the same time")
    if "optional" in result:
        result["required"] = False
        del result["optional"]
    if "noindels" in result:
        result["indels"] = False
        del result["noindels"]
    return result


def expand_braces(sequence: str) -> str:
    """
    Replace all occurrences of ``x{n}`` (where x is any character) with n
    occurrences of x. Raise ValueError if the expression cannot be parsed.

    >>> expand_braces('TGA{5}CT')
    'TGAAAAACT'
    """
    # Simple DFA with four states, encoded in prev
    result = ""
    prev = None
    for s in re.split("([{}])", sequence):
        if s == "":
            continue
        if prev is None:
            if s == "{":
                raise ValueError('"{" must be used after a character')
            if s == "}":
                raise ValueError('"}" cannot be used here')
            prev = s
            result += s
        elif prev == "{":
            prev = int(s)
            if not 0 <= prev <= 10000:
                raise ValueError(f"Value {prev} invalid")
        elif isinstance(prev, int):
            if s != "}":
                raise ValueError('"}" expected')
            result = result[:-1] + result[-1] * prev
            prev = None
        else:
            if s != "{":
                raise ValueError('Expected "{"')
            prev = "{"
    # Check if we are in a non-terminating state
    if isinstance(prev, int) or prev == "{":
        raise ValueError("Unterminated expression")
    return result


def _normalize_ellipsis(spec1: str, spec2: str, adapter_type) -> Tuple[str, str]:
    if adapter_type == "anywhere":
        raise ValueError('No ellipsis ("...") allowed in "anywhere" adapters')
    if not spec1:
        if adapter_type == "back":
            # -a ...ADAPTER
            spec = spec2
        else:
            # -g ...ADAPTER
            raise ValueError("Invalid adapter specification")
    elif not spec2:
        if adapter_type == "back":
            # -a ADAPTER...
            adapter_type = "front"
            spec = spec1
        else:
            # -g ADAPTER...
            spec = spec1
    else:
        raise ValueError("Expected either spec1 or spec2")
    return spec, adapter_type


class AdapterSpecification:
    """# noqa: E501
    Description of a single non-linked adapter.

    These are the attributes:

    - name (None or str)
    - restriction (None, 'anchored', or 'noninternal')
    - sequence (nucleotide sequence as string)
    - search parameters (dict with keys such as 'max_errors', 'min_overlap')
    - adapter_type ('front' for -a, 'back' for -g and 'anywhere' for -b)

    >>> AdapterSpecification.parse('a_name=ACGT;anywhere', 'back')
    AdapterSpecification(name='a_name', restriction=None, sequence='ACGT', parameters={'anywhere': True}, adapter_type='back')
    """

    def __init__(
        self,
        name: Optional[str],
        restriction: Optional[str],
        sequence: str,
        parameters,
        adapter_type: str,
        rightmost: bool,
    ):
        assert restriction in (None, "anchored", "noninternal")
        assert adapter_type in ("front", "back", "anywhere")
        self.name = name
        self.restriction = restriction
        self.sequence = sequence
        self.parameters = parameters
        self.adapter_type = adapter_type
        self.rightmost = rightmost

    def __repr__(self):
        return "{}(name={!r}, restriction={!r}, sequence={!r}, parameters={!r}, adapter_type={!r})".format(
            self.__class__.__name__,
            self.name,
            self.restriction,
            self.sequence,
            self.parameters,
            self.adapter_type,
        )

    def __eq__(self, other):
        return (
            self.name == other.name
            and self.restriction == other.restriction
            and self.sequence == other.sequence
            and self.parameters == other.parameters
            and self.adapter_type == other.adapter_type
        )

    @staticmethod
    def _extract_name(spec: str) -> Tuple[Optional[str], str]:
        """
        Parse an adapter specification given as 'name=adapt' into 'name' and 'adapt'.
        """
        fields = spec.split("=", 1)
        name: Optional[str] = None
        if len(fields) > 1:
            name, spec = fields
            name = name.strip()
        spec = spec.strip()
        return name, spec

    @classmethod
    def parse(cls, spec: str, adapter_type: str) -> "AdapterSpecification":
        """
        Parse an adapter specification for a non-linked adapter (without '...')
        and return an AdapterSpecification instance.

        Allow:
        'back' and ADAPTER
        'back' and ADAPTERX
        'back' and ADAPTER$
        'front' and ADAPTER
        'front' and XADAPTER
        'front' and ^ADAPTER
        'anywhere' and ADAPTER
        """
        if adapter_type not in ("front", "back", "anywhere"):
            raise ValueError("adapter_type must be front, back or anywhere")

        spec, middle, parameters_spec = spec.partition(";")
        name, spec = cls._extract_name(spec)
        spec = spec.strip()
        parameters = parse_search_parameters(parameters_spec)
        spec = expand_braces(spec)
        rightmost = parameters.pop("rightmost", False)

        # Special case for adapters consisting of only X characters:
        # This needs to be supported for backwards-compatibilitity
        if len(spec.strip("X")) == 0:
            return cls(name, None, spec, {}, adapter_type, False)

        try:
            front_restriction, back_restriction, spec = cls._parse_restrictions(spec)
        except ValueError:
            raise ValueError(
                "You cannot use multiple placement restrictions for an adapter at the same time. "
                "Choose one of ^ADAPTER, ADAPTER$, XADAPTER or ADAPTERX"
            ) from None

        if adapter_type == "front" and back_restriction:
            raise ValueError(
                "Allowed placement restrictions for a 5' adapter are XADAPTER and ^ADAPTER"
            )
        if adapter_type == "back" and front_restriction:
            raise ValueError(
                "Allowed placement restrictions for a 3' adapter are ADAPTERX and ADAPTER$"
            )

        if front_restriction is not None:
            restriction: Optional[str] = front_restriction
        else:
            restriction = back_restriction

        if adapter_type == "anywhere" and restriction is not None:
            raise ValueError(
                "Placement restrictions (with X, ^, $) not supported for 'anywhere' (-b) adapters"
            )

        if "min_overlap" in parameters and restriction == "anchored":
            raise ValueError(
                "Setting 'min_overlap=' (or 'o=') for anchored adapters is not possible because "
                "anchored adapters always need to match in full."
            )

        if parameters.get("min_overlap", 0) > len(spec):
            raise ValueError(
                f"min_overlap={parameters['min_overlap']}"
                f" exceeds length of adapter {spec}"
            )

        if rightmost and (adapter_type != "front" or restriction is not None):
            raise ValueError("'rightmost' only allowed with regular 5' adapters")

        return cls(name, restriction, spec, parameters, adapter_type, rightmost)

    @staticmethod
    def _parse_restrictions(spec: str) -> Tuple[Optional[str], Optional[str], str]:
        front_restriction = None
        if spec.startswith("^"):
            front_restriction = "anchored"
            spec = spec[1:]
        if spec.upper().startswith("X"):
            if front_restriction is not None:
                raise ValueError("two front restrictions")
            front_restriction = "noninternal"
            spec = spec.lstrip("xX")

        back_restriction = None
        if spec.endswith("$"):
            back_restriction = "anchored"
            spec = spec[:-1]
        if spec.upper().endswith("X"):
            if back_restriction is not None:
                raise ValueError("two back restrictions")
            back_restriction = "noninternal"
            spec = spec.rstrip("xX")

        n_placement_restrictions = int(bool(front_restriction)) + int(
            bool(back_restriction)
        )
        if n_placement_restrictions > 1:
            raise ValueError("front and back restrictions")
        assert front_restriction is None or back_restriction is None
        return front_restriction, back_restriction, spec

    @staticmethod
    def _restriction_to_class(adapter_type, restriction, rightmost):
        """
        restriction: None, "anchored", or "noninternal"
        """
        if adapter_type == "front":
            if rightmost:
                assert restriction is None
                return RightmostFrontAdapter
            elif restriction is None:
                return FrontAdapter
            elif restriction == "anchored":
                return PrefixAdapter
            elif restriction == "noninternal":
                return NonInternalFrontAdapter
            else:
                raise ValueError(
                    f"Value {restriction} for a front restriction not allowed"
                )
        elif adapter_type == "back":
            if restriction is None:
                return BackAdapter
            elif restriction == "anchored":
                return SuffixAdapter
            elif restriction == "noninternal":
                return NonInternalBackAdapter
            else:
                raise ValueError(
                    f"Value {restriction} for a back restriction not allowed"
                )
        else:
            assert adapter_type == "anywhere"
            if restriction is None:
                return AnywhereAdapter
            else:
                raise ValueError(
                    'No placement may be specified for "anywhere" adapters'
                )

    def adapter_class(self):
        return self._restriction_to_class(
            self.adapter_type, self.restriction, self.rightmost
        )


def make_adapters_from_specifications(
    type_spec_pairs: List[Tuple[str, str]],
    search_parameters: Dict[str, Any],
) -> List[Adapter]:
    """
    Create a list of Adapter classes from specification strings and adapter types.

    type_spec_pairs -- a list of (str, str) pairs, where the first is
      the adapter type (either 'front', 'back' or 'anywhere') and the second is the
      adapter specification string, such as "ACGT;o=3" or "file:adapters.fasta"

    search_parameters -- A dict with default search parameters. These can be overriden by the
      adapter specifications. They are passed as **kwargs when instantiating the
      adapter classes.
      Possible keys: max_error_rate, min_overlap, read_wildcards, adapter_wildcards, indels

    Return a list of appropriate Adapter instances.
    """
    adapters: List[Adapter] = []
    for adapter_type, spec in type_spec_pairs:
        adapters.extend(
            make_adapters_from_one_specification(spec, adapter_type, search_parameters)
        )
    return adapters


def make_adapters_from_one_specification(
    spec: str,
    adapter_type: str,
    search_parameters: Dict[str, Any],
) -> Iterable[Adapter]:
    """
    Parse an adapter specification and yield appropriate Adapter classes.
    """
    if (
        spec.startswith("file:")
        or spec.startswith("^file:")
        or spec.startswith("file$:")
    ):
        anchoring_prefix = ""
        anchoring_suffix = ""
        if spec.startswith("^"):
            spec = spec[1:]
            anchoring_prefix = "^"
        elif spec.startswith("file$:"):
            spec = "file:" + spec[6:]
            anchoring_suffix = "$"
        path, _, parameters_spec = spec[5:].partition(";")
        parameters = search_parameters.copy()
        parameters.update(parse_search_parameters(parameters_spec))
        for name, spec in read_adapters_fasta(path):
            yield make_adapter(
                anchoring_prefix + spec + anchoring_suffix,
                adapter_type,
                parameters,
                name=name,
            )
    else:
        try:
            yield make_adapter(spec, adapter_type, search_parameters)
        except InvalidCharacter as e:
            if Path(spec).exists():
                extra_message = (
                    f"A file exists named '{spec}'. "
                    "To use the sequences in that file as adapter sequences, write 'file:' "
                    f"before the path, as in 'file:{spec}'."
                )
                raise InvalidCharacter(e.args[0] + "\n" + extra_message)
            else:
                raise


def make_adapter(
    spec: str,
    adapter_type: str,
    search_parameters: Dict[str, Any],
    name: Optional[str] = None,
) -> Adapter:
    """
    Parse an adapter specification not using ``file:`` notation and return
    an object of an appropriate Adapter class.

    name -- Adapter name if not included as part of the spec. (If spec is
    'name=ADAPTER', name will be 'name'.)

    adapter_type -- describes which commandline parameter was used (``-a``
    is 'back', ``-b`` is 'anywhere', and ``-g`` is 'front').

    search_parameters -- dict with default search parameters
    """
    if adapter_type not in ("front", "back", "anywhere"):
        raise ValueError("adapter_type must be front, back or anywhere")
    spec1, middle, spec2 = spec.partition("...")
    if middle == "..." and spec1 and spec2:
        return _make_linked_adapter(spec1, spec2, name, adapter_type, search_parameters)

    if middle == "...":
        spec, adapter_type = _normalize_ellipsis(spec1, spec2, adapter_type)
    else:
        spec = spec1
    return _make_not_linked_adapter(spec, name, adapter_type, search_parameters)


def _make_linked_adapter(
    spec1: str,
    spec2: str,
    name: Optional[str],
    adapter_type: str,
    search_parameters: Dict[str, Any],
) -> LinkedAdapter:
    """Return a linked adapter from two specification strings"""

    if adapter_type == "anywhere":
        raise ValueError("'anywhere' (-b) adapters may not be linked")
    front_spec = AdapterSpecification.parse(spec1, "front")
    back_spec = AdapterSpecification.parse(spec2, "back")
    if name is None:
        name = front_spec.name

    front_anchored = front_spec.restriction is not None
    back_anchored = back_spec.restriction is not None

    front_parameters = search_parameters.copy()
    front_parameters.update(front_spec.parameters)
    back_parameters = search_parameters.copy()
    back_parameters.update(back_spec.parameters)

    if adapter_type == "front":
        # -g requires both adapters to be present
        front_required = True
        back_required = True
    else:
        # -a requires only the anchored adapters to be present
        front_required = front_anchored
        back_required = back_anchored

    # Handle parameters overriding whether an adapter is required
    front_required = front_parameters.pop("required", front_required)
    back_required = back_parameters.pop("required", back_required)

    front_adapter = front_spec.adapter_class()(
        front_spec.sequence, name="linked_front", **front_parameters
    )
    back_adapter = back_spec.adapter_class()(
        back_spec.sequence, name="linked_back", **back_parameters
    )

    return LinkedAdapter(
        front_adapter=front_adapter,
        back_adapter=back_adapter,
        front_required=front_required,
        back_required=back_required,
        name=name,
    )


def _make_not_linked_adapter(
    spec: str,
    name: Optional[str],
    adapter_type: str,
    search_parameters: Dict[str, Any],
) -> Adapter:
    aspec = AdapterSpecification.parse(spec, adapter_type)
    adapter_class: Type[Adapter] = aspec.adapter_class()

    if aspec.parameters.pop("anywhere", False) and adapter_class in (
        FrontAdapter,
        BackAdapter,
        RightmostFrontAdapter,
    ):
        aspec.parameters["force_anywhere"] = True
    if "required" in aspec.parameters:
        raise ValueError(
            "'required' and 'optional' can only be used within linked adapters"
        )
    parameters = search_parameters.copy()
    parameters.update(aspec.parameters)
    return adapter_class(
        sequence=aspec.sequence,
        name=aspec.name if name is None else name,
        **parameters,
    )


def read_adapters_fasta(path):
    """
    Read adapter sequences from a FASTA file
    """
    with xopen(path, mode="rb", threads=0) as f:
        fasta = FastaReader(f)  # type: ignore
        for record in fasta:
            header = record.name.split(None, 1)
            name = header[0] if header else None
            yield name, record.sequence
