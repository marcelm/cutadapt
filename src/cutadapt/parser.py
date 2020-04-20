"""
Parse adapter specifications
"""
import re
import logging
from typing import Type, Optional, List, Tuple, Iterator, Any, Dict
from xopen import xopen
from dnaio.readers import FastaReader
from .adapters import (
    Adapter, FrontAdapter, NonInternalFrontAdapter, BackAdapter, NonInternalBackAdapter,
    AnywhereAdapter, PrefixAdapter, SuffixAdapter, LinkedAdapter
)

logger = logging.getLogger(__name__)


class AdapterSpecification:
    """  # noqa: E501
    Description of a single non-linked adapter.

    These are the attributes:

    - name (None or str)
    - restriction (None, 'anchored', or 'noninternal')
    - sequence (nucleotide sequence as string)
    - parameters (dict with extra parameters such as 'max_error_rate', 'min_overlap')
    - cmdline_type ('front' for -a, 'back' for -g and 'anywhere' for -b)

    >>> AdapterSpecification.parse('a_name=ACGT;anywhere', 'back')
    AdapterSpecification(name='a_name', restriction=None, sequence='ACGT', parameters={'anywhere': True}, cmdline_type='back')
    """

    def __init__(
        self,
        name: str,
        restriction: Optional[str],
        sequence: str,
        parameters,
        cmdline_type: str,
    ):
        assert restriction in (None, "anchored", "noninternal")
        assert cmdline_type in ("front", "back", "anywhere")
        self.name = name
        self.restriction = restriction
        self.sequence = sequence
        self.parameters = parameters
        self.cmdline_type = cmdline_type

    @classmethod
    def parse(cls, spec: str, cmdline_type: str):
        """Factory for creating an instance from a string specification"""
        name, restriction, sequence, parameters = cls._parse(spec, cmdline_type)
        return cls(name, restriction, sequence, parameters, cmdline_type)

    def __repr__(self):
        return '{}(name={!r}, restriction={!r}, sequence={!r}, parameters={!r}, cmdline_type={!r})'.format(
            self.__class__.__name__, self.name, self.restriction, self.sequence, self.parameters, self.cmdline_type)

    def __eq__(self, other):
        return (
            self.name == other.name
            and self.restriction == other.restriction
            and self.sequence == other.sequence
            and self.parameters == other.parameters
            and self.cmdline_type == other.cmdline_type
        )

    @staticmethod
    def expand_braces(sequence: str) -> str:
        """
        Replace all occurrences of ``x{n}`` (where x is any character) with n
        occurrences of x. Raise ValueError if the expression cannot be parsed.

        >>> AdapterSpecification.expand_braces('TGA{5}CT')
        'TGAAAAACT'
        """
        # Simple DFA with four states, encoded in prev
        result = ''
        prev = None
        for s in re.split('([{}])', sequence):
            if s == '':
                continue
            if prev is None:
                if s == '{':
                    raise ValueError('"{" must be used after a character')
                if s == '}':
                    raise ValueError('"}" cannot be used here')
                prev = s
                result += s
            elif prev == '{':
                prev = int(s)
                if not 0 <= prev <= 10000:
                    raise ValueError('Value {} invalid'.format(prev))
            elif isinstance(prev, int):
                if s != '}':
                    raise ValueError('"}" expected')
                result = result[:-1] + result[-1] * prev
                prev = None
            else:
                if s != '{':
                    raise ValueError('Expected "{"')
                prev = '{'
        # Check if we are in a non-terminating state
        if isinstance(prev, int) or prev == '{':
            raise ValueError("Unterminated expression")
        return result

    @staticmethod
    def _extract_name(spec: str) -> Tuple[Optional[str], str]:
        """
        Parse an adapter specification given as 'name=adapt' into 'name' and 'adapt'.
        """
        fields = spec.split('=', 1)
        name = None  # type: Optional[str]
        if len(fields) > 1:
            name, spec = fields
            name = name.strip()
        spec = spec.strip()
        return name, spec

    allowed_parameters = {
        # abbreviations
        'e': 'max_error_rate',
        'error_rate': 'max_error_rate',
        'o': 'min_overlap',

        # allowed parameters
        'max_error_rate': None,
        'min_overlap': None,
        'anywhere': None,
        'required': None,
        'optional': None,  # If this is specified, 'required' will be set to False
    }

    @classmethod
    def _parse_parameters(cls, spec: str):
        """Parse key=value;key=value;key=value into a dict"""

        fields = spec.split(';')
        result = dict()  # type: Dict[str,Any]
        for field in fields:
            field = field.strip()
            if not field:
                continue
            key, equals, value = field.partition('=')  # type: (str, str, Any)
            if equals == '=' and value == '':
                raise ValueError('No value given')
            key = key.strip()
            if key not in cls.allowed_parameters:
                raise KeyError('Unknown parameter {}'.format(key))
            # unabbreviate
            while cls.allowed_parameters[key] is not None:
                key = cls.allowed_parameters[key]  # type: ignore
            value = value.strip()
            if value == '':
                value = True
            else:
                try:
                    value = int(value)
                except ValueError:
                    value = float(value)
            if key in result:
                raise KeyError('Key {} specified twice'.format(key))
            result[key] = value
        if 'optional' in result and 'required' in result:
            raise ValueError("'optional' and 'required' cannot be specified at the same time")
        if 'optional' in result:
            result['required'] = False
            del result['optional']
        return result

    @classmethod
    def _parse(cls, spec, cmdline_type):
        """
        Parse an adapter specification for a non-linked adapter (without '...')

        Allow:
        'back' and ADAPTER
        'back' and ADAPTERX
        'back' and ADAPTER$
        'front' and ADAPTER
        'front' and XADAPTER
        'front' and ^ADAPTER
        'anywhere' and ADAPTER
        """
        if cmdline_type not in ("front", "back", "anywhere"):
            raise ValueError("cmdline_type must be front, back or anywhere")
        error = ValueError(
            "You cannot use multiple placement restrictions for an adapter at the same time. "
            "Choose one of ^ADAPTER, ADAPTER$, XADAPTER or ADAPTERX")
        spec, middle, parameters_spec = spec.partition(';')
        name, spec = cls._extract_name(spec)
        spec = spec.strip()
        parameters = cls._parse_parameters(parameters_spec)
        spec = cls.expand_braces(spec)

        # Special case for adapters consisting of only X characters:
        # This needs to be supported for backwards-compatibilitity
        if len(spec.strip('X')) == 0:
            return name, None, spec, {}

        front_restriction = None
        if spec.startswith('^'):
            front_restriction = 'anchored'
            spec = spec[1:]
        if spec.upper().startswith('X'):
            if front_restriction is not None:
                raise error
            front_restriction = 'noninternal'
            spec = spec.lstrip('xX')

        back_restriction = None
        if spec.endswith('$'):
            back_restriction = 'anchored'
            spec = spec[:-1]
        if spec.upper().endswith('X'):
            if back_restriction is not None:
                raise error
            back_restriction = 'noninternal'
            spec = spec.rstrip('xX')

        n_placement_restrictions = int(bool(front_restriction)) + int(bool(back_restriction))
        if n_placement_restrictions > 1:
            raise error

        if cmdline_type == 'front' and back_restriction:
            raise ValueError(
                "Allowed placement restrictions for a 5' adapter are XADAPTER and ^ADAPTER")
        if cmdline_type == 'back' and front_restriction:
            raise ValueError(
                "Allowed placement restrictions for a 3' adapter are ADAPTERX and ADAPTER$")

        assert front_restriction is None or back_restriction is None
        if front_restriction is not None:
            restriction = front_restriction
        else:
            restriction = back_restriction

        if cmdline_type == 'anywhere' and restriction is not None:
            raise ValueError(
                "Placement restrictions (with X, ^, $) not supported for 'anywhere' (-b) adapters")

        return name, restriction, spec, parameters

    @staticmethod
    def _restriction_to_class(cmdline_type, restriction):
        if cmdline_type == 'front':
            if restriction is None:
                return FrontAdapter
            elif restriction == 'anchored':
                return PrefixAdapter
            elif restriction == 'noninternal':
                return NonInternalFrontAdapter
            else:
                raise ValueError(
                    'Value {} for a front restriction not allowed'.format(restriction))
        elif cmdline_type == 'back':
            if restriction is None:
                return BackAdapter
            elif restriction == 'anchored':
                return SuffixAdapter
            elif restriction == 'noninternal':
                return NonInternalBackAdapter
            else:
                raise ValueError(
                    'Value {} for a back restriction not allowed'.format(restriction))
        else:
            assert cmdline_type == 'anywhere'
            if restriction is None:
                return AnywhereAdapter
            else:
                raise ValueError('No placement may be specified for "anywhere" adapters')

    def adapter_class(self):
        return self._restriction_to_class(self.cmdline_type, self.restriction)


class AdapterParser:
    """
    Factory for Adapter classes that all use the same default parameters (error rate,
    indels etc.). The given **kwargs will be passed to the Adapter constructors.
    """
    def __init__(self, **kwargs):
        # kwargs: max_error_rate, min_overlap, read_wildcards, adapter_wildcards, indels
        self.default_parameters = kwargs

    def _parse(self, spec: str, cmdline_type: str = "back", name: Optional[str] = None) -> Adapter:
        """
        Parse an adapter specification not using ``file:`` notation and return
        an object of an appropriate Adapter class.

        name -- Adapter name if not included as part of the spec. (If spec is
        'name=ADAPTER', name will be 'name'.)

        cmdline_type -- describes which commandline parameter was used (``-a``
        is 'back', ``-b`` is 'anywhere', and ``-g`` is 'front').
        """
        if cmdline_type not in ('front', 'back', 'anywhere'):
            raise ValueError('cmdline_type cannot be {!r}'.format(cmdline_type))
        spec1, middle, spec2 = spec.partition('...')
        if middle == '...' and spec1 and spec2:
            return self._parse_linked(spec1, spec2, name, cmdline_type)

        if middle == '...':
            spec, cmdline_type = self._normalize_ellipsis(spec1, spec2, cmdline_type)
        else:
            spec = spec1
        return self._parse_not_linked(spec, name, cmdline_type)

    @staticmethod
    def _normalize_ellipsis(spec1: str, spec2: str, cmdline_type) -> Tuple[str, str]:
        if not spec1:
            if cmdline_type == 'back':
                # -a ...ADAPTER
                spec = spec2
            else:
                # -g ...ADAPTER
                raise ValueError('Invalid adapter specification')
        elif not spec2:
            if cmdline_type == 'back':
                # -a ADAPTER...
                cmdline_type = 'front'
                spec = spec1
            else:
                # -g ADAPTER...
                spec = spec1
        else:
            raise ValueError("Expected either spec1 or spec2")
        return spec, cmdline_type

    def _parse_not_linked(self, spec: str, name: Optional[str], cmdline_type: str) -> Adapter:
        aspec = AdapterSpecification.parse(spec, cmdline_type)
        adapter_class = aspec.adapter_class()  # type: Type[Adapter]
        if not name:
            name = aspec.name
        if aspec.parameters.pop('anywhere', False) and adapter_class in (FrontAdapter, BackAdapter):
            aspec.parameters['force_anywhere'] = True
        parameters = self.default_parameters.copy()
        parameters.update(aspec.parameters)

        return adapter_class(sequence=aspec.sequence, name=name, **parameters)

    def _parse_linked(self, spec1: str, spec2: str, name: Optional[str], cmdline_type: str) -> LinkedAdapter:
        """Return a linked adapter from two specification strings"""

        if cmdline_type == 'anywhere':
            raise ValueError("'anywhere' (-b) adapters may not be linked")
        front_spec = AdapterSpecification.parse(spec1, 'front')
        back_spec = AdapterSpecification.parse(spec2, 'back')
        if not name:
            name = front_spec.name

        if cmdline_type == 'back' and front_spec.restriction is None:
            import textwrap
            logger.warning('\n'.join(textwrap.wrap(
                "You specified a linked adapter as '-a ADAPTER1...ADAPTER2'. "
                "The interpretation of what this means has changed in Cutadapt 2.0. "
                "(The 5' adapter is now no longer anchored by default.) "
                "To get results consist with the old behavior, you need to anchor "
                "the 5' adapter explicitly as in '-a ^ADAPTER1...ADAPTER2'."
            )))

        front_anchored = front_spec.restriction is not None
        back_anchored = back_spec.restriction is not None

        front_parameters = self.default_parameters.copy()
        front_parameters.update(front_spec.parameters)
        back_parameters = self.default_parameters.copy()
        back_parameters.update(back_spec.parameters)

        if cmdline_type == 'front':
            # -g requires both adapters to be present
            front_required = True
            back_required = True
        else:
            # -a requires only the anchored adapters to be present
            front_required = front_anchored
            back_required = back_anchored

        # Handle parameters overriding whether an adapter is required
        front_required = front_parameters.pop('required', front_required)
        back_required = back_parameters.pop('required', back_required)

        front_adapter = front_spec.adapter_class()(front_spec.sequence, name=None,
            **front_parameters)
        back_adapter = back_spec.adapter_class()(back_spec.sequence, name=None,
            **back_parameters)

        return LinkedAdapter(
            front_adapter=front_adapter,
            back_adapter=back_adapter,
            front_required=front_required,
            back_required=back_required,
            name=name,
        )

    def parse(self, spec: str, cmdline_type: str = 'back') -> Iterator[Adapter]:
        """
        Parse an adapter specification and yield appropriate Adapter classes.
        This works like the _parse_no_file() function above, but also supports the
        ``file:`` notation for reading adapters from an external FASTA
        file. Since a file can contain multiple adapters, this
        function is a generator.
        """
        if spec.startswith('file:'):
            # read adapter sequences from a file
            with xopen(spec[5:], mode="rb", threads=0) as f:
                fasta = FastaReader(f)
                for record in fasta:
                    name = record.name.split(None, 1)
                    name = name[0] if name else None
                    yield self._parse(record.sequence, cmdline_type, name=name)
        else:
            yield self._parse(spec, cmdline_type, name=None)

    def parse_multi(self, type_spec_pairs: List[Tuple[str, str]]) -> List[Adapter]:
        """
        Parse all three types of commandline options that can be used to
        specify adapters. adapters must be a list of (str, str) pairs, where the first is
        the adapter type (either 'front', 'back' or 'anywhere') and the second is the
        adapter specification given on the commandline

        Return a list of appropriate Adapter classes.
        """
        adapters = []  # type: List[Adapter]
        for cmdline_type, spec in type_spec_pairs:
            if cmdline_type not in {'front', 'back', 'anywhere'}:
                raise ValueError('adapter type must be front, back or anywhere')
            adapters.extend(self.parse(spec, cmdline_type))
        return adapters
