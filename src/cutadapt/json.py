import json


class OneLine:
    """Wrap any value in this class to print it on one line in the JSON file"""

    def __init__(self, value):
        self.value = value


def dumps(obj, indent: int = 2, _level: int = 0) -> str:
    """
    Encode an object hierarchy as JSON string. In addition to
    what json.dumps in the standard library provides, this function
    allows disabling indentation for selected parts of the hierarchy
    by marking lists or dicts with the "OneLine" class.

    Arguments:
        obj: object to encode
        indent: indentation level

    >>> print(dumps({"a": [1, 2], "b": OneLine([3, 4]), "c": dict(x=5, y=6), "d": OneLine(dict(x=7, y=8))}))
    {
      "a": [
        1,
        2
      ],
      "b": [3, 4],
      "c": {
        "x": 5,
        "y": 6
      },
      "d": {"x": 7, "y": 8}
    }
    >>> print(dumps({"a": []}))
    {
      "a": []
    }
    """
    if isinstance(obj, (float, int, str, bool, OneLine)) or obj is None:
        if isinstance(obj, OneLine):
            obj = obj.value
        return json.dumps(obj)

    start = "\n" + (_level + 1) * indent * " "
    sep = "," + start
    end = "\n" + _level * indent * " "
    if isinstance(obj, (tuple, list)):
        if not obj:
            return "[]"
        return (
            "["
            + start
            + sep.join(dumps(elem, indent, _level + 1) for elem in obj)
            + end
            + "]"
        )
    elif isinstance(obj, dict):
        if not obj:
            return "{}"
        return (
            "{"
            + start
            + sep.join(
                json.dumps(k) + ": " + dumps(v, indent, _level + 1)
                for k, v in obj.items()
            )
            + end
            + "}"
        )
    else:
        raise ValueError(f"cannot serialize type {obj.__class__.__name__}")
