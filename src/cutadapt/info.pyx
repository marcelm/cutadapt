cdef class ModificationInfo:
    """
    An object of this class is created for each read that passes through the pipeline.
    Any information (except the read itself) that needs to be passed from one modifier
    to one later in the pipeline or from one modifier to the filters is recorded here.
    """
    cdef:
        public object matches
        public object original_read
        public object cut_prefix
        public object cut_suffix
        public object is_rc

    def __init__(self, read):
        self.matches = []
        self.original_read = read
        self.cut_prefix = None
        self.cut_suffix = None
        self.is_rc = None

    def __repr__(self):
        return (
            "ModificationInfo("
            f"matches={self.matches!r}, "
            f"original_read={self.original_read}, "
            f"cut_prefix={self.cut_prefix}, "
            f"cut_suffix={self.cut_suffix}, "
            f"is_rc={self.is_rc})"
        )
