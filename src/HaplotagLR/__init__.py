__version__ = '1.1.9'
try:
    from ._HaplotagLR import longest  # noqa
except ImportError:
    def longest(args):
        """
        Args:
            args:
        """
        return max(args, key=len)
