__version__ = '1.1.0'
try:
    from ._LRphase import longest  # noqa
except ImportError:
    def longest(args):
        """
        Args:
            args:
        """
        return max(args, key=len)
