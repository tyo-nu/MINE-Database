import psutil

print(f"{psutil.Process().memory_info().rss / 1024 ** 2} MB")
print(f"n_cpds: {len(local_cpds)}\nn_rxns: {len(local_rxns)}")


# Get list of vars
locals()
zero_depth_bases = (str, bytes, Number, range, bytearray)
iteritems = 'items'

import sys
from numbers import Number
from collections import Set, Mapping, deque
def getsize(obj_0):
    """Recursively iterate to sum size of object & members."""
    _seen_ids = set()
    def inner(obj):
        obj_id = id(obj)
        if obj_id in _seen_ids:
            return 0
        _seen_ids.add(obj_id)
        size = sys.getsizeof(obj)
        if isinstance(obj, zero_depth_bases):
            pass # bypass remaining control flow and return
        elif isinstance(obj, (tuple, list, Set, deque)):
            size += sum(inner(i) for i in obj)
        elif isinstance(obj, Mapping) or hasattr(obj, iteritems):
            size += sum(inner(k) + inner(v) for k, v in getattr(obj, iteritems)())
        # Check for custom object instances - may subclass above too
        if hasattr(obj, '__dict__'):
            size += inner(vars(obj))
        if hasattr(obj, '__slots__'): # can have __slots__ with __dict__
            size += sum(inner(getattr(obj, s)) for s in obj.__slots__ if hasattr(obj, s))
        return size
    return inner(obj_0)


a = {i: getsize(i) for i in locals()}
print(a)
print(sum(a.values()))
