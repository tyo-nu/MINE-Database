"""Utils.py: contains basic functions reused in various contexts in other modules"""
__author__ = 'JGJeffryes'

def convert_sets_to_lists(obj):
    """Recursively converts dictionaries that contain sets to lists"""
    if isinstance(obj, set):
        try:
            obj = sorted(list(obj), key=lambda x: len(x))  # this brings short names to the top of the list
        except TypeError:
            obj = list(obj)
    elif isinstance(obj, dict):
        for key in obj:
            obj[key] = convert_sets_to_lists(obj[key])
    return obj
