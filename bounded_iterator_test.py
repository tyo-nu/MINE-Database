from threading import BoundedSemaphore
from collections import Iterator

class BoundedIterator(Iterator):
    def __init__(self, it, bound):
        self._it = it

        self._sem = BoundedSemaphore(bound)
    
    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self, timeout):
        if not self._sem.acquire(timeout=timeout):
            raise TimeoutError('Too many values un-acknowledged.')

        return next(self._it)

    def processed(self):
        self.sem.release()
