from collections import MutableMapping
import time


class TimeoutDict(MutableMapping):
    def _expire_check(self):
        call_time = time.time()
        for key, exp_time in list(self._expires.items()):
            if call_time >= exp_time:
                del self._data[key]
                del self._expires[key]

    __slots__ = ('_data', '_expires', '_timeout')

    def __init__(self, *args, **kwargs):
        self._timeout = 600

        self._data = dict(*args, **kwargs)
        self._expires = {k: time.time() + self._timeout for k in self._data}

    def __getitem__(self, item):
        self._expire_check()
        return self._data[item]

    def __setitem__(self, item, value):
        self._expire_check()
        self._data[item] = value
        self._expires[item] = time.time() + self._timeout

    def __delitem__(self, item):
        self._expire_check()
        del self._data[item]
        del self._expires[item]

    def __len__(self):
        self._expire_check()
        return len(self._data)

    def __contains__(self, item):
        self._expire_check()
        return item in self._data

    def __iter__(self):
        call_time = time.time()
        for key in list(self._data.keys()):
            if call_time < self._expires.get(key):
                yield key
            else:
                del self._data[key]
                del self._expires[key]

    def keys(self):
        call_time = time.time()
        for key, exp_time in list(self._expires.items()):
            if call_time < exp_time:
                yield key
            else:
                del self._data[key]
                del self._expires[key]

    def items(self):
        call_time = time.time()
        for key, data in list(self._data.items()):
            if call_time < self._expires[key]:
                yield key, data
            else:
                del self._data[key]
                del self._expires[key]

    def values(self):
        call_time = time.time()
        for key, data in list(self._data.items()):
            if call_time < self._expires[key]:
                yield data
            else:
                del self._data[key]
                del self._expires[key]