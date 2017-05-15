from collections import MutableMapping
import time


class TimeoutDict(MutableMapping):
    '''Dictionary that expires entries after _timeout seconds unless the dictionary is touched (set/get)'''

    def _expire_check(self):
        '''Check whether all keys are expired, if so, delete the key
            Note: this evaluates all the keys, and deletes expired ones, not just the item queried
            probably a bit redundant, and less efficient than it could be if we passed item and just expired that'''
        call_time = time.time()
        for key, exp_time in list(self._expires.items()):
            if call_time >= exp_time:
                del self._data[key]
                del self._expires[key]

    __slots__ = ('_data', '_expires', '_timeout')

    def __init__(self, *args, **kwargs):
        '''Initialize timeout, dictionary and expiries'''
        self._timeout = 600

        self._data = dict(*args, **kwargs)
        self._expires = {k: time.time() + self._timeout for k in self._data}

    def __getitem__(self, item):
        '''get with expire check'''
        self._expire_check()
        return self._data[item]

    def __setitem__(self, item, value):
        '''set item and expiry'''
        self._expire_check()
        self._data[item] = value
        self._expires[item] = time.time() + self._timeout

    def __delitem__(self, item):
        '''delete item, this also expires other keys if needed'''
        self._expire_check()
        del self._data[item]
        del self._expires[item]

    def __len__(self):
        '''len and expire check'''
        self._expire_check()
        return len(self._data)

    def __contains__(self, item):
        '''contain and expire check'''
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


import unittest


class TestTimeoutDict(unittest.TestCase):
    def test_timeout(self):
        '''Test to ensure that timeout works as expected'''
        d = TimeoutDict()
        d._timeout = 6  # make test run faster!
        d[0] = 0  # assign d[0], set its time
        self.assertEqual(len(d), 1)  # there is one item in dict
        time.sleep(1)  # after 1 second
        d[1] = 1  # assign d[1], set its time
        self.assertEqual(d[1], 1)  # d[0] is valid
        self.assertEqual(d[0], 0)  # d[1] is valid
        self.assertEqual(len(d), 2)  # there are two item in dict
        time.sleep(5)  # after 5 additional seconds
        self.assertEqual(d[1], 1)  # d[1] is still valid
        # there is only one item in dict, even though we haven't asked for d[0]
        # yet
        self.assertEqual(len(d), 1)
        self.assertRaises(KeyError, d.__getitem__, 0)  # d[0] has expired


suite = unittest.TestLoader().loadTestsFromTestCase(TestTimeoutDict)
unittest.TextTestRunner(verbosity=2).run(suite)
