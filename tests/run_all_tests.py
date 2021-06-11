"""
This file runs all the tests defined in the /tests/ subsirectoris

"""
import time
import sys
import pytest

if __name__ == "__main__":
    print("Running all Tweezepy tests.")
    start = time.time()
    pytest.main()
    end = time.time()
    print("-------------------------")
    print("All tests done in %2.3f s" % (end-start))    