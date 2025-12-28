# Python Tests

All Python tests in this repository use the standard library `unittest` framework.

Run the full suite from the repo root:

```sh
python -m unittest discover -s python/tests -p "test_*.py"
```

Run a single test module:

```sh
python -m unittest python.tests.test_compare_blsolv_strict
```
