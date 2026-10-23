# Contributing

Issues and pull requests are welcome, especially:

- results from running the pipeline on another area (open an issue with the `bbox` you used and what you had to change);
- fixes for Windows or Linux (the pipeline has only been run on macOS, Apple silicon);
- clearer wording in the README.

## Before you open a pull request

```bash
mamba env create -f environment.yml    # or environment.lock.yml for the exact versions used here
mamba activate andrews
python -m pytest tests -q              # offline, takes a few seconds
```

Please keep these conventions:

- Every parameter lives in `config.yml`; scripts read it through `common.load_config()`.
- Shared logic goes in `scripts/common.py` and gets a test in `tests/`.
- Do not commit anything from `data/`, and no personal paths, emails or tokens: the repository is public.
