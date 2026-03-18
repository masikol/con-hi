# AGENTS.md - Agent Coding Guidelines for con-hi

This file provides guidelines for agentic coding agents working on the con-hi codebase.

## Project Overview

- **Project name**: con-hi (consensus-highlighter)
- **Language**: Python 3.6+
- **Entry point**: `./con-hi.py`
- **Purpose**: Annotates low-coverage and high-coverage regions in FASTA files using read mapping in BAM format
- **Main dependencies**: Biopython, samtools 1.11+

## Build, Lint, and Test Commands

### Running Tests

Run all tests with pytest:
```bash
python3 -m pytest tests/
```

Run a single test file:
```bash
python3 -m pytest tests/test_arguments.py
```

Run a single test by name:
```bash
python3 -m pytest tests/test_arguments.py::test_function_name
```

Run tests with verbose output:
```bash
python3 -m pytest tests/ -v
```

### Dependencies

Install runtime dependencies:
```bash
pip3 install biopython
```

Install test dependencies:
```bash
pip3 install pytest
```

### Running the Application

Basic usage:
```bash
./con-hi.py -f <TARGET_FASTA> -b <MAPPING_BAM>
```

With custom coverage thresholds:
```bash
./con-hi.py -f my_sequence.fasta -b my_mapping.sorted.bam -c 25,55 -X 1.5,2.0
```

## Code Style Guidelines

### General Principles

- Python 3.6+ compatible code only
- Use type hints from the `typing` module
- Use `logging` module for all logging (not print statements)
- Use f-strings for string formatting where appropriate, or `.format()` method
- Keep lines under 100 characters when reasonable

### Naming Conventions

- **Functions/variables**: `snake_case` (e.g., `parse_arguments`, `target_fasta_fpath`)
- **Classes**: `PascalCase` (e.g., `HighlighterArgs`, `CoverageThreshold`)
- **Constants**: `UPPER_SNAKE_CASE` (e.g., `_ARG_SEP`)
- **Private functions**: prefix with underscore (e.g., `_parse_options`)

### Import Organization

Organize imports in the following order with a blank line between groups:

1. Standard library imports
2. Third-party imports
3. Local project imports

Example:
```python
import os
import sys
import logging
from typing import List, Sequence

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

import src.output as out
import src.obtain_coverage as oc
from src.main import main
```

### Type Hints

- Always use type hints for function parameters and return types
- Use `typing` module for complex types: `List`, `Sequence`, `Set`, `Optional`, etc.
- Inline variable type annotations are encouraged:
  ```python
  args: HighlighterArgs = parse_arguments()
  opts: List[List[str]] = getopt.gnu_getopt(...)
  ```

### Code Block End Markers

Always use end markers for blocks (if, for, def, class):
```python
if condition:
    do_something()
# end if

for item in items:
    process(item)
# end for

def function():
    return value
# end def

class MyClass:
    pass
# end class
```

### Error Handling

- Use `logging` for errors, warnings, and info messages
- Exit with appropriate error codes using `platf_depend_exit(code)`
- Catch specific exceptions rather than using bare except
- Example:
  ```python
  try:
      os.makedirs(outdpath)
  except OSError as err:
      logging.error(f'Cannot create output directory `{outdpath}`.')
      logging.error(str(err))
      platf_depend_exit(1)
  # end try
  ```

### Docstrings

- Use simple docstrings for documentation
- Include parameter descriptions and return values for complex functions

### File Organization

- Main source code in `src/` directory
- Tests in `tests/` directory
- Test fixtures in `tests/fixtures.py`
- Coverage threshold classes in `src/coverage_thresholds/`

### Testing Guidelines

- Use pytest fixtures (prefix with `@pytest.fixture`)
- Test file naming: `test_<module_name>.py`
- Group related tests in classes if appropriate
- Use descriptive test function names: `test_<what_is_being_tested>`

### Git Conventions

- Do not commit directly unless explicitly requested
- Run tests before committing
- Do not commit generated files, caches, or temporary files
- Follow existing commit message style (check `git log` for examples)

## Architecture Notes

- `con-hi.py`: Entry point, handles version check and dependency validation
- `src/main.py`: Main orchestration logic
- `src/arguments.py`: Command-line argument parsing using `getopt`
- `src/coverage_thresholds/`: Coverage threshold classes (inherit from `CoverageThreshold`)
- `src/obtain_coverage.py`: Uses samtools to obtain coverage data
- `src/highlight_features.py`: Identifies regions matching coverage criteria
- `src/output.py`: Writes GenBank output files
