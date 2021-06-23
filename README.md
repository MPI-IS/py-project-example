# py-project-example

This repository serves as a minimal example on how to structure you python project.

## Requirements

This project requires `python >= 3.5`.

## Installation

You can perform a user-specific installation by running

    $ python -m pip install .

from the root of the project. We strongly advise an installation in a virtual environment. You can create and activate one by executing the following two commands from the root of the project

    $ python -m venv venv
    $ . venv/bin/activate

and then performing the installation as usual by running

    (venv) $ python -m pip install .

If you plan to extend the code, then you should perform an editable installation with

    (venv) $ python -m pip install -e .

## Testing

You can run the unit-tests by executing

    $ python -m unittest nlse.tests

from the root of the project. The tests include some standard problems of pulse propagation in nonlinear media. During the tests an interactive plotter demonstrating the integration results will be shown. Unfortunately, at the moment it is not possible to disable it, so running tests in a headless setup is not straightforward.

## Useful References

- [.gitignore](https://git-scm.com/docs/gitignore) -- how to ignore files in git. Github has some pretty good defaults when creating the repository.
- [Markdown](https://www.markdownguide.org/basic-syntax/) -- how to use Markdown format to write `README.md`.
- [Packaging tutorial](https://packaging.python.org/tutorials/packaging-projects/) -- an official tutorial from Python Packaging Authority.
