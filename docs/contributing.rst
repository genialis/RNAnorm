============
Contributing
============

Preparing environment
=====================

If you don't have Git installed on your system, follow `these
instructions <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.

`Fork <https://help.github.com/articles/fork-a-repo>`_ the main repository.
Clone your fork (replace ``<username>`` with your GitHub account name) and
change directory::

   git clone git@github.com:<username>/RNAnorm.git
   cd RNAnorm

Install the project for development::

   pip install -e '.[docs,package,test]'

.. note::
   We recommend using `venv <http://docs.python.org/3/library/venv.html>`_
   to create an isolated Python environment.


Ensuring reproducibility
========================

Project package dependencies should be defined in ``pyproject.toml``.
Recommended is specifying a minimal version with some ability to update, i.e.
using a tilde requirement (e.g., ``pandas~=1.4.2``).


Running tests
=============

Run unit tests and linters::

   pytest --verbose --cov=rnanorm
   black --check src tests
   flake8 src tests
   pydocstyle src
   isort --check-only --diff src tests

To run the tests with Tox, use::

   tox -r


Building documentation
======================

.. code::

   sphinx-build docs build/sphinx


Preparing release
=================


Describe the new features in ``changelog.rst``. Replace the Unreleased heading
with the new version, followed by the release date (e.g.,
``13.2.0 - 2018-10-23``).

Add the new dependencies to ``pyproject.toml`` and update the package version.
Use `semantic versioning <https://semver.org>`_.

Commit changes to Git::

   git commit -a -m "Prepare release <new-version>"

Push changes to your fork and open a pull request::

   git push --set-upstream <RNAnorm-fork-name> release-<new-version>

Wait for tests to pass and pull request to be approved. Merge code to main::

   git checkout main
   git merge --ff-only release-<new-version>
   git push <RNAnorm-upstream-name> main <new-version>

Tag the new release from the latest commit::

   git checkout main
   git tag -sm "Version <new-version>" <new-version>

Push the tag to the upstream repository::

   git push <RNAnorm-upstream-name> main <new-version>
