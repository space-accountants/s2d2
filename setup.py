#!/usr/bin/env python
import os
from setuptools import setup

extras_require={
 "test": ["coverage", "mypy", "pycodestyle", "pytest", "pytest-cov",
 "pytest-mock", "pytest-asyncio"]}

# see setup.cfg
setup()
