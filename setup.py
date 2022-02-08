from setuptools import find_packages, setup  # pragma: no cover

setup(  # pragma: no cover
    name="binaryshift",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    author="Peter Smith",
    email="smith.peter.902@gmail.com",
    uri="https://github.com/pjs902/honors-thesis",
    license="MIT",
    description="Shift mass bins to emulate binary stars",
    copyright="Copyright 2021 Peter Smith",
    python_requires=">=3.8",
    version="0.2.0",
    include_package_data=True,
    package_data={"binaryshift": ["resources/*.dat"]},
)
