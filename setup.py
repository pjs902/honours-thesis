from setuptools import find_packages, setup # pragma: no cover

setup( # pragma: no cover
    name="binaryshift",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
)
