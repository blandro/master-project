from setuptools import setup

setup(
    name="blast_transporters",
    version="1.0.0",
    py_modules=["blast_transporters"],
    install_requires=["pandas"],
    entry_points={
        "console_scripts": [
            "blast_transporters = blast_transporters:main",
        ],
    },
)