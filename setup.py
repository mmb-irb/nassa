from setuptools import setup, find_packages
setup(
    name="Nucleic_Acid_Sequence_Structural_Analysis",
    version="0.0.2",
    python_requires=">=3.6.0",
    packages=find_packages(),
    install_requires=[],
    entry_points={
        "console_scripts": [
            "nassa = nassa.scripts.cli:entry_point",
        ]
    },)
