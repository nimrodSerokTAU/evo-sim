from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="indel-simulator",
    version="1.0.0",
    author="Elya Wygoda, Asher Moshe, Nimrod Serok, Edo Dotan, Noa Ecker, Omer Israeli, Itsik Pe'er, Tal Pupko",
    author_email="talp@tauex.tau.ac.il",
    description="Efficient CLI tool for simulating indel events along phylogenetic trees",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nimrodSerokTAU/evo-sim",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Academic Free License (AFL)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21.0",
        "matplotlib>=3.5.0",
        "ete3>=3.1.2",
        "llist>=0.7.1",
        "typing_extensions>=4.14.0",
    ],
    license="Academic Free License v. 3.0",
    entry_points={
        "console_scripts": [
            "indel-simulator=indelsim.indel_simulator:main",
        ],
    },
)

