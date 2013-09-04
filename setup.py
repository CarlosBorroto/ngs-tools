from setuptools import setup

setup(
    name='ngs-tools',
    version='0.1.4',
    author='Carlos Borroto',
    author_email='carlos.borroto@gmail.com',
    url="https://github.com/cjav/ngs-tools",
    packages=['ngs', 'ngs.test'],
    scripts=["scripts/ngs-tools"],
    description='Compilation of tools for processing NGS sequencing data..',
    long_description=open('README.md').read(),
    install_requires=[
        "docopt",
        "biopython",
        "python-levenshtein"
    ],
)
