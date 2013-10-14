from setuptools import setup

version = __import__('ngstools.version').get_version()

setup(
    name='ngs-tools',
    version=version,
    author='Carlos Borroto',
    author_email='carlos.borroto@gmail.com',
    url="https://github.com/cjav/ngs-tools",
    packages=['ngstools', 'ngstools.test'],
    scripts=["scripts/ngs-tools"],
    description='Compilation of tools for processing NGS sequencing data..',
    long_description=open('README.md').read(),
    install_requires=[
        "docopt",
        "biopython",
        "python-levenshtein"
    ],
)
