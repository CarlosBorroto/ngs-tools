from distutils.core import setup

setup(
    name='454-tools',
    version='0.1.2',
    author='Carlos Borroto',
    author_email='carlos.borroto@gmail.com',
    url="https://github.com/cjav/454-tools",
    packages=['fourfivefour', 'fourfivefour.test'],
    scripts=["scripts/454-tools"],
    description='Compilation of tools for processing 454 sequencing data..',
    long_description=open('README.md').read(),
    install_requires=[
        "docopt",
        "biopython",
        "python-levenshtein"
    ],
)
