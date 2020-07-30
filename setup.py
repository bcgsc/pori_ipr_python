from setuptools import find_packages, setup

# Dependencies required to use your package
INSTALL_REQS = [
    'graphkb>=1.3.6, <2',
    'argparse-env==0.1.0',
    'biopython==1.76',
    'progressbar2==3.51.0',
]

# Dependencies required for development
DEV_REQS = ['flake8', 'black', 'flake8-annotations', 'isort', 'mypy']

# Dependencies required only for running tests
TEST_REQS = ['pytest', 'pytest-cov']

# Dependencies required for deploying to an index server
DEPLOYMENT_REQS = ['twine', 'wheel', 'm2r']

DOC_REQS = ['mkdocs', 'mkdocs-material', 'markdown-refdocs']

long_description = ''
long_description_content_type = 'text/markdown'

try:
    import m2r
    import re

    long_description = m2r.parse_from_file('README.md')
    long_description = re.sub(
        r'.. code-block::.*', '.. code::', long_description
    )  # pyshop has issues with fenced code blocks
    long_description_content_type = 'text/rst'
except ImportError:
    with open('README.md', 'r') as fh:
        long_description = fh.read()


setup(
    name='genomic_report',
    version='2.0.0',
    packages=find_packages(),
    install_requires=INSTALL_REQS,
    extras_require={
        'dev': TEST_REQS + DEPLOYMENT_REQS + DEV_REQS + DOC_REQS,
        'deploy': DEPLOYMENT_REQS,
        'test': TEST_REQS,
        'doc': DOC_REQS,
    },
    package_data={'genomic_report': ['py.typed']},
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    python_requires='>=3.6',
    author_email='creisle',
    author='creisle@bcgsc.ca',
    maintainer='creisle',
    maintainer_email='creisle@bcgsc.ca',
    dependency_links=[],
    test_suite='tests',
    tests_require=TEST_REQS,
    entry_points={'console_scripts': ['genomic_report = genomic_report.main:command_interface']},
)
