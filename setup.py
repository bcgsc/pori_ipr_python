from setuptools import find_packages, setup

# Dependencies required to use your package
INSTALL_REQS = [
    'graphkb>=1.5.1, <2',
    'biopython==1.76',
    'progressbar2>=3.51.0, <4',
    'pandas>=1.1.0, <2',
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
    with open('README.md', 'r') as fh:
        long_description = fh.read()
except Exception:
    pass

setup(
    name='ipr',
    version='2.0.2',
    packages=find_packages(),
    install_requires=INSTALL_REQS,
    extras_require={
        'dev': TEST_REQS + DEPLOYMENT_REQS + DEV_REQS + DOC_REQS,
        'deploy': DEPLOYMENT_REQS,
        'test': TEST_REQS,
        'doc': DOC_REQS,
    },
    package_data={'ipr': ['py.typed']},
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    python_requires='>=3.6',
    author_email='graphkb@bcgsc.ca',
    author='graphkb',
    maintainer='graphkb',
    maintainer_email='graphkb@bcgsc.ca',
    dependency_links=[],
    test_suite='tests',
    tests_require=TEST_REQS,
    entry_points={'console_scripts': ['ipr = ipr.main:command_interface']},
)
