from setuptools import setup, find_packages

setup(
    name='mmonitor',
    version='0.3.0',
    description='MMonitor - Metagenomic Monitoring Tool',
    author='MMonitor Team',
    author_email='',
    url='https://github.com/lucast122/MMonitor',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    package_data={
        'mmonitor': [
            # Resources
            'resources/*.json',
            'resources/images/*.png',
            'resources/icons/*.icns',
            'resources/icons/*.png',
            'resources/icons/**/*.png',
            'resources/themes/*.json',
            'resources/r/*',
            # Workflows
            'workflows/Snakefile',
            'workflows/Dockerfile',
            'workflows/rules/*.smk',
            'workflows/envs/*.yaml',
            'workflows/envs/*.yml',
            'workflows/config/*.yaml',
            'workflows/config/*.yml',
            'workflows/scripts/*',
            'workflows/scripts/**/*',
        ],
    },
    python_requires='>=3.10',
    install_requires=[
        'click>=8.1',
        'pandas>=2.0',
        'numpy>=1.24',
        'requests>=2.31',
        'keyring>=24.0',
        'osfclient>=0.0.5',
        'watchdog>=3.0',
        'pyyaml>=6.0',
        'tqdm>=4.65',
        'biopython>=1.81',
        'snakemake>=8.0',
    ],
    extras_require={
        'gui': [
            'customtkinter>=5.2',
            'CTkMessagebox>=2.5',
            'pillow>=10.0',
        ],
        'dev': [
            'pytest>=7.0',
            'pytest-cov>=4.0',
            'black>=23.0',
            'isort>=5.12',
            'flake8>=6.0',
        ],
    },
    entry_points={
        'console_scripts': [
            # New CLI (v0.3)
            'mmonitor=mmonitor.cli.main:main',
            # Legacy CLI (deprecated, shows migration message)
            'mmonitor-cmd=mmonitor.userside.run_cmd:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
