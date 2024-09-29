from setuptools import setup, find_packages

setup(
    # ... (other setup parameters)
    entry_points={
        'console_scripts': [
            'mmonitor-cmd=mmonitor.userside.run_cmd:main',
        ],
    },
)