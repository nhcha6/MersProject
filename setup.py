from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()


setup(
    name = 'MersProject',
    version = '0.1.0',
    description = 'Software to splice proteins into peptides',
    long_description = readme,
    author = 'Arpit Bajaj, Nic Chapman',
    author_email = 'arpitbajaj98@gmail.com',
    url = 'https://github.com/arpitbajaj98/UROP/',
    license = license,
    packages = find_packages(exclude=('tests', 'docs'))
)