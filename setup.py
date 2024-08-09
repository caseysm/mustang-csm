from setuptools import setup, find_packages
import os
from setuptools.command.install import install
import stat

class CustomInstall(install):
    def run(self):
        install.run(self)
        # Set executable permissions for the MUSTANG binary
        mustang_path = os.path.join(self.install_lib, 'mustang_align', 'bin', 'mustang-3.2.4')
        if os.path.exists(mustang_path):
            st = os.stat(mustang_path)
            os.chmod(mustang_path, st.st_mode | stat.S_IEXEC)

# Read the contents of your README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="mustang_align",
    version="0.1.2",  # Incremented version number
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for performing MUSTANG analysis on PDB files and generating pairwise alignments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/mustang_align",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.7",
    install_requires=[
        "tqdm",
    ],
    entry_points={
        "console_scripts": [
            "mustang_align=mustang_align.mustang_align:main",
        ],
    },
    package_data={
        "mustang_align": ["bin/mustang-3.2.4"],
    },
    cmdclass={
        'install': CustomInstall,
    },
)