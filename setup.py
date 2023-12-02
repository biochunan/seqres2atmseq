# setup.py
import os 
from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.install import install

class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        # Calling install.run(self) first will run the standard installation steps
        install.run(self)
        
        # Grant executable permissions to the clustalo executable
        Path(__file__).parent.joinpath('seqres2atmseq', 'assets', 'clustalo').chmod(0o755)
        Path(__file__).parent.joinpath('seqres2atmseq', 'assets', 'clustaloMac').chmod(0o755)
        Path(__file__).parent.joinpath('seqres2atmseq', 'assets', 'clustalo-1.2.4-Ubuntu-x86_64').chmod(0o755)
        Path(__file__).parent.joinpath('seqres2atmseq', 'assets', 'clustal-omega-1.2.3-macosx').chmod(0o755)

setup(
    name='seqres2atmseq',
    version='0.1',
    packages=find_packages(),
    # include non-python files 
    package_data={
        'seqres2atmseq': [
            'assets/clustalo', 
            'assets/clustalo-1.2.4-Ubuntu-x86_64',
            'assets/clustaloMac',
            'assets/clustal-omega-1.2.3-macosx',
            
        ]
    },  
    entry_points={
        'console_scripts': [
            'seqres2atmseq=seqres2atmseq.app:main',  # invoke the main function in seqres2atmseq/app.py
        ],
    },
    # Add other necessary information like author, URL, dependencies, etc.
    cmdclass={
        'install': PostInstallCommand,
    },
    author='ChuNan Liu',
    author_email='chunan.liu@outlook.com',
    description='Align PDB ATOM sequence (ATMSEQ) to SEQRES, and output a mask string mapping residues from SEQRES to ATMSEQ.',
    # add dependency packages
    install_requires=[
        'loguru>=0.7.2',
        'biopython>=1.01',
    ],
)
