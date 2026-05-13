from setuptools import setup, find_packages

setup(
    name="modes",
    version="0.1.0",
    description="MoDES: Multi-Omics Discordance/Event State inference",
    long_description=open("README.md").read() if __import__("os").path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    author="MoDES Team",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21",
        "scipy>=1.7",
        "pandas>=1.3",
        "statsmodels>=0.13",
        "anndata>=0.8",
        "matplotlib>=3.5",
        "seaborn>=0.11",
        "networkx>=2.6",
    ],
    python_requires=">=3.7",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
