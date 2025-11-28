# from setuptools import setup, find_packages

# setup(
#     name="aaby",
#     version="0.1.0",
#     description="Automated AMBER System Builder Workflow",
#     author="Your Name",
#     packages=find_packages(),
#     include_package_data=True,
#     package_data={"aaby": ["data/tleap.in"]},
#     install_requires=[
#         "pyyaml",
#         # "importlib_resources",  # For Python <3.9
#     ],
#     entry_points={
#         "console_scripts": [
#             "aaby = aaby.main:main",
#         ],
#     },
# )
from setuptools import setup, find_packages

setup(
    name="aaby",
    version="0.1.0",
    description="Automated AMBER System Builder Workflow",
    author="Amanda D. Stange",
    packages=find_packages(),
    include_package_data=True,  # needed so [tool.setuptools.package-data] is honoured
    install_requires=[
        "pyyaml",
        "openbabel",
        "parmed",
        "propka",
        "numpy<1.24",
    ],
    entry_points={
        "console_scripts": [
            "aaby = aaby.main:main",
        ],
    },
)