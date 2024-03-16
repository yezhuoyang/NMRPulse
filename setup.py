import setuptools

with open("README.md", "r") as fh:
    description = fh.read()

setuptools.setup(
    name="NMRPulse",
    version="0.0.1",
    author="John Ye",
    author_email="yezhuoyang98@g.ucla.edu",
    packages=["Simulator"],
    description="This is the package for pulse level NMR quantum algorithm simulation.",
    long_description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/yezhuoyang/NMRPulse",
    license='MIT',
    python_requires='>=3.8',
    install_requires=[]
)


