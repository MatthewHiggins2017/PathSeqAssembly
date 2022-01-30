import setuptools
import os

setuptools.setup(

    name="PathSeqAssembly",
    version="1.0.0",
    packages=["DeNovoCore"],
    license="MIT",
    long_description="Hybrid assembly pipeline.",
    scripts= ["scripts/%s" % x for x in os.listdir("scripts")],

)
