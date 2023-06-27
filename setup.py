from distutils.core import setup

setup(name='photonpilot',
      version='0.1',
      description='Tool for calculating cavity field strengths in simple realistic cavity setups',
      author='Mark Kamper Svendsen',
      author_email='mark-kamper.svendsen@mpsd.mpg.de',
      url='https://github.com/flickgroup/photonpilot',
      license="GNU",
      install_requires=['numpy'],
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GNU License",
        "Operating System :: OS Independent",
    ],
)