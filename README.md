# Scripts and Notebooks for Electrical methods

I've created this repository to store, manage, and share some Python scripts and notebooks (coming soon) for testing, modelling & understanding VES and ERT modelling and inversion. All of the geophysical processine elements use the [pyGIMLi](https://www.pygimli.org/) package, and many of the scripts and notebooks draw extensively from examples published by the pyGIMLi team.

With the exception of the 'For_Colab' folder, all files are intended to be used from a locally installed python environment with pyGIMLi and its dependencies installed (as well as Jupyter notebook where relevant). Instructions for installation can be found at "[pyGIMLi Installation](https://www.pygimli.org/installation.html)" while Jupyter notebook can be installed with ```conda install conda-forge::notebook```. As far as I can remember, Python scripts in "[scripts](https://github.com/NoteboomM/VES-ERT_Acacia/tree/master/scripts)" and Jupyter notebooks in "[Notebooks](https://github.com/NoteboomM/VES-ERT_Acacia/tree/master/notebooks)" have been tested and work with a local installation. The For_Colab folder contains some adapted versions of notebooks for use with Google Colab, however that implementation has issues. As of 17/10/2024 only the 1D_VES notebook runs on Colab without critical issues.

Otherwise, "[exampledata](https://github.com/NoteboomM/VES-ERT_Acacia/tree/master/exampledata)" stores some data files used as inputs and references, and "[working](https://github.com/NoteboomM/VES-ERT_Acacia/tree/master/working)" contains various scratch files and works-in-progress.

Matt Noteboom

Gouda, 17/10/2024
