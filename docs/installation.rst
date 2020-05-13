Installation
=================

The ferro package is developed on top of Python3, one of the most
`popular open-source languages for scientific computing
<https://spectrum.ieee.org/at-work/innovation/the-2018-top-programming-languages>`_.
This is in part due to its extensibility, with a vast ecosystem of packages available
that add functionality onto the core language (the "standard library").



One of the easiest ways to get started with Python3 is to use the
`distribution released by Anaconda Inc. <https://www.anaconda.com/distribution/>`_,
which includes the standard library as well as several packages
that are popular in scientific computing.

In order to get started, download the appropriate Anaconda distribution
for your OS from the link above.
Once it is installed, open the Anaconda Prompt that was installed with admin
rights and type

..

    conda update conda

followed by enter. Next execute the following command:

..

    pip install --upgrade pip

After that is complete, run

..

    conda create -n ferro

and press enter when it prompts you for confirmation.
This creates a virtual environment for ferro and its dependencies, allowing you
to install them without changing your base installation at all.
This helps prevent the headache that can arise when a new package and an
installed package having conflicting requirements by minimizing the number
of installed packages to only those needed for a particular task.
A separate virtual environment can then be created for each distinct task the
user is working on. To activate the environment after creation, you may enter
the command

..

    conda activate ferro

at the prompt. If this is successful, you should see (ferro) appear to the left of
your directory prompt. You may now install the ferro package by typing

..

    pip install ferro

After this has finished, your environment should be set up. You may now create a script
in your editor of choice. Visual Studio Code and Pycharm are two popular options,
but the Spyder IDE that is installed with Anaconda is perfectly fine to get started with.
At the prompt, type

..

    spyder

to launch the Spyder IDE. Ff you are on Windows and this does not work,
try searching for Spyder in the start menu and launching it from there.

At the beginning of your script, type

..

    | from ferro import data as hd
    | from ferro import models as lf

This tells the python interpreter to import the data and modeling submodules
from ferro and make them accessible with the two-letter abbreviation given.
A simple script to plot a hysteresis measurement looks something like this:

..

   | from ferro import data as hd
   | from ferro import models as lf

   | sampledata = hd.HysteresisData()
   | sampledata.tsv_read('PathToFileOnDiskHere')
   | sampledata.hyst_plot()


You are now set up to create and run scripts for your own data analysis.
Feel free to look at the examples in the test and bin directories of the
GitHub repository for some inspiration and let me know if you have any questions.

Good Luck!