# Requirements

- Linux or Mac OS X
- Python 2.7+ (have not tested on Python 3+)
- GCC 4.6 or higher, with C++ support (G++)


# Recommended installation (self-contained, no admin required)
This method will install SWGA into a Python virtual environment (virtualenv) on your machine. This is the best choice if you are using a computer in which you do not have administrative rights (e.g. a server or cluster). You will still need **GCC 4.6** or higher installed. 

0. Choose a folder to work in. I recommend just creating a new folder; for these instructions I'll pretend it's called `swga_workspace` but you can use whatever folder you want. From a terminal window:
```sh
$ mkdir ~/swga_workspace
```

1. Install `virtualenv` locally:
    ```sh
    $ curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-12.0.tar.gz
    $ tar xvfz virtualenv-12.0.tar.gz
    $ cd virtualenv-12.0
    $ python virtualenv.py ~/swga_workspace
    $ cd ~/swga_workspace
    $ source bin/activate
    ```
    We're now working in a virtual environment, so all the Python packages we need will be installed here and not system-wide (and so we don't need admin privileges). To exit this virtualenv, type `deactivate`. To return to it (e.g. in a new terminal session), simply run `cd ~/swga_workspace; source bin/activate`

2. Download SWGA
    ```sh
    $ curl -LOk https://github.com/eclarke/swga/archive/master.tar.gz
    $ tar xvfz master.tar.gz
    $ cd swga-master
    ```

3. Build and install SWGA

    This part should be easy: just type `make`. 
    The compiler should detect if you're using a Mac or Linux computer automatically. If you encounter problems, make sure you have a relatively up-to-date version of GCC installed, or have installed the Xcode tools on Mac OS X.

4. (optional) Run the tests
    ```sh
    $ py.test 
    ```

5. Return to the SWGA working directory `~/swga_workspace` and run SWGA init:
    ```
    cd ~/swga_workspace
    swga init
    ```
    And just follow the prompts!
