Thank you for your interest in contributing to `linkinglines`. Please feel free to open up issues with bugs or requested features. Any contributions you make will benefit everybody else and are greatly appreciated.

We recommend using a virutal environment to manage packages `virtualenv`, see [here](https://virtualenv.pypa.io/en/latest/). 

```
# Install virtualenv if you haven't already
pip install virtualenv

# Navigate to the project directory
cd path/to/LinkingLines

# Create a virtual environment
virtualenv venv

# Activate the virtual environment
# On Windows
venv\Scripts\activate
# On Unix or MacOS
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

If you would like to contribute code please do so in a seperate branch and open up an issue describing your contribution.

```
git clone git@github.com:aikubo/LinkingLines.git
git checkout my-development-branch
```

Before submitting your pull request please verify the following:

1. Code is documented in [NumPy Docstring Style](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html)
2. Code is tested and passes test 
    - To run the tests please go to "/tests" and run `pytest`
    - Add your test code to any file with the name `test`
    - More here on [pytest and testing practices](https://docs.pytest.org/en/8.0.x/)
3. Open an issue and pull request 

After your pull request the code will be reviewed by maintainers. After the sucessful merge the documentation will be regenerated.
