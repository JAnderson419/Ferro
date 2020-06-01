Tutorials in this directory are versioned as py files and synced to md and jupyter notebook representations 
with jupytext. If viewing on binder, jupytext will sync the python files to notebooks on build. 
If running locally, execute the following from shell in this directory
to generate the notebook files:

```
pip install -r tut-requirements.txt
jupytext --sync ./pyscript/*.py
```