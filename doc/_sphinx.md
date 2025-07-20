# Sphinx usage

## Install

```
$ conda install sphinx
```

According to the document <https://www.sphinx-doc.org/en/master/usage/markdown.html>, we need [MyST-Parser](https://myst-parser.readthedocs.io) for MarkDown.

```
$ conda install -c conda-forge myst-parser
```

For the Read-the-Docs theme,

```
$ conda install -c conda-forge sphinx_rtd_theme
```

Autobuild <https://qiita.com/Nenshu_agetai/items/3dc7b2d7717e7bb07c47>

```
$ conda install -c conda-forge sphinx-autobuild
```

`linkify-it-py` is needed by `linkify` extension in `myst-parser`.

<https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#linkify>


### One command for mamba

```
$ mamba install sphinx myst-parser sphinx_rtd_theme sphinx-autobuild linkify-it-py
```

# Build the document

## HTML

```
$ make html
```

or 

```
$ sphinx-build -b html . ./_build/html
```

To update the HTML files automatically when source files are edited, enable autobuild by

```
sphinx-autobuild -b html . ./_build/html
```

## PDF via latex

```
$ make latexpdf
```

or 

```
$ sphinx-build -b latex . ./_build/latex
```
