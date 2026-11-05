# BOLOS documentation

- [Collision reconstruction sign audit](COLLISION_SIGN_VALIDATION.md): independent
  integral and matrix checks, the adopted sign correction, and its numerical
  impact compared with the legacy reconstruction.
- [BOLOS vs BOLSIG+ validation notebook](../validate_bolsig.ipynb): executed
  comparisons of EEDFs, selected reaction rates, and electron transport,
  including mesh sensitivity checks and exported tables and figures.
- [Source code responsibilities and execution flow](CODE_GUIDE.md): an English
  guide to all six modules in `bolos/`, their data flow, numerical processing,
  and complete calculation examples.
- [Tutorial source](source/tutorial.rst): the existing Sphinx tutorial.
- [API documentation source](source/api.rst): the existing Sphinx API reference.
- [Sphinx index](source/index.rst): the entry point for the existing
  reStructuredText documentation build.

`CODE_GUIDE.md` is intended to be read directly in a Markdown viewer. It is not
currently included in the Sphinx build, which uses reStructuredText sources.
