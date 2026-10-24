"""Analyses — concrete rule sets for specific languages.

Each module here instantiates the framework for one analysis and
builds a :class:`phonology.profile.AnalysisProfile` for it. A new
analysis supplies its language-specific steps through the profile's
hooks, so it needs no edits outside its own file beyond any new
label it adds to the shared ``Normalisation`` / ``RefinedCategory``
vocabularies.
"""
