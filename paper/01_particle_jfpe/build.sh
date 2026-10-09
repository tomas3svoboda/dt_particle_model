#!/bin/sh
# Same sequence as build.bat, for Unix-like systems.
set -e
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode supplementary_submission.tex
bibtex supplementary_submission
pdflatex -interaction=nonstopmode supplementary_submission.tex
pdflatex -interaction=nonstopmode supplementary_submission.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode supplementary_submission.tex
