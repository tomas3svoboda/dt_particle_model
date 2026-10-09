@echo off
rem Build the manuscript (main.pdf) and the Supplementary Material
rem (supplementary_submission.pdf). The two documents cross-reference each
rem other through the xr package, so each is compiled after the other's .aux
rem exists. Requires pdflatex, bibtex and the apacite, natbib, siunitx and xr
rem packages (any current TeX Live or MiKTeX).
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
