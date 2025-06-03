ALL: DesignDocument_ClassStructure.html

%.html: %.qmd
	quarto render $<
