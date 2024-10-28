How to use the lyx thesis class:

	Note: in this document, replace <path> with the path to the directory that lyx was installed to. 

1) Copy thesis.layout to the <path>lyx\share\lyx\layouts directory.

2) Copy thesis.cls to the <path>lyx\share\lyx\tex directory.

3) Add the following line to the textclass.lst file stored in the <path>lyx\share\lyx directory: 
	"thesis" "<path>/thesis" "thesis"

4) Start lyx - the thesis class should be available by selecting it under "document class" by clicking Layout->Document. The recommended font is pslatex.

5) If you'd like to see an example of the thesis class, open the thesis_template.lyx file. To edit the title page/abstract/etc, modify the preamble by clicking Layout->Document->Preamble and inserting the desired text.

Eric Sampson,
University of Saskatchewan
2004/10/19 