# need to know for exam

### module 1

![alt text](exam/imgs/module1.png)

### module 2

![alt text](exam/imgs/module2.png)

### module 3

![alt text](exam/imgs/module3.png)

### module 4

![alt text](exam/imgs/module4.png)

### global alignemnt

that global thing
[text](vscode-local:/c%3A/Users/Henrik/Downloads/L4_SequenceAlignments_V24.pdf)

### semiglobal thing
terminal gaps shjould not be penalized
align short seq at its position within a large sequence

yea you better do this thing
file:///C:/Users/Henrik/Downloads/slides_BINF100_05_06.pdf

### local aligbnmenbt

with local cannot go negative, so keep 0 as max

that backtrack from cels with max score nbot rugth bottom

ends at cells with 0

### transformations

levenshtein distance is min num of operatoins to make one the other

you can extend gap penalites sometimes, cuz like if u have one gap you may have more probably

![alt text](exam/imgs/module7.png)

### Blast is kinda important

e-value is chance of a lsignment of uinrelated sequences

bnecause we use k-tupkes (n-grams) and smalleer k-tuples makes it so we find less matching sequences (what if query is small?)
higher k will fiond very acuratte 

query coverage means what % iof query matched some part of the match

max ident: % identical neucloetides in match region

file:///C:/Users/Henrik/Downloads/slides_BINF100_07%20(1).pdf
do this excrsie on 16.

![alt text](exam/imgs/module8.png)

mopdule 8!

evaulation of alignmwent scorwes

e val is number of hits with scre higher than s that one can expcect to see by chance when searching database of given size
if database size too big, then e score gaint

s big then gumbel

![alt text](exam/imgs/confusion_matrix.png)