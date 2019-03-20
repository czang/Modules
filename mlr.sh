#!/bin/sh

DNASEFILE=DNase_matrix_$1
EXPRFILE=expression_$1
OUTFILE=coefficients_$1

/usr/local/bin/R --no-save -q << EOF

dnase <- read.table("$DNASEFILE")
x <- data.frame(t(dnase[,6:45]))
expre <- read.table("$EXPRFILE")
y <- t(expre[,2:41])
m <- data.frame(x,y)
mlr <- lm(y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20, data=m)
write.table(coefficients(mlr),"$OUTFILE")

EOF
