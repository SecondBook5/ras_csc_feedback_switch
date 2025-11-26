library( DESeq )
library( genefilter )
library( EBImage )
library( statmod )

setwd("~/Documents/Leo_scRNAseq_SCC_data/SALMON_data_with_Hras_rtTA_mCherry/")

## Read file and split genes and spike-ins

dataMouse <- read.csv( "all_tumor_K14morethan7_Counts_geneCollapsed_filtered.csv", row.names=1, check.names = FALSE )
setdataMouse <- round(dataMouse, digits = 0)
dataMouse[ 1:10, 1:7 ]
dim(dataMouse)

geneTypes <- factor( c( EN="ENSMUSG", ER="ERCC" )[substr( rownames(dataMouse), 1, 2 ) ] )
countsMmus <- dataMouse[ which( geneTypes=="ENSMUSG" ),  ]
countsERCC <- dataMouse[ which( geneTypes=="ERCC" ),  ]

## Normalization (DEseq, sizeFactor)

sfMmus <- estimateSizeFactorsForMatrix( countsMmus )
sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
rbind( sfMmus, sfERCC )

nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsMmus <- t( t(countsMmus) / sfMmus )

## calculate the sample moments

meansERCC <- rowMeans( nCountsERCC )
varsERCC <- rowVars( nCountsERCC )
cv2ERCC <- varsERCC / meansERCC^2

meansMmus <- rowMeans( nCountsMmus )
varsMmus <- rowVars( nCountsMmus )
cv2Mmus <- varsMmus / meansMmus^2

## Fit technical noise

minMeanForFit <- unname( quantile( meansERCC[ which( cv2ERCC > 0.2 ) ], .80 ) )
useForFit <- meansERCC >= minMeanForFit
minMeanForFit
table( useForFit )

fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFit] ),cv2ERCC[useForFit] )

## Test for high variance

minBiolDisp <- 0.25^2

xi <- mean( 1 / sfERCC )
m <- ncol(countsMmus)
psia1theta <- mean( 1 / sfERCC ) + ( coefficients(fit)["a1tilde"] - xi ) * mean( sfERCC / sfMmus )
cv2th <- coefficients(fit)["a0"] + minBiolDisp + coefficients(fit)["a0"] * minBiolDisp
testDenom <- ( meansMmus * psia1theta + meansMmus^2 * cv2th ) / ( 1 + cv2th/m )
p <- 1 - pchisq( varsMmus * (m-1) / testDenom, m-1 )
padj <- p.adjust( p, "BH" )
sig <- padj < .1

sig[is.na(sig)] <- FALSE
table( sig )

## Plot of results

plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-2, 2e4 ), ylim = c( .05, 1100 ),
      xlab = "Average normalized read count", ylab = "Squared coefficient of variation (CV^2)" )
axis( 1, 10^(-2:4), c( "0.01","0.1", "1", "10", "100", "1000",expression(10^4) ))
axis( 2, 10^(-1:3), c( "0.1", "1", "10" ,"100", "1000"), las=2 )
abline( h=10^(-2:2), v=10^(-2:5), col="#D0D0D0", lwd=2 )

# Plot the genes, use a different color if they are highly variable
points( meansMmus, cv2Mmus, pch=20, cex=0.5, col = ifelse( padj < .1, "#C0007090", "#70500040" ) )
# Add the technical noise fit, as before
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, coefficients(fit)["a1tilde"] / xg + coefficients(fit)["a0"], col="#0060B8A0", lwd=3 )
# Add a curve showing the expectation for the chosen biological CV^2 threshold
lines( xg, psia1theta/xg + coefficients(fit)["a0"] + minBiolDisp,lty="dashed", col="#C0007090", lwd=3 )

# Add the normalized ERCC points
points( meansERCC, cv2ERCC, pch=20, cex=1, col="#0060B8A0" )

## table of highly variable genes

log2RelExprMmus <- log2( nCountsMmus / meansMmus )
highVarTable <- data.frame(
  row.names = NULL,
  geneID = rownames(countsMmus)[ sig ],
  meanNormCount = meansMmus[ sig ],
  meanCV2 = cv2Mmus[ sig ],
  strongest = factor( colnames( log2RelExprMmus )[apply( log2RelExprMmus[ sig, ], 1, which.max ) ] ),
  log2RelExprMmus[ sig, ],
  check.names=FALSE )

write.csv(sig,"all_tumor_K14morethan7_Counts_geneCollapsed_filtered_variable.csv")
