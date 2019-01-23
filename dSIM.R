####### Try removing or imputing the  missing data from both data matrices before analysing them with dSIM.
####### dep.data : Outcome data matrix with column names as variable names (for example gene names) and row names as patient or sample IDs. Can be copy number, gene expreseeion or any other genomic dataset. The samples (patients) should be on the rows and the variables (for example: genes) on the columns.
####### indep.data : Independent data matrix with column names as variable names (for example gene names) and row names as patient or sample IDs. Can be copy number, gene expreseeion or any other genomic dataset. The samples (patients) should be on the rows and the variables (for example: genes) on the columns.
####### chrarm : Chromosome arm for each variable (for example: genes) in the dep.data matrix. For example: 1p,1q,2a,2p etcetra for each variable. Should be a character vector.
####### chrpos : Chromosome position for each variable (for example: genes) in the dep.data matrix. should be a numeric vector.
####### groups : a numeric vector defining the grouping factor for the samples (for example: a vector with 1s for samples in group 1 and -1s for samples in group 2).
####### minl2 : minimum lambda2 for the ridge penalty (keep it small but not 0, for example start with 0.00000001).
####### maxl2 : maximum lambda2 for the ridge penalty (keep it large, for example: 100000. minl2 and maxl2 defines the boudaries for ridge penalty).
####### nperm : number of permutations to be performed. This number depends on how many samples are there. If not sure how to calculate this then start with a 1000 permutations.
####### alevel : alpha level for Meinshausen permutation correction
####### plot : if T plots the dSIM results in a pdf file with the name "dSIM_results.pdf". default is F. 
####### input.region : this variable can take three values: 
#######         1. "all arms" for plotting the dSIM p-values per chromosome arm for all arms. Default
#######         2. "all arms auto" for plotting the dSIM p-values per chromosome arm for all arms, excluding X and Y chromosomes.  
#######         3. particular chromosome arm for plotting the dSIM p-values only for the chromosome arm. for example, if input.regions = c("12p","12q"), then it will plot the p-values only for chromosome arms 12p and 12q.


######### example: run = dSIM (dep.data = y,indep.data = x,chrarm = chrarm,chrpos = chrpos,groups = groups,minl2 = 0.0000001,maxl2 = 100000,nperm = 1000,alevel=0.05,plot=T,input.region = "all arms")
######### y is the outcome data matrix and x is the dependent data matrix
######### run will be a list with these values:

######### 1. permuted p-value matrix where the rows are the variables and columns are p-values for each permutation (run$pval.mat), 
######### 2. permutation sample indices (run$perm.mat) 
######### 3. Meinshausen selection (run$Mcorr.pval) with either 1 or 0 for each variable. If the value is 0 then the variable is selected by Meinshausen multiple testing correction and hence have a significant associational difference between two groups. 

library(penalized)
library(globaltest)
library(SIM)
library(cherry)

dSIM <- function(dep.data,indep.data,chrarm,chrpos,groups,minl2,maxl2,nperm,alevel=0.05,plot=F,input.region = "all arms"){

        ### dep.data and indep.data have variables on the columns and samples on the rows
        
        ### Initialize
        p.val<- matrix(0,ncol(dep.data),nperm+1)

        ### mean 0 for the Grouping variable
        Ginfo = groups
        Ginfo[groups==as.numeric(names(table(groups))[1])] = (sum(groups==as.numeric(names(table(groups))[2])))/(length(Ginfo))
        Ginfo[groups==as.numeric(names(table(groups))[2])] = -1*(sum(groups==as.numeric(names(table(groups))[1])))/(length(Ginfo))

        ### generate permutation matrix
        
            permT = matrix(0,nrow(dep.data),nperm)
            for(p in 1:ncol(permT)){
              permT[,p] = sample(c(-1,1),nrow(permT),prob = c(0.5,0.5),replace=T)
            }

        for (i in 1:ncol(dep.data)){
            
            response = as.vector(as.matrix(dep.data[,i]))

            prof = optL2(response,penalized=indep.data,minlambda2=minl2,maxlambda2=maxl2)
            pred = (cvl(response,penalized=indep.data,lambda2=prof$lambda))$predictions[,1]
            rm(prof)
            gc()
            resid = response - pred

            C = matrix(0,nrow(indep.data),nrow(indep.data))
            diag(C)= Ginfo
            gtX = C%*%indep.data
            gtX = t(gtX)
            colnames(gtX) = rownames(indep.data)
            rownames(gtX) = colnames(indep.data)
            globj = gt(response=resid,alternative=gtX,null=Ginfo)
            p.val[i,1] = p.value(globj)
            
            ###### Permutations start from here
            
            for (j in 1:nperm){
                  perm = permT[,j]
                  C_perm = perm
        
                  C_perm[C_perm==1] = (sum(perm == -1))/(length(perm))
                  C_perm[C_perm==-1] = -1*(sum(perm == 1))/(length(perm))

                  C = matrix(0,nrow(indep.data),nrow(indep.data))

                  diag(C) = C_perm
    
                  gtX = C%*%indep.data
                  gtX = t(gtX)
                  colnames(gtX) = rownames(indep.data)
                  rownames(gtX) = colnames(indep.data)
                  globj = gt(response=resid,alternative=gtX,null=C_perm)
                  p.val[i,j+1] = p.value(globj)

             }

        }

colnames(p.val) = c("p.obs",paste("p.perm",1:nperm,sep=""))
rownames(p.val) = colnames(dep.data)

#### Meinshausen permutation correction for FDP

p = p.val[,1]
PM = p.val[,2:ncol(p.val)]
pl = unique(curveMeinshausen(p,PM,plot=F,alpha=alevel))
adj_p = rep(1,length(p))
   
  if(length(pl)>1){
     op = order(p)
     adj_p[op[1:pl[length(pl)]]] = 0
  }
      
names(adj_p) = rownames(p.val)
        
if(plot==T){
	pdf("dSIM_results.pdf")
    par(las=2)
    plot(adj_p, main=" ", ylim=c(0,1), pch=20, col="blue", 
    xlab="variables", ylab=paste("Meinshausen selection"),axes="F")
    axis(1,at = 1:length(adj_p),labels = names(adj_p),cex.axis=0.6)	
    axis(2,at = 0:1,labels = 0:1,cex.axis=0.6)
    par(las = 0)
    chrarm = as.character(chrarm)
    chrpos = as.numeric(chrpos)
    dSIMp = adj_p
    ###############################################################################
# function: getGenomicRegion()
# description: converts userdefined region to start and end position
###############################################################################
getGenomicRegion <- function(input.region, rescale=1e9) 
{	
		if(input.region == "whole genome"){
	    region1 <- getGenomicRegion("1p")
	    regionY <- getGenomicRegion("Yq")
	    region <-  data.frame(absolute.start=region1$absolute.start, absolute.end=regionY$absolute.end)
			return(region)
    }
    else if(input.region == "whole genome auto"){
	    region1 <- getGenomicRegion("1p")
	    regionY <- getGenomicRegion("22q")
	    region <-  data.frame(absolute.start=region1$absolute.start, absolute.end=regionY$absolute.end)
			return(region)
    } 
    
    input.region <- as.character(input.region)
    
    extractChr <- function(x) {chr <- unlist(strsplit(x, "[^0-9|xy|XY]{1,2}"))[1]; ifelse(nchar(chr) > 0, chr, NA)}
    extractArm <- function(x) {loc <- regexpr("[p|q]", x); ifelse(loc > 0, substr(x, loc, loc), NA)}
    extractBand <- function(x) {loc <- regexpr("[p|q]", x); ifelse(loc > 0 & loc != nchar(x), substr(x, loc+1, nchar(x)), NA)}
    extractBases <- function(x) {loc <- regexpr(" [0-9]+-[0-9]+", x); unlist(ifelse(loc > 0, strsplit(substr(x, loc+1, nchar(x)),"-"), NA))}
    
    chr <- extractChr(input.region)
    arm <- extractArm(input.region) 
    band <- extractBand(input.region)
    bases <- extractBases(input.region)
    region <- NA
    if(!is.na(chr) & !is.na(arm) & !is.na(band))
        region <- convertGenomicRegion(chr, arm, band) 
    else if(!is.na(chr) & !is.na(arm) & is.na(band))
        region <- convertGenomicRegion(chr, arm)
    else if(!is.na(chr) && !is.na(bases))
        region <- Bases(chr, bases[1], bases[2])
    else if(!is.na(chr) & is.na(arm) & is.na(band))
        region <- convertGenomicRegion(chr)	
    else
        stop("Unknown input region: ", input.region)		
    
    #add absolute genomic scale
    region$absolute.start <- as.numeric(region$start) + as.numeric(region$chr)*rescale 
    region$absolute.end <- as.numeric(region$end) + as.numeric(region$chr)*rescale
    region
}



###############################################################################
# function: predefiendRegions()
# description: transforms predefined input region to chromosomal regions
#              some chromosomal region are less interesting and are excluded
###############################################################################
predefinedRegions <- function(x)
{ 
    x <- as.character(x)
    #convert UCSC format to default
    x <- gsub("^chr", "", x)
    x <- gsub(":", " ", x)	
    switch(x, #`whole genome`=c("1p", "Yq"),  #`whole genome auto`=c("1p", "22q"),
            `all chrs`=c(1:22, "X", "Y"),
            `all chrs auto`=c(1:22),
            `all arms`=c(paste(c(1:22, "X", "Y"), "p", sep=""), paste(c(1:22, "X", "Y"), "q", sep="")) ,
            `all arms auto`=c(paste(c(1:22), "p", sep=""), paste(c(1:22), "q", sep="")),
            x)
}

###############################################################################
# function: convertGenomicRegion()
# description: converts the a genomic region to start and end base pair position
#              using a chromosome table.
###############################################################################
convertGenomicRegion <- function(...) 
{
    data("chrom.table", package="SIM")
    n <- nargs()
    args <- list(...)
    arguments <- list(chr=unique(chrom.table$chr), arm=unique(chrom.table$arm), band=unique(chrom.table$band))
    Columns <- c("start", "end")
    Rows <- TRUE
    for(i in 1:n)
    {
        args[i] <- match.arg(as.character(toupper(args[i])), toupper(arguments[[i]]))		
        Columns <- c(Columns, names(arguments[i]))
        Rows <- Rows & toupper(chrom.table[,names(arguments[i])]) == toupper(args[i])      
    }
    
    select <- chrom.table[Rows, Columns]	
    
    if(nrow(select) < 1)
        stop("Invalid combination of either ", paste(names(arguments), sep=", "))
    
    select[1, "end"] <- select[nrow(select),"end"]
    select[1,]
}




################################################


    plot.dSIM <- function(dSIMp,chrarm,chrpos,input.regions = "all arms"){
    
          
          data("chrom.table", package = "SIM")

          
          #### check the lengths of all inputs
      
          if(!((length(dSIMp)==length(chrarm)) && (length(dSIMp)==length(chrpos)) && (length(chrarm)==length(chrpos)))){
              stop("lengths of dSIMp and annotation data does not match! please provide annotation data per gene")
          }

          input.regions <- unlist(sapply(input.regions, predefinedRegions, 
          USE.NAMES = FALSE))


          nchrs <- length(levels(chrom.table$chr))
          chrsLength <- sapply(levels(chrom.table$chr), function(x) convertGenomicRegion(x)$end)
          xlines <- matrix(c(rep(0, nchrs), chrsLength), nrow = 2, 
          byrow = TRUE)
          ylines <- matrix(c(nchrs:1, nchrs:1), nrow = 2, byrow = TRUE)
          xpoints <- sapply(levels(chrom.table$chr), function(x) convertGenomicRegion(x, 
          "p")$end)
          ypoints <- nchrs:1
          hsegments <- 0.4
          col.chromosomes <- rep("black", nchrs)
          col.centromers <- "purple"
          col.legend <- c("#0088ff", "#9f9f9f")
          par(mar = c(0.1, 4.1, 4.1, 0.1))
          plot.new()
          plot.window(c(0, max(chrsLength)), c(0, nchrs))
          matlines(xlines, ylines, col = col.chromosomes, lty = 1)
   
          for (input.region in input.regions) {
             region <- getGenomicRegion(input.region)
             chrom <- as.numeric(region$chr)
             raw.pvals <- dSIMp[chrarm==input.region]
             abs.start <- chrpos[chrarm==input.region]


             adjpval <- raw.pvals
             xsegments <- matrix(c(abs.start, abs.start), nrow = 2, 
             byrow = TRUE)
             ysegments <- matrix(rep(c(nchrs + 1 - chrom - hsegments, 
             nchrs + 1 - chrom + hsegments), each = length(abs.start)), 
             nrow = 2, byrow = TRUE)
             indices <- which(adjpval == 1)
             if (length(indices) > 0){ 
                 matlines(xsegments[, indices], ysegments[, indices], 
                 col = "grey", lty = 1, pch = 1)}

             indices <- which(adjpval == 0)
             if (length(indices) > 0){ 
                 matlines(xsegments[, indices], ysegments[, indices], 
                 col = "red", lty = 1, pch = 1)}
                 
             lines(c(region$start, region$end), ylines[, chrom], col = "orange", 
             lty = 1)
           }
             points(xpoints, ypoints, col = col.centromers, bg = col.centromers, 
             pch = 19)
             axis(2, at = 1:nchrs, labels = c("Y", "X", (nchrs - 2):1), 
             las = 2)
             title(main = "dSIM selection", ylab = "chromosomes")
             legendText1 <- "selected"
             legendText2 <- "not selected" 
             legend("right", do.call("expression", list(legendText1, legendText2)), fill = c("red","dark grey"), bty = "n")
           }

    plot.dSIM(dSIMp,chrarm,chrpos, input.regions = input.region)
    dev.off()
  }

return(list(pval.mat = p.val,perm.mat = permT,Mcorr.pval = adj_p))

}

