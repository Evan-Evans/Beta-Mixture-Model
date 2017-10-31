# change directory 
setwd("/path/to/working/directory")
source("beta_mixture_model.r")

# read the input file
# for you own analysis, you will just need to provide the input file with mandatory columns
# mandatory values: CHR_ID, START, END, CELL_POSI_INFO in first three and the last columns, respectively
x <- read.table("test.candidate.txt")


mat.out <- c()
# read the information of one region each time (each row stores one candidate region) 
for(i in 1:nrow(x)) {
    # read the last column "CELL_POSI_INFO" and store each cell into one element of the vector 
	cells <- strsplit(as.character(x[i,ncol(x)]), ";")[[1]] 

	mat.list <- list(NULL); cell.names <- c()
	# read the information of each cell
	for(j in 1:length(cells)) {
		# read the methylation calls of each CpG site in cell j
		meth.count <- strsplit(strsplit(cells[j], ":")[[1]], "_")
		meth.count <- unlist(meth.count)

		# the cell.names will store the cell names (5 in this case)
		cell.names <- c(cell.names, meth.count[1])
		meth.count <- meth.count[2:length(meth.count)]
		meth.count <- matrix(as.numeric(meth.count), length(meth.count)/3, 3, byrow=T)[,2:3]
		# meth.count includes the methylated calls and total calls for each CpG site in region i of cell j
		# rows represent all CpG sites in region i of cell j, columns represent the methylated calls and total calls

		# the mat.list will store the meth.count for all cells 
		if(is.null(mat.list[[1]])) mat.list <- list(meth.count) else mat.list <- c(mat.list, list(meth.count))
	}
	names(mat.list) <- cell.names
  
# analyze the region i by performing beta_mixture_model function
  output <- unlist(beta_mixture_model(mat.list))
# calculate the confidence interval of methylation variance and combined to output
  ci <- quantile(variance_boot_ci(mat.list), c(0.025, 0.975)) 
  mat.out <- rbind(mat.out, c(output, ci)) 
}

# formatting and outputting
label <- colnames(mat.out)
rownames(mat.out) <- rownames(x)
# calculate adjusted p-value and combined to output, and add the first three columns "CHR_ID", "START", and "END" into output
mat.out <- cbind(x[,1:3], mat.out, p.adjust(as.numeric(paste(mat.out[,14])), method="BH")) 
colnames(mat.out) <- c("chr", "start", "end", label, "LRT.pval.adjusted") 
# for your own analysis, replace the output file name
write.table(mat.out, file="test.candidate.BetaMixtureModelResult.txt", quote=F, sep="\t", row.names=F, col.names=T)

# get csm with the following cutoffs:
# 1,valid cluster1 and cluster2: avg.m1 != -1; avg.m2 != -1
# 2,observed minimum methylation difference of two methylation states: min.delta >0.1
# 3,threshold of methylation difference of two methylation states: theta1–theta2 >=0.3 
# 4,required number of cells: cell.num >=8
# 5,adjusted p-value: LRT.pval.adjusted <0.05
mat.out <- read.table("test.candidate.BetaMixtureModelResult.txt",h=T,stringsAsFactors=FALSE)
x.out = mat.out[mat.out$min.delta > 0.1 & mat.out$avg.m1 != -1 & mat.out$avg.m2 != -1 & mat.out$cell.num >= 8 & mat.out$LRT.pval.adjusted < 0.05 & mat.out$theta1 - mat.out$theta2 >= 0.3,]
# for your own analysis, replace the output file name
write.table(x.out, file="test.csm.txt", quote=F, sep="\t", row.names=F, col.names=T)
