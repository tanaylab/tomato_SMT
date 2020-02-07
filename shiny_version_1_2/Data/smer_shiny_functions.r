

########################################################################
### this function gets two vectors of at least one gene in  each,
### and plot their normalized expression per 100K UMIs in each meristem, colored by its stage
###
###
smer_shiny_plot_xy <- function(
						 norm_mat = NULL,
						 rel_metadata = NULL,
						 bg_plot=FALSE,
						 genes_a = NULL,
						 genes_b = NULL,
						 genes_a_label = "label A",
						 genes_b_label = "label B",						 
						 log2_exp = FALSE,
						 color_by_dev=FALSE,
						 color_na=FALSE,
						 reg_log=1,
						 test_meristems=FALSE,
						 show_legend=FALSE,
						 show_grid=FALSE,
						 legend_pos="topright",
						 my_xlim=NULL,
						 my_ylim=NULL,
						 my_pch=19,
						 my_cex=0.8
						 )
{
	if (bg_plot){
		par(mar=c(5,5,1,1))
		plot(c(1:5),c(1:5),
			col="#00000080",
			pch=my_pch,
			cex.axis=1.5,cex.lab=1.5,
			xlab="gene list A",
			ylab="gene list B",
			cex=my_cex)
	} else {

if (color_na){rel_metadata$color=NA}		

	if (!log2_exp){
		if(length(genes_a)==1 & length(genes_b)==1) {
		par(mar=c(5,5,1,1))
		plot(
			norm_mat[genes_a,],
			norm_mat[genes_b,],
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex.axis=1.5,cex.lab=1.5,
			xlab=genes_a_label,
			ylab=genes_b_label,
			xlim=list( c(min(c(norm_mat[genes_a,])), max(c(norm_mat[genes_a,]))),
				my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
			ylim=list( c(min(c(norm_mat[genes_b,])), max(c(norm_mat[genes_b,]))),
				my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
			cex=my_cex)
		} else if (length(genes_a) > 1 & length(genes_b) == 1){
		par(mar=c(5,5,1,1))
		plot(
			colSums(norm_mat[genes_a,]),
			norm_mat[genes_b,],
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex.axis=1.5,cex.lab=1.5,
			xlab=genes_a_label,
			ylab=genes_b_label,
				xlim=list( c(min(c(colSums(norm_mat[genes_a,]))), max(c(colSums(norm_mat[genes_a,])))),
			my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
				ylim=list( c(min(c(norm_mat[genes_b,])), max(c(norm_mat[genes_b,]))),
			my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
			cex=my_cex)
		} else if (length(genes_a) == 1 & length(genes_b) > 1){
		par(mar=c(5,5,1,1))
			plot(
			norm_mat[genes_a,],
			colSums(norm_mat[genes_b,]),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex.axis=1.5,cex.lab=1.5,
			xlab=genes_a_label,
			ylab=genes_b_label,
			xlim=list( c(min(c(norm_mat[genes_a,])), max(c(norm_mat[genes_a,]))),
				my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
			ylim=list( c(min(c(colSums(norm_mat[genes_b,]))), max(c(colSums(norm_mat[genes_b,])))),
				my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
			cex=my_cex)
		} else {
			par(mar=c(5,5,1,1))
			plot(
			colSums(norm_mat[genes_a,]),
			colSums(norm_mat[genes_b,]),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex.axis=1.5,cex.lab=1.5,
			xlab=genes_a_label,
			ylab=genes_b_label,
			xlim=list( c(min(c(colSums(norm_mat[genes_a,]))), max(c(colSums(norm_mat[genes_a,])))),
				my_xlim)[[ifelse(is.null(my_xlim),1,2)]], 
			ylim = list( c(min(c(colSums(norm_mat[genes_a,]))), max(c(colSums(norm_mat[genes_a,])))),
				my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
			cex=my_cex)	
		}
		
		if (show_grid){
			grid()
			}
			
		if (show_legend & color_by_dev ) {
			legend(x=legend_pos,
					x.intersp=0.2,
					y.intersp=0.8,
					bty="n",
					fill=list(c("green","darkgreen","yellow","orange","red"),c("green","darkgreen","yellow","orange","grey","red"))[[ifelse(test_meristems,2,1)]],
					legend=list(c("VM","TM0","TM1","TM2","FM"),c("VM","TM0","TM1","TM2","TM","FM"))[[ifelse(test_meristems,2,1)]],
					cex=1)}
	
	} else {					 
	
		if(length(genes_a)==1 & length(genes_b)==1) {
		par(mar=c(5,5,1,1))
		plot(
			log2(norm_mat[genes_a,]+reg_log),
			log2(norm_mat[genes_b,]+reg_log),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex.axis=1.5,cex.lab=1.5,
			xlab=genes_a_label,
			ylab=genes_b_label,
			xlim=list( c(min(c(log2(norm_mat[genes_a,]+reg_log))), max(c(log2(norm_mat[genes_a,]+reg_log)))),
				my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
			ylim=list( c(min(c(log2(norm_mat[genes_b,]+reg_log))), max(c(log2(norm_mat[genes_b,]+reg_log)))),
				my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
			cex=my_cex)
		} else if (length(genes_a) > 1 & length(genes_b) == 1){
		par(mar=c(5,5,1,1))
		plot(
			log2(colSums(norm_mat[genes_a,])+reg_log),
			log2(norm_mat[genes_b,]+reg_log),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex.axis=1.5,cex.lab=1.5,
			xlab=genes_a_label,
			ylab=genes_b_label,
			xlim=list( c(min(c(log2(colSums(norm_mat[genes_a,])+reg_log))), max(c(log2(colSums(norm_mat[genes_a,])+reg_log)))),
				my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
			ylim=list( c(min(c(log2(norm_mat[genes_b,]+reg_log))), max(c(log2(norm_mat[genes_b,]+reg_log)))),
				my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
			cex=my_cex)	
		} else if (length(genes_a) == 1 & length(genes_b) > 1){
			par(mar=c(5,5,1,1))
			plot(
			log2(norm_mat[genes_a,]+reg_log),
			log2(colSums(norm_mat[genes_b,])+reg_log),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex.axis=1.5,cex.lab=1.5,
			xlab=genes_a_label,
			ylab=genes_b_label,
			xlim=list( c(min(c(log2(norm_mat[genes_a,]+reg_log))), max(c(log2(norm_mat[genes_a,]+reg_log)))),
				my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
			ylim=list( c(min(c(log2(colSums(norm_mat[genes_b,])+reg_log))), max(c(log2(colSums(norm_mat[genes_b,])+reg_log)))),
				my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
			cex=my_cex)
		} else {
			par(mar=c(5,5,1,1))
			plot(
			log2(colSums(norm_mat[genes_a,])+reg_log),
			log2(colSums(norm_mat[genes_b,])+reg_log),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex.axis=1.5,cex.lab=1.5,
			xlab=genes_a_label,
			ylab=genes_b_label,
			cex=my_cex)
		}
		
		if (show_grid){
		grid()
		}
	
		if (show_legend & color_by_dev) {
			legend(x=legend_pos,
					x.intersp=0.2,
					y.intersp=0.8,
					bty="n",
					fill=list(c("green","darkgreen","yellow","orange","red"),c("green","darkgreen","yellow","orange","grey","red"))[[ifelse(test_meristems,2,1)]],
					legend=list(c("VM","TM0","TM1","TM2","FM"),c("VM","TM0","TM1","TM2","TM","FM"))[[ifelse(test_meristems,2,1)]],
					cex=1)
		}
	
	}
 }
}


# This functions adds points (of different genotypes), to a current WT scatter
####################################################################################################################################################################
smer_shiny_add_points_xy <- function(
						 norm_mat = NULL,
						 rel_metadata = NULL,
						 genes_a = NULL,
						 genes_b = NULL,
						 genes_a_label = "label A",
						 genes_b_label = "label B",						 
						 log2_exp = FALSE,
						 color_by_dev=FALSE,
						 reg_log=1,
						 my_pch=19,
						 my_cex=0.8){
						 

	if (!log2_exp){
	
		if(length(genes_a)==1 & length(genes_b)==1) {
		points(
			norm_mat[genes_a,],
			norm_mat[genes_b,],
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex=my_cex)
		} else if (length(genes_a) > 1 & length(genes_b) == 1){
		points(
			colSums(norm_mat[genes_a,]),
			norm_mat[genes_b,],
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex=my_cex)
		} else if (length(genes_a) == 1 & length(genes_b) > 1){
		points(
			norm_mat[genes_a,],
			colSums(norm_mat[genes_b,]),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex=my_cex)
		} else {
		points(
			colSums(norm_mat[genes_a,]),
			colSums(norm_mat[genes_b,]),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex=my_cex)	
		}
		grid()

	} else {					 
	
		if(length(genes_a)==1 & length(genes_b)==1) {
		points(
			log2(norm_mat[genes_a,]+reg_log),
			log2(norm_mat[genes_b,]+reg_log),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex=my_cex)
		} else if (length(genes_a) > 1 & length(genes_b) == 1){
		points(
			log2(colSums(norm_mat[genes_a,])+reg_log),
			log2(norm_mat[genes_b,]+reg_log),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex=my_cex)	
		} else if (length(genes_a) == 1 & length(genes_b) > 1){
			points(
			log2(norm_mat[genes_a,]+reg_log),
			log2(colSums(norm_mat[genes_b,])+reg_log),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex=my_cex)
		} else {
			points(
			log2(colSums(norm_mat[genes_a,])+reg_log),
			log2(colSums(norm_mat[genes_b,])+reg_log),
			col=list(rel_metadata$color,"#00000080")[[ifelse(color_by_dev,1,2)]],
			pch=my_pch,
			cex=my_cex)
		}
		grid()
	
		}
	
	}
####################################################################################################################################################################
####################################################################################################################################################################
####################################################################################################################################################################
 
# This function draws ordered expression of a single gene
####################################################################################################################################################################
smer_shiny_plot_single_gene_by_order <- function (expression_vec,
												  show_trendline = FALSE,
												  k_trendline = 11,
												  gene_name,
												  my_ylim = NULL,
												  meristems_metadata = smer_md,
												  bg_group_colors = NULL,
												  bg_group_fracs = NULL,
												  break_by_groups = NULL)
{

if (!is.null(bg_group_colors) & !is.null(bg_group_fracs)){
bg_color = bg_group_colors
bg_frac = bg_group_fracs
}

if (is.null(break_by_groups)) {
xs = 1:length(expression_vec)
	} else {
	xgroups = length(break_by_groups)+1
	space_group = length(expression_vec) / xgroups
	xpos = c();
		for (i in 1:xgroups){
	
			if (i==1) {group_n = break_by_groups[1]; current_range= c(0,space_group)}
			else if (i==xgroups) {group_n = length(expression_vec) - break_by_groups[xgroups-1];
								current_range = c((i*space_group)-space_group,i*space_group); 
			} else {
			current_range = c((i*space_group)-space_group,i*space_group); 
			group_n = break_by_groups[i] - break_by_groups[i-1]
			}
		xpos = c(xpos,seq(current_range[1],current_range[2],length.out = group_n))
		}
	xs = xpos;
	}

plot(
        xs,  
		expression_vec, 
		pch=19,cex=1.5,cex.main=1.5,cex.lab=1.5,
		xlab="ordered meristems",
		ylab="expression",
		main = gene_name,
		xlim = c(1,length(expression_vec)),
		ylim = list(range(expression_vec),my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
		col = meristems_metadata[names(expression_vec),"color"])
		
if (show_trendline) {
	  cat(paste0(pryr::mem_used(),"..PLOTTING TREND1..\n"))
	lines(xs,
	rollmean(expression_vec, k = k_trendline, na.pad=TRUE),
		   lwd=4, col="black") }

if( (!is.null(bg_group_colors) & !is.null(bg_group_fracs))){
	  cat(paste0(pryr::mem_used(),"..PLOTTING BACKGROUND #1..\n"))
	bg_lims= cumsum(c(par("usr")[1],cumsum(par("usr")[2]-par("usr")[1])*bg_frac))
			for (i in 1:length(bg_color)){
			rect(bg_lims[i], par("usr")[3], bg_lims[i+1], par("usr")[4], col = bg_color[i], border = bg_color[i])
			}
	}
			   

}
####################################################################################################################################################################
####################################################################################################################################################################
####################################################################################################################################################################
 
 
 
 
# This function draws barplots of gene averaged expression per genotype
####################################################################################################################################################################
smer_shiny_gene_barplots <- function(gene = NULL,
									 order_by_genotypes = FALSE,
									 show_se = FALSE,
									 text_cex=1.5,
									 raw_umi_matrix = NULL,
									 metadata_table = NULL,
									 bar_width=1,
									 se_lwd=3,
									 ci_width=0.3)
	#								 maximal_umi_per_meristem = 2e6)
{

	if (is.null (gene)) {stop ("Undefined gene!")}
	if (is.null (raw_umi_matrix)) {stop ("Undefined UMI matrix!")}
	if (is.null (metadata_table)) {stop ("Undefined Metadata table!")}
	
# remove ufs
raw_umi_matrix = raw_umi_matrix[,!metadata_table$genotype=="uf"]
metadata_table = metadata_table[!metadata_table$genotype=="uf",]
	

umis = raw_umi_matrix
stages = metadata_table$group

if (!order_by_genotypes){
	levels_ordered = c("dst sft_veget.",
					   "WT_veget.","dst_veget.","sft_veget.",
					   "WT_trans1","dst_trans1","sft_trans1",
					   "WT_trans2","dst_trans2","sft_trans2",
					   "WT_flower","dst_flower","sft_flower")
	levels_color = c(rep("green2",4),
					 rep("yellow2",3),
					 rep("orange",3),
					 rep("red",3))
	levels_nms = c("dst;sft",
					"wt-veg.","dst-veg.","sft-veg.",
					"wt-tr1", "dst-tr1", "sft-tr1",
					"wt-tr2", "dst-tr2", "sft-tr2",
					"wt-flo.","dst-flo.","sft-flo.") 
	levels_spaces = c(0.2,
					rep(c(0.5,0.2,0.2),4))
					
} else {

	levels_ordered = c(
				  "WT_veget.", "WT_trans1", "WT_trans2", "WT_flower",
				  "dst_veget.","dst_trans1","dst_trans2","dst_flower",
				  "sft_veget.","sft_trans1","sft_trans2","sft_flower",
				  "dst sft_veget.")
	levels_color = c(rep(c("green2","yellow2","orange","red"),4),"green2")
	levels_nms = c("wt-veg.","wt-tr1","wt-tr2","wt-flo.",
					"dst-veg.","dst-tr1","dst-tr2","dst-flo.",
					"sft-veg.","sft-tr1","sft-tr2","sft-flo.",
					"dst;sft")
	levels_spaces=c(0.2,rep(c(0.2,0.2,0.2,0.5),3))
} 

stages = factor(stages, levels = levels_ordered)
#	if (!is.null (maximal_umi_per_meristem)) {
#	umis_comb = cbind(
#				 umis[,colSums(umis)< maximal_umi_per_meristem],
#				 smer_shiny_ds(umis[,colSums(umis)>=maximal_umi_per_meristem],maximal_umi_per_meristem)
#				 )
#	umis = umis_comb[,colnames(umis)]
#	}
stages_total = tapply(colSums(umis),stages,sum)
gene_total = tapply(umis[gene,],stages,sum)
cat(paste0("total gene UMIs: ",sum(gene_total),"\n"))

stages_norm = (gene_total / stages_total) * 1e5
p_stage = gene_total / stages_total
se_p = ((p_stage*(1-p_stage)) / sqrt (stages_total))

stages_norm_sqrtK = (sqrt(gene_total) / stages_total) * 1e5

par(mgp = c(3,0.7,0),mar=c(6,4,4,2), bg = NA)
barplot2(stages_norm,
		plot.ci = show_se,
		ci.u = stages_norm + stages_norm_sqrtK,
		ci.l = stages_norm - stages_norm_sqrtK,
		col = levels_color,
		names = levels_nms,
		space = levels_spaces,
		cex.lab = text_cex,
		cex.axis = text_cex,
		cex.names = text_cex,
		width = bar_width,
		ci.lwd = se_lwd,
		ci.width = ci_width,
		ylab = "",
		main = gene,
		cex.main = text_cex,
		yaxt = "n",
		las = 2)
magaxis(side=2, cex = text_cex, cex.axis = text_cex, cex.lab = text_cex, line=0, las=2)
title(ylab = "UMIs per 100K UMIs",cex.lab = text_cex,line=2)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################






## This function downsample a umi_tab to ds_val
#################################################################################################################
## smer_shiny_ds <- function (umi_tab, ds_val)
## {
## umi_tab_filt = umi_tab[,colSums(umi_tab)>=ds_val]
## umi_tab_ds = apply(umi_tab_filt, 2, .downsamp_one, ds_val)
## rownames(umi_tab_ds)=rownames(umi_tab);
## return(umi_tab_ds)
## }
## #################################################################################################################
## .downsamp_one=function(v,n)
## {
##   hist(sample(rep(1:length(v),times=v),replace=F,size=n),0.5+0:length(v),plot=F)$counts
## }
##############################################################################################
##############################################################################################

