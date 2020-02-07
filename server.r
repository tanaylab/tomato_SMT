cat(paste0(pryr::mem_used(),"..\n"))

require(shiny)
cat(paste0(pryr::mem_used(),"..\n"))

require(shinyWidgets)
cat(paste0(pryr::mem_used(),"..\n"))

require(zoo)
cat(paste0(pryr::mem_used(),"..\n"))

require(pryr)
cat(paste0(pryr::mem_used(),"..\n"))

require(magicaxis)
cat(paste0(pryr::mem_used(),"..\n"))

require(gplots)
cat(paste0(pryr::mem_used(),"..\n"))


# load source functions
source("./Data/smer_shiny_functions.r")

# load expression matrices
cat(paste0(pryr::mem_used(),"..\n"))
cat("loading expression matrices (part1)..\n")
load("./Data/expression_matrices_part1.Rda")
cat(paste0(pryr::mem_used(),"..\n"))
cat("loading expression matrices (part2)..\n")
load("./Data/expression_matrices_part2.Rda")
cat(paste0(pryr::mem_used(),"..\n"))
#
cat("loading semi-downsampled matrix..\n")
load("./Data/semiDownsampled_matrix.Rda")
cat(paste0(pryr::mem_used(),"..\n"))
##objects loaded:
#smer_n = read.csv("Data/..)
#smer_umi = read.csv("Data/..)
#smer_md = read.csv("Data/..)
#smer_squeezed_meristems (maximal 2e6 UMIs per meristem)
###################################


# Loading tables with gene-gene correlations
cat("loading gene correlations (part1)..\n")
load("./Data/corr_objects_part1.Rda")
cat(paste0(pryr::mem_used(),"..\n"))
#
cat("loading gene correlations (part2)..\n")
load("./Data/corr_objects_part2.Rda")
cat(paste0(pryr::mem_used(),"..\n"))
####

flag_plot1=FALSE 
cat("DONE!\n")
##### finsihed loading #####



shinyServer(function(input, output) {

  meristem_genotype_Input1 <- reactive({
    switch(input$genotype_meristems1,
           "WT" = rownames(smer_md)[smer_md$genotype=="WT"],
           "dst" = rownames(smer_md)[smer_md$genotype=="dst"],
           "sft" = rownames(smer_md)[smer_md$genotype=="sft"],
           "dst;sft" = rownames(smer_md)[smer_md$genotype=="dst sft"],
           "uf" = rownames(smer_md)[smer_md$genotype=="uf"])
  })  
  
  meristem_genotype_Input2 <- reactive({
    switch(input$genotype_meristems2,
           "WT" = rownames(smer_md)[smer_md$genotype=="WT"],
           "dst" = rownames(smer_md)[smer_md$genotype=="dst"],
           "sft" = rownames(smer_md)[smer_md$genotype=="sft"],
           "dst;sft" = rownames(smer_md)[smer_md$genotype=="dst sft"],
           "uf" = rownames(smer_md)[smer_md$genotype=="uf"])
  })  
  
  meristem_genotype_Input3 <- reactive({
    switch(input$genotype_meristems3,
           "WT" = rownames(smer_md)[smer_md$genotype=="WT"],
           "dst" = rownames(smer_md)[smer_md$genotype=="dst"],
           "sft" = rownames(smer_md)[smer_md$genotype=="sft"],
           "dst;sft" = rownames(smer_md)[smer_md$genotype=="dst sft"],
           "uf" = rownames(smer_md)[smer_md$genotype=="uf"])
  })  
   # dependent on expression type, select matrix
    datasetInput <- reactive({
			switch(input$exp_type,
			       "KNN-norm." =  as.matrix(smer_knn),
			       "norm. counts" = as.matrix(smer_n),
			       "log-norm. counts" = as.matrix(log2((smer_n*7)+1)),
			       "raw counts" = as.matrix(smer_umi))
				   	})
	
    # dependent on genotype, select genes-correlation matrix

    corr_datasetInput <- reactive({
      switch(input$cor_genotype,
             "WT" = wt_pos_vals,
             "dst" = dst_pos_vals,
             "sft" = sft_pos_vals)
    })    				
    corr_datasetInput_nms <- reactive({
      switch(input$cor_genotype,
             "WT" = wt_pos_genes,
             "dst" = dst_pos_genes,
             "sft" = sft_pos_genes)
    })
    
    anticorr_datasetInput <- reactive({
      switch(input$cor_genotype,
            "WT" = wt_neg_vals,
            "dst" = dst_neg_vals,
            "sft" = sft_neg_vals)
    })    				
    anticorr_datasetInput_nms <- reactive({
      switch(input$cor_genotype,
             "WT"  = wt_neg_genes,
             "dst" = dst_neg_genes,
             "sft" = sft_neg_genes)
    })
    
    
# screen for 3 gene names inputs  
   formulaText <- reactive({
    rownames(datasetInput())[grepl(input$gene,rownames(datasetInput()))]
  })
  
    formulaText2 <- reactive({
    rownames(datasetInput())[grepl(input$gene2,rownames(datasetInput()))]
  })
  
    formulaText3 <- reactive({
    rownames(datasetInput())[grepl(input$gene3,rownames(datasetInput()))]
  })
  
#gene1 
  output$caption_1 <- renderText({
   if (length(formulaText())!=1) {
   "gene #1 is not a valid Solyc gene name"} else {
   formulaText()}
  })
  
    output$flag_1 <- reactive({
    (length(formulaText())==1)
  })
      outputOptions(output, "flag_1", suspendWhenHidden = FALSE) 

  output$ordered_gene_Plot_1 <- renderPlot({
   if (length(formulaText())==1) {
     cat(paste0(pryr::mem_used(),"..PLOTTING GENE1..\n"))
     covered_meristems = colnames(smer_umi)[colSums(smer_umi)>=(input$min_mer*1e3)]
     relevant_meristems_plot1 = intersect(meristem_genotype_Input1(),covered_meristems)
     
     all_genes_expression = c(datasetInput()[formulaText(),relevant_meristems_plot1])
     if (length(formulaText2())==1 & length(formulaText3())!=1)   {
       relevant_meristems_plot2 = intersect(meristem_genotype_Input2(),covered_meristems)
       all_genes_expression = c(all_genes_expression,datasetInput()[formulaText2(),relevant_meristems_plot2])} 
     else if (length(formulaText3())==1 & length(formulaText2())!=1)   {
       relevant_meristems_plot3 = intersect(meristem_genotype_Input3(),covered_meristems)
       all_genes_expression = c(all_genes_expression,datasetInput()[formulaText3(),relevant_meristems_plot3])}         
     else if (length(formulaText3())==1 & length(formulaText2())==1)   {
       relevant_meristems_plot2 = intersect(meristem_genotype_Input2(),covered_meristems)
       relevant_meristems_plot3 = intersect(meristem_genotype_Input3(),covered_meristems)
       all_genes_expression = c(all_genes_expression,datasetInput()[formulaText2(),relevant_meristems_plot2],datasetInput()[formulaText3(),relevant_meristems_plot3]) }
     fixed_expression = range(all_genes_expression)
     if (sum(input$order_plot_feats %in% "plot_bg")==1) {
      ordered_stages = list( c("veget.","trans1","trans2"), c("veget.","uf"))[[ifelse(input$genotype_meristems1=="uf",2,1)]]
      groups_table = table(smer_md[relevant_meristems_plot1,"group"])[table(smer_md[relevant_meristems_plot1,"group"])>0]
      group_sz = cumsum(groups_table[match(ordered_stages,gsub("[a-z]+_","",tolower(names(groups_table))))])
        if (input$genotype_meristems1=="dst;sft"){
          group_sz = rep(length(relevant_meristems_plot1),4)
        }
      uf_bg_colors = c("#00ff0035","#ffdd0035","#ff000035")
      non_uf_bg_colors = c("#00ff0035","#ffff0035","#ffA50035","#ff000035")
      
      smer_shiny_plot_single_gene_by_order (expression_vec = datasetInput()[formulaText(),relevant_meristems_plot1],
                                            show_trendline = (sum(input$order_plot_feats %in% "trend_bool")==1) & length(formulaText())==1,
                                            gene_name = formulaText(),
                                            meristems_metadata = smer_md,
                                            k_trendline = input$smooth_n,
                                            my_ylim = list(NULL,fixed_expression)[[ifelse((sum(input$order_plot_feats %in% "fixed_ylim")==1),2,1)]],
                                            bg_group_colors = list(uf_bg_colors,non_uf_bg_colors)[[ifelse(input$genotype_meristems1=="uf",1,2)]],
                                            bg_group_fracs = list(c(1/4,1/2,1/4),rep(1/4,4))[[ifelse(input$genotype_meristems1=="uf",1,2)]],
                                            break_by_groups = group_sz)
  
    } else {
      smer_shiny_plot_single_gene_by_order (expression_vec = datasetInput()[formulaText(),relevant_meristems_plot1],
                                             show_trendline = (sum(input$order_plot_feats %in% "trend_bool")==1) & length(formulaText())==1,
                                             gene_name = formulaText(),
                                             my_ylim = list(NULL,fixed_expression)[[ifelse((sum(input$order_plot_feats %in% "fixed_ylim")==1),2,1)]],                                            
                                             meristems_metadata = smer_md,
                                             k_trendline = input$smooth_n)
     }
    }
    
	
	flag_plot1=TRUE;
	})#, width=900,height=300)

#gene2
	  output$caption_2 <- renderText({
   if (length(formulaText2())!=1)  {
   "gene #2 is not a valid Solyc gene name"} else {
   formulaText2()}
  })

    output$flag_2 <- reactive({
    (length(formulaText2())==1)
  })
      outputOptions(output, "flag_2", suspendWhenHidden = FALSE)

 #   output$gene_2_removed <- reactive({
 #   (input$gene_2 & input$gene_2_off)
 # })	  
 # observEvent(input$gene_2_,{input$gene_2=FALSE})
 #     outputOptions(output, "gene_2_removed", suspendWhenHidden = FALSE)
  
  output$ordered_gene_Plot_2 <- renderPlot({
   if (length(formulaText2())==1) {
     cat(paste0(pryr::mem_used(),"..PLOTTING GENE2..\n"))
     covered_meristems = colnames(smer_umi)[colSums(smer_umi)>=(input$min_mer*1e3)]
     relevant_meristems_plot2 = intersect(meristem_genotype_Input2(),covered_meristems)

     fixed_expression = c()
     if (length(formulaText())==1 & (sum(input$order_plot_feats %in% "fixed_ylim")==1))   {
       relevant_meristems_plot1 = intersect(meristem_genotype_Input1(),covered_meristems)
       all_genes_expression = c(datasetInput()[formulaText(),relevant_meristems_plot1])
       if (length(formulaText2())==1 & length(formulaText3())!=1)   {
         relevant_meristems_plot2 = intersect(meristem_genotype_Input2(),covered_meristems)
         all_genes_expression = c(all_genes_expression,datasetInput()[formulaText2(),relevant_meristems_plot2])} 
       else if (length(formulaText3())==1 & length(formulaText2())!=1)   {
         relevant_meristems_plot3 = intersect(meristem_genotype_Input3(),covered_meristems)
         all_genes_expression = c(all_genes_expression,datasetInput()[formulaText3(),relevant_meristems_plot3])}         
       else if (length(formulaText3())==1 & length(formulaText2())==1)   {
         relevant_meristems_plot2 = intersect(meristem_genotype_Input2(),covered_meristems)
         relevant_meristems_plot3 = intersect(meristem_genotype_Input3(),covered_meristems)
         all_genes_expression = c(all_genes_expression,datasetInput()[formulaText2(),relevant_meristems_plot2],datasetInput()[formulaText3(),relevant_meristems_plot3]) }
       fixed_expression = range(all_genes_expression)
     }
          
     if (sum(input$order_plot_feats %in% "plot_bg")==1) {
       ordered_stages = list( c("veget.","trans1","trans2"), c("veget.","uf"))[[ifelse(input$genotype_meristems2=="uf",2,1)]]
       groups_table = table(smer_md[relevant_meristems_plot2,"group"])[table(smer_md[relevant_meristems_plot2,"group"])>0]
       group_sz = cumsum(groups_table[match(ordered_stages,gsub("[a-z]+_","",tolower(names(groups_table))))])
       if (input$genotype_meristems2=="dst;sft"){
         group_sz = rep(length(relevant_meristems_plot2),4)
       }
       uf_bg_colors = c("#00ff0035","#ffdd0035","#ff000035")
       non_uf_bg_colors = c("#00ff0035","#ffff0035","#ffA50035","#ff000035")

       smer_shiny_plot_single_gene_by_order (expression_vec = datasetInput()[formulaText2(),relevant_meristems_plot2],
                                             show_trendline = (sum(input$order_plot_feats %in% "trend_bool")==1) & length(formulaText2())==1,
                                             gene_name = formulaText2(),
                                             meristems_metadata = smer_md,
                                             k_trendline = input$smooth_n,
                                             my_ylim = list(NULL,fixed_expression)[[ifelse((sum(input$order_plot_feats %in% "fixed_ylim")==1),2,1)]],
                                             bg_group_colors = list(uf_bg_colors,non_uf_bg_colors)[[ifelse(input$genotype_meristems2=="uf",1,2)]],
                                             bg_group_fracs = list(c(1/4,1/2,1/4),rep(1/4,4))[[ifelse(input$genotype_meristems2=="uf",1,2)]],
                                             break_by_groups = group_sz)
       
     } else {
       smer_shiny_plot_single_gene_by_order (expression_vec = datasetInput()[formulaText2(),relevant_meristems_plot2],
                                             show_trendline = (sum(input$order_plot_feats %in% "trend_bool")==1) & length(formulaText2())==1,
                                             gene_name = formulaText2(),
                                             my_ylim = list(NULL,fixed_expression)[[ifelse((sum(input$order_plot_feats %in% "fixed_ylim")==1),2,1)]],                                            
                                             meristems_metadata = smer_md,
                                             k_trendline = input$smooth_n)
     }
	}
	
  })#, width=900,height=300)
 
#gene3
	  output$caption_3 <- renderText({
   if (length(formulaText3())!=1)  {
   "gene #3 is not a valid Solyc gene name"} else {
   formulaText3()}
  })

  output$flag_3 <- reactive({
    (length(formulaText3())==1)
  })
      outputOptions(output, "flag_3", suspendWhenHidden = FALSE) 

  output$ordered_gene_Plot_3 <- renderPlot({
   if (length(formulaText3())==1) {
     cat(paste0(pryr::mem_used(),"..PLOTTING GENE3..\n"))
     covered_meristems = colnames(smer_umi)[colSums(smer_umi)>=(input$min_mer*1e3)]
     relevant_meristems_plot3 = intersect(meristem_genotype_Input3(),covered_meristems)

     fixed_expression = c()
     if (length(formulaText())==1 & (sum(input$order_plot_feats %in% "fixed_ylim")==1))   {
       relevant_meristems_plot1 = intersect(meristem_genotype_Input1(),covered_meristems)
       all_genes_expression = c(datasetInput()[formulaText(),relevant_meristems_plot1])
       if (length(formulaText2())==1 & length(formulaText3())!=1)   {
         relevant_meristems_plot2 = intersect(meristem_genotype_Input2(),covered_meristems)
         all_genes_expression = c(all_genes_expression,datasetInput()[formulaText2(),relevant_meristems_plot2])} 
       else if (length(formulaText3())==1 & length(formulaText2())!=1)   {
         relevant_meristems_plot3 = intersect(meristem_genotype_Input3(),covered_meristems)
         all_genes_expression = c(all_genes_expression,datasetInput()[formulaText3(),relevant_meristems_plot3])}         
       else if (length(formulaText3())==1 & length(formulaText2())==1)   {
         relevant_meristems_plot2 = intersect(meristem_genotype_Input2(),covered_meristems)
         relevant_meristems_plot3 = intersect(meristem_genotype_Input3(),covered_meristems)
         all_genes_expression = c(all_genes_expression,datasetInput()[formulaText2(),relevant_meristems_plot2],datasetInput()[formulaText3(),relevant_meristems_plot3]) }
       fixed_expression = range(all_genes_expression)
     }
     
   if (sum(input$order_plot_feats %in% "plot_bg")==1) {
       ordered_stages = list( c("veget.","trans1","trans2"), c("veget.","uf"))[[ifelse(input$genotype_meristems3=="uf",2,1)]]
       groups_table = table(smer_md[relevant_meristems_plot3,"group"])[table(smer_md[relevant_meristems_plot3,"group"])>0]
       group_sz = cumsum(groups_table[match(ordered_stages,gsub("[a-z]+_","",tolower(names(groups_table))))])
       if (input$genotype_meristems3=="dst;sft"){
         group_sz = rep(length(relevant_meristems_plot3),4)
       }
       uf_bg_colors = c("#00ff0035","#ffdd0035","#ff000035")
       non_uf_bg_colors = c("#00ff0035","#ffff0035","#ffA50035","#ff000035")
 
       smer_shiny_plot_single_gene_by_order (expression_vec = datasetInput()[formulaText3(),relevant_meristems_plot3],
                                             show_trendline = (sum(input$order_plot_feats %in% "trend_bool")==1) & length(formulaText3())==1,
                                             gene_name = formulaText3(),
                                             meristems_metadata = smer_md,
                                             k_trendline = input$smooth_n,
                                             my_ylim = list(NULL,fixed_expression)[[ifelse((sum(input$order_plot_feats %in% "fixed_ylim")==1),2,1)]],
                                             bg_group_colors = list(uf_bg_colors,non_uf_bg_colors)[[ifelse(input$genotype_meristems3=="uf",1,2)]],
                                             bg_group_fracs = list(c(1/4,1/2,1/4),rep(1/4,4))[[ifelse(input$genotype_meristems3=="uf",1,2)]],
                                             break_by_groups = group_sz)
       
     } else {
       smer_shiny_plot_single_gene_by_order (expression_vec = datasetInput()[formulaText3(),relevant_meristems_plot3],
                                             show_trendline = (sum(input$order_plot_feats %in% "trend_bool")==1) & length(formulaText3())==1,
                                             gene_name = formulaText3(),
                                             my_ylim = list(NULL,fixed_expression)[[ifelse((sum(input$order_plot_feats %in% "fixed_ylim")==1),2,1)]],
                                             meristems_metadata = smer_md,
                                             k_trendline = input$smooth_n)
     }
    }


  })#, width=900,height=300)

#screen for name of gene for correlation computation:
  formulaText_corr <- reactive({
    rownames(datasetInput())[grepl(input$gene_corr,rownames(datasetInput()))]
  })
  
  formulaText_corr_coverage <- reactive({
    sum(grepl(input$gene_corr,rownames(corr_datasetInput_nms())))
  })
  
# create caption for gene (correlation)
  output$caption_corr <- renderText({
    if (length(formulaText_corr())!=1) {
      "not a valid Solyc gene name"} 
        else {
          formulaText_corr();
        }
  })

# create flag for gene_corr    
  output$flag_corr <- reactive({
    (length(formulaText_corr())==1) #& sum(formulaText_corr_coverage())>=1)
  })
  outputOptions(output, "flag_corr", suspendWhenHidden = FALSE) 
  
### subset up-correlation  
  output$corr_table <- renderTable(expr = {
    if (length(formulaText_corr())==1 & formulaText_corr_coverage()==0) {
      "No enough stats to support corr.";
    } else if (length(formulaText_corr())==1) {
    temp_cor = rev(corr_datasetInput()[formulaText_corr(),])
    temp_cor_nms = rev(corr_datasetInput_nms()[formulaText_corr(),] )
	gtot_genes = list(wt_gtot,dst_gtot,sft_gtot)[[ ifelse(input$cor_genotype=="sft",3,ifelse(input$cor_genotype=="dst",2,1))]]
    data.frame(gene = substring(temp_cor_nms[1:input$up_corr],1,28), 
               corr =  temp_cor[1:input$up_corr], 
               total_UMI = paste0( round(as.numeric(gsub("\\..*","",gtot_genes[temp_cor_nms][1:input$up_corr]))/1000),"k") )   
    }
  }, include.colnames=FALSE)
  
### subset anti-correlation  
  output$corr_table_anti <- renderTable(expr = {
    #output$corr_table <- renderTable({
      if (length(formulaText_corr())==1 & formulaText_corr_coverage()==0) {
        expr = "(gene covered by less than 100 UMIs)"# colnames=FALSE;
      } else if (length(formulaText_corr())==1) {
    temp_cor = anticorr_datasetInput()[formulaText_corr(),]  
    temp_cor_nms = anticorr_datasetInput_nms()[formulaText_corr(),]  
	gtot_genes = list(wt_gtot,dst_gtot,sft_gtot)[[ ifelse(input$cor_genotype=="sft",3,ifelse(input$cor_genotype=="dst",2,1))]]
    data.frame(gene = substring(temp_cor_nms[1:input$down_corr],1,28), 
			   corr =  temp_cor[1:input$down_corr], 
			   total_UMI = paste0(round(as.numeric(gsub("\\..*","",gtot_genes[temp_cor_nms][1:input$down_corr]))/1000,digits=1),"k") ) 
      }                        
  }, include.colnames=FALSE)
    


### gene list A text
  formulaText_listA <- reactive({
    text_lista = unlist(strsplit(input$genes_a,split="\n"))
    if (length(text_lista)==0) { uniq_lista =FALSE }else{
    uniq_lista = sapply(text_lista,function(x){sum(grepl(x,rownames(datasetInput())))==1})}
  if (sum(!uniq_lista)>0 | sum(uniq_lista)==0)
    {"there is a non-valid gene name in list A"} else{
    genes_to_plot_a = unique(rownames(datasetInput())[grepl(paste0(text_lista,collapse="|"),rownames(datasetInput()))]);
    #paste0(unique(rownames(datasetInput())[grepl(paste0(text_lista,collapse="|"),rownames(datasetInput()))]),collapse="\n");
    #paste(paste0(genes_to_plot_a,collapse="\n"),"tot: ",length(genes_to_plot_a))}
    paste0("total genes A: ",length(genes_to_plot_a))}
  })

### gene list B text
  formulaText_listB <- reactive({
    text_listb = unlist(strsplit(input$genes_b,split="\n"))
    if (length(text_listb)==0) { uniq_listb =FALSE }else{
      uniq_listb = sapply(text_listb,function(x){sum(grepl(x,rownames(datasetInput())))==1})}
    if (sum(!uniq_listb)>0 | sum(uniq_listb)==0)
    {"there is a non-valid gene name in list B"} else{
      genes_to_plot_b = unique(rownames(datasetInput())[grepl(paste0(text_listb,collapse="|"),rownames(datasetInput()))]);
      paste0("total genes B: ",length(genes_to_plot_b))}
  })  

### gene list B text for plot
  formulaText_listB_plot <- reactive({
    text_listb = unlist(strsplit(input$genes_b,split="\n"))
    if (length(text_listb)==0) { uniq_listb =FALSE }else{
      uniq_listb = sapply(text_listb,function(x){sum(grepl(x,rownames(datasetInput())))==1})}
    if (sum(!uniq_listb)>0 | sum(uniq_listb)==0)
    {NULL} else{
      genes_to_plot_b = unique(rownames(datasetInput())[grepl(paste0(text_listb,collapse="|"),rownames(datasetInput()))]);
      genes_to_plot_b;}
  })

### gene list A text for plot
  formulaText_listA_plot <- reactive({
    text_lista = unlist(strsplit(input$genes_a,split="\n"))
    if (length(text_lista)==0) { uniq_lista =FALSE }else{
      uniq_lista = sapply(text_lista,function(x){sum(grepl(x,rownames(datasetInput())))==1})}
    if (sum(!uniq_lista)>0 | sum(uniq_lista)==0)
    {NULL} else{
      genes_to_plot_a = unique(rownames(datasetInput())[grepl(paste0(text_lista,collapse="|"),rownames(datasetInput()))]);
      genes_to_plot_a;}
  })
  
  
  
# create caption for gene list A
output$caption_listA <- renderText({
      formulaText_listA()
})

# create caption for gene list B
output$caption_listB <- renderText({
  formulaText_listB()
})


### scatterplot of norm. expression
output$expression_scatter <-renderPlot ({
  if (length(formulaText_listA_plot())==0 | length(formulaText_listB_plot())==0 ){  
    smer_shiny_plot_xy(
      bg_plot = TRUE)} else {
        all_n = smer_n[,colSums(smer_umi)>=input$min_xy*1e3]
        if (length(formulaText_listA_plot())==1){
          exp_x = all_n[formulaText_listA_plot(),]} else {
          exp_x = colSums(all_n[formulaText_listA_plot(),])}      
        if (length(formulaText_listB_plot())==1){
          exp_y = all_n[formulaText_listB_plot(),]} else {
          exp_y = colSums(all_n[formulaText_listB_plot(),])}        
        if (sum(input$scatter_features %in% "xy_log")>0){
          exp_x = log2(exp_x+1);
          exp_y = log2(exp_y+1);
        }

    
      smer_shiny_plot_xy(norm_mat = smer_n[,colSums(smer_umi)>=input$min_xy*1e3 & smer_md$genotype=="WT"],
                         bg_plot=FALSE,
                         rel_metadata = smer_md[colSums(smer_umi)>=input$min_xy*1e3 & smer_md$genotype=="WT",],
                         genes_a = formulaText_listA_plot(),
                         genes_b = formulaText_listB_plot(),
                         genes_a_label = input$genes_a_label,
                         genes_b_label = input$genes_b_label,
                         my_xlim= c(min(exp_x),max(exp_x)),
                         my_ylim= c(min(exp_y),max(exp_y)),
                         color_by_dev = TRUE,
                         color_na = sum(input$genotype_xy2 %in% "wt_xy")==0,
                         my_pch = 19,
                         my_cex=1,
                         #show_grid = sum(input$scatter_features %in% "xy_grid")>0,
						 show_grid = TRUE,
                         log2_exp = sum(input$scatter_features %in% "xy_log")>0,
                         reg_log=1)
						 
      if (sum(input$genotype_xy2 %in% "dst_xy")>0) {
        smer_shiny_add_points_xy(norm_mat = smer_n[,smer_md$total_umi>=input$min_xy*1e3 & smer_md$genotype=="dst"],
								 rel_metadata = smer_md[smer_md$total_umi>=input$min_xy*1e3 & smer_md$genotype=="dst",],
                                 genes_a = formulaText_listA_plot(),
                                 genes_b = formulaText_listB_plot(),
                                 color_by_dev = TRUE,
                                 my_cex=1,
                                 log2_exp = sum(input$scatter_features %in% "xy_log")>0,
								 reg_log = 1,
                                 #my_pch = 8,
								 my_pch = 19)
                                
      }
      
      if (sum(input$genotype_xy2 %in% "sft_xy")>0) {
        smer_shiny_add_points_xy(norm_mat = smer_n[,smer_md$total_umi>=input$min_xy*1e3 & smer_md$genotype=="sft"],
								rel_metadata = smer_md[smer_md$total_umi>=input$min_xy*1e3 & smer_md$genotype=="sft",],
                                 genes_a = formulaText_listA_plot(),
                                 genes_b = formulaText_listB_plot(),
                                 color_by_dev = TRUE,
                                 my_cex=1,
                                 log2_exp = sum(input$scatter_features %in% "xy_log")>0,
                                 #my_pch = 11,
								 my_pch=19,
								 reg_log = 1)
        
      }
      
      if (sum(input$genotype_xy2 %in% "dstsft_xy")>0) {
        smer_shiny_add_points_xy(norm_mat = smer_n[,smer_md$total_umi>=input$min_xy*1e3 & smer_md$genotype=="dst sft"],
								 rel_metadata = smer_md[smer_md$total_umi>=input$min_xy*1e3 & smer_md$genotype=="dst sft",],
                                 genes_a = formulaText_listA_plot(),
                                 genes_b = formulaText_listB_plot(),
                                 color_by_dev = TRUE,
                                 my_cex=1,
                                 log2_exp = sum(input$scatter_features %in% "xy_log")>0,
                                 #my_pch = 1,
								 my_pch=19,
								 reg_log = 1)
        
      }
      
      if (sum(input$genotype_xy2 %in% "uf_xy")>0) {
        smer_shiny_add_points_xy(norm_mat = smer_n[,smer_md$total_umi>=input$min_xy*1e3 & smer_md$genotype=="uf"],
								 rel_metadata = smer_md[smer_md$total_umi>=input$min_xy*1e3 & smer_md$genotype=="uf",],
                                 genes_a = formulaText_listA_plot(),
                                 genes_b = formulaText_listB_plot(),
                                 color_by_dev = TRUE,
                                 my_cex = 1,
                                 log2_exp = sum(input$scatter_features %in% "xy_log")>0,
                                 #my_pch = 17,
								 my_pch = 19,
								 reg_log = 1)
        
      }
  }
                         
                         
})

### scatterplot of norm. expression
output$barplot_gene1 <-renderPlot ({
  if (length(formulaText())==1){  
  covered_meristems = colnames(smer_umi)[colSums(smer_umi)>=(input$min_mer*1e3)]
  cat("generating barplots..")
    smer_shiny_gene_barplots(
		gene = formulaText(),
		raw_umi_matrix = smer_squeezed_meristems[,covered_meristems],
		metadata_table = smer_md[covered_meristems,],
		show_se = TRUE,
		order_by_genotype = input$order_bars_by_genotypes)
    }                     
                         
})

})
