require(shiny)
require(shinyWidgets)

# Define UI for dataset viewer application
	shinyUI(
#pageWithSidebar(
		fluidPage(
# Application title
			headerPanel(title = "The tomato Single-Meristem-Transcriptome (SMT) database	 v1.2", windowTitle = "the tomato SMT database v1.2"),
			img(src="mainPanel_image_labeled_big.png", align = "center"),

    
####################################
#### Dynamics page
##
#
  
tabsetPanel(
	tabPanel("Temporal expression dynamics",
  sidebarPanel(
  
  radioButtons(inputId = "assembly_version", 
                 choices=c("SL3.0, ITAG3.2" = "SL3",
			   "SL4.0, ITAG4.0" = "SL4"),
                 inline=TRUE,
                 label = "Assembly version:"),
     
  textInput(inputId = "gene",
			  label = "Enter gene name:",
			  value = "",
			  placeholder = "SolycXXgYYYYYY"),

			  
  radioButtons(inputId = "genotype_meristems1",
               choices = c("WT" = "WT",
                           "dst" = "dst",
                           "sft" = "sft",
                           "dst;sft" = "dst;sft",
                           "uf" = "uf"),
                inline = TRUE,
                label = "Genotype:"),

	conditionalPanel(
				condition = "output.flag_1 == true", 
					plotOutput("barplot_gene1",height = "300px"),
					
					  checkboxInput(inputId = "order_bars_by_genotypes",
				    			label = "Order bars by genotypes",
				  			value = FALSE),
					),
	
	br(),
	
	radioButtons(inputId = "genes_to_plot", 
                 choices=c("1" = "one", "2" = "two", "3" = "three"),
                 inline=TRUE,
                 label = "# Genes to compare:"),
				

	conditionalPanel(
					condition = "input.genes_to_plot=='two' || input.genes_to_plot=='three'", 
						textInput(inputId = "gene2",
						label="",
						value = "",
						placeholder = "SolycXXgYYYYYY")
					),

  conditionalPanel(
    condition = "input.genes_to_plot=='two' || input.genes_to_plot=='three'", 
    radioButtons(inputId = "genotype_meristems2",
                 choices = c("WT" = "WT",
                             "dst" = "dst",
                             "sft" = "sft",
                             "dst;sft" = "dst;sft",
                             "uf" = "uf"),
                 inline = TRUE,
                 label = "genotype")),
	
  conditionalPanel(
					condition = "input.genes_to_plot=='three'" ,
						textInput(inputId = "gene3",
						label="",
						value = "",
						placeholder = "SolycXXgYYYYYY")),
	
  conditionalPanel(
    condition = "input.genes_to_plot=='three'" ,
    radioButtons(inputId = "genotype_meristems3",
                 choices = c("WT" = "WT",
                             "dst" = "dst",
                             "sft" = "sft",
                             "dst;sft" = "dst;sft",
                             "uf" = "uf"),
                 inline = TRUE,
                 label = "genotype"))
  

  ),


################ show here ordered expresion of 1-3 genes  
mainPanel(
					 
	radioButtons(inputId = "exp_type",
	                        choices = c("KNN-norm." = "KNN-norm.",
									    "norm. counts" = "norm. counts",
									    "log-norm. counts" = "log-norm. counts",
	                                    "raw counts" = "raw counts"),
	                        selected = "KNN-norm.",
	                        label = "Expression values:",
				inline = TRUE
	   ),
	
	fluidRow(
	column(3, sliderInput(inputId = "min_mer",
			      post = "k UMIs",
	            	      label = "Min. meristem coverage:",
	            	      min=100,
	            	      max=1000,
	            	      value = 100)),
	
	column(3, sliderInput(inputId = "smooth_n",
				 label = "# Meristems to smooth:",
				 min=5,
				 max=50,
				 value = 13))
	),
	
	awesomeCheckboxGroup(inputId = "order_plot_feats",
						choices = c("Comparative mode" = "plot_bg",
							    "Show trend" = "trend_bool",
							    "Fixed scales" = "fixed_ylim"),
						label="",
						status = "primary",
						inline = TRUE),
								   				 
						conditionalPanel(
						condition = "output.flag_1 != true",
						h5(textOutput("caption_1"),align="center")),
						plotOutput("ordered_gene_Plot_1",height="250px"),
						conditionalPanel(
										condition = "(input.genes_to_plot=='two' || input.genes_to_plot=='three') && output.flag_2 != true",
										h5(textOutput("caption_2"),align="center")), 
						plotOutput("ordered_gene_Plot_2",height="250px"),			
						conditionalPanel(
										condition = "input.genes_to_plot=='three' && output.flag_3 != true",
										h5(textOutput("caption_3"),align="center")), 
						plotOutput("ordered_gene_Plot_3",height="250px")
						)
	),

######################################################
#### Co-expression page
##
#

	tabPanel("Co-expression",
		sidebarPanel(
			
		  textInput(inputId = "gene_corr",
			    label = "Insert gene name to compute correlation:",
			    value = "",
			    placeholder = "SolycXXgYYYYYY"),
			
			h5(textOutput("caption_corr"),
			   align="center"),
			br(),
			

			sliderInput(inputId = "up_corr",
			            label = "# coexpressed genes:",
			            min=1,
			            max=100,
			            value = 10),
			
			sliderInput(inputId = "down_corr",
			            label = "# 'anti' genes:",
			            min=1,
			            max=100,
			            value = 5),
			
			radioButtons(inputId = "cor_genotype", 
			             choices=c("WT" = "WT",
			                       "dst" = "dst", 
			                       "sft" = "sft"),
			             inline=TRUE,
			             label = "Choose genotype:"),
			
			conditionalPanel(
			  condition = "output.flag_corr == true",
		    tableOutput("corr_table"),
			  tableOutput("corr_table_anti"),
			  tags$head(tags$style("#corr_table{color: blue;
                                 font-size: 10px;
                                 font-style: italic;
                                 }"
			                      )
			            ),
			  tags$head(tags$style("#corr_table_anti{color: red;
                                 font-size: 10px;
                                 font-style: italic;
                                 }"
			                      )
			           )
			 )#conditionalPanel
			
						
			
			
	 ),#sidebarPanel
	 mainPanel(
		
		
		fluidRow(
	    column(6, plotOutput("expression_scatter",width = 400, height = 400), align = "center"),
			fluidRow(
			column(4, awesomeCheckboxGroup(inputId = "genotype_xy2",
										   choices=c("WT" = "wt_xy",
											     "dst" = "dst_xy",
											     "sft" = "sft_xy",
											     "dst;sft" = "dstsft_xy",
											     "uf" = "uf_xy"),
											selected="wt_xy",
											label = "Add meristems to plot:",
						       					status = "primary",
						       					inline = TRUE), 
				   align = "left"),
	   
			column(4, sliderInput(inputId = "min_xy",post = "k UMIs",
								 label = "Min. meristem coverage:",
								 min=100,
								 max=1000,
								 value = 100), 
				   align = "left"),
	   
			column(4,awesomeCheckboxGroup(inputId = "scatter_features",
										  choices=c(#"grid" = "xy_grid",
											    "log-norm" = "xy_log"),
										  selected=NULL,
										  label = "",status = "info",inline = TRUE),
				   align = "left")
		)),		
				
				
     	   
	   div(style="display: inline-block;vertical-align:top; width: 35%;",textInput(inputId="genes_a_label", label= "Enter label for gene list A:",value="",placeholder="gene list A")),
	   div(style="display: inline-block;vertical-align:top; width: 35%;",textInput(inputId="genes_b_label", label= "Enter label for gene list B:",value="", placeholder="gene list B")),
	   div(style="display: inline-block;vertical-align:top; width: 35%;",textAreaInput(inputId="genes_a", label= "",value="", rows=10,placeholder="SolycXXgYYYYYY")),
	   div(style="display: inline-block;vertical-align:top; width: 35%;",textAreaInput(inputId="genes_b", label= "",value="", rows=10,placeholder="SolycXXgYYYYYY")),
	   
	   div(style="display: inline-block;vertical-align:top; width: 35%;",h5(textOutput("caption_listA"))),
	   div(style="display: inline-block;vertical-align:top; width: 35%;",h5(textOutput("caption_listB")))
	 ) #mainPanel
	) #tabPanel
)
  
  
)

)
