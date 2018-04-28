library(shiny)

# TODO:
# really show RNA and meth together
# coordinates manually

### init
library(ExperimentHub)
library(ggplot2)
library(reshape2)
library(SummarizedExperiment)
library(bsseq)
library(rtracklayer)
eh <- ExperimentHub()
se_rpkms <- eh[["EH1074"]]
BS.obj.ex.fit <- eh[["EH1072"]]
# annotation
annotation_genome = "mm10" # mm10
mm10_refseq_file <- "/Users/charlie/Desktop/Michael/annotation/RefSeq_Sept_2013_from_annovar_Genes_plain.bed.gz"
mm10_cgi_file <- "/Users/charlie/Desktop/Michael/annotation/mm10.CGI.bed"
mm10_promotor_file <- "/Users/charlie/Desktop/Michael/annotation/mm10_RefSeq_UPSTREAM2000DOWNSTREAM500.csv"
mm10_enhancer_file <- "/Users/charlie/Desktop/Michael/annotation/all_enhancer_mm10.bed"
if (annotation_genome == "mm10") {
  refseq <- import(mm10_refseq_file, format="bed")
  cgi <- import(mm10_cgi_file, format="bed")
  seqlevels(cgi)<- sub('chr','',seqlevels(cgi))
  promoter <- import(mm10_promotor_file, format="bed")
  seqlevels(promoter)<- sub('chr','',seqlevels(promoter))
  enhancer <- import(mm10_enhancer_file, format="bed")
  seqlevels(enhancer)<- sub('chr','',seqlevels(enhancer))
  annoTrack=c(GENE=refseq, CGI=cgi, PROM=promoter, ENHAN=enhancer)
}
sample.colors <- rep(c("red", "blue", "yellow", "green", "orange"), rep(3,5))
names(sample.colors) <- c("Fat-Treg-R1", "Fat-Treg-R2", "Fat-Treg-R3",
                         "Liver-Treg-R1", "Liver-Treg-R2", "Liver-Treg-R3",
                         "Skin-Treg-R1", "Skin-Treg-R2", "Skin-Treg-R3",
                         "Lymph-N-Tcon-R1", "Lymph-N-Tcon-R2", "Lymph-N-Tcon-R3",
                         "Lymph-N-Treg-R1", "Lymph-N-Treg-R2", "Lymph-N-Treg-R3")
group.colors <- rep(c("red", "blue", "yellow", "green", "orange"), rep(3,5))
names(group.colors) <- rep(c("Fat-Treg", "Liver-Treg", "Skin-Treg", "Lymph-N-Tcon", "Lymph-N-Treg"), rep(3,5))

# Define UI for application that draws a histogram
ui <- fluidPage(
   sidebarLayout(
      sidebarPanel(
        textInput("gene", "Gene", "Foxp3"),
        verbatimTextOutput("value")
      ),
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("methylation"),
         plotOutput("expression")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$expression <- renderPlot({
    foxp3_rpkm <- assay(se_rpkms)[rowData(se_rpkms)$id_symbol == input$gene,]
    foxp3_rpkm_molten <- melt(foxp3_rpkm)
    # actual plot
    ggplot(data=foxp3_rpkm_molten, aes(x=rownames(foxp3_rpkm_molten), y=value, fill=colData(se_rpkms)$tissue_cell)) +
      geom_bar(stat="identity") +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      xlab("Sample") +
      ylab("RPKM") +
      ggtitle(input$caption) +
      guides(fill=guide_legend(title="tissue / cell type")) +
      scale_fill_manual(values = group.colors)
  })
   output$methylation <- renderPlot({
      resolveGenesymbolToEntrez <- function(genesysmbol = NA) {
        library(biomaRt)
        mart <- useMart("ensembl")
        mart <- useDataset("mmusculus_gene_ensembl", mart = mart)
        genes.table <- getBM(filters= "external_gene_name",
                             attributes= c("entrezgene", "chromosome_name", "start_position", "end_position",
                                           "external_gene_name", "description"), values=genesysmbol , mart=mart)
        genesysmbol <- data.frame(external_gene_name = genesysmbol)
        results <- merge(x=as.data.frame(genesysmbol), y=genes.table, by.x="external_gene_name", by.y="external_gene_name", all.x = T)
        return(results)
      }
      genes_of_interest_anno <- resolveGenesymbolToEntrez(input$gene)
      regions <- GRanges(
        seqnames = genes_of_interest_anno$chromosome_name,
        ranges = IRanges(start = genes_of_interest_anno$start_position,
                         end = genes_of_interest_anno$end_position
        )
      )
      pData <- pData(BS.obj.ex.fit)
      pData$col <- unname(sample.colors)
      pData(BS.obj.ex.fit) <- pData
      # samples_subset <- c("Lymph-N-Tcon-R1", "Lymph-N-Tcon-R2", "Lymph-N-Tcon-R3", "Lymph-N-Treg-R1", "Lymph-N-Treg-R2", "Lymph-N-Treg-R3")
      i <- 1
      plotRegion(BS.obj.ex.fit, region=regions[i,],
                 extend = 10000,
                 main=genes_of_interest_anno[i,]$external_gene_name,
                 addRegions = regions,
                 annoTrack = annoTrack)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

