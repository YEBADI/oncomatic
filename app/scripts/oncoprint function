#' Creating a oncoprint
#' @param mut A dataframe from the mutation annotation file (see TCGAquery_maf from TCGAbiolinks)
#' @param genes Gene list
#' @param filename name of the pdf
#' @param color named vector for the plot
#' @param height pdf height
#' @param width pdf width
#' @param rm.empty.columns If there is no alteration in that sample, whether remove it on the oncoprint
#' @param show.row.barplot  Show barplot annotation on rows?
#' @param show.column.names Show column names? Default: FALSE
#' @param rows.font.size Size of the fonts
#' @param column.names.size Size of the fonts of the columns names
#' @param dist.col distance between columns in the plot
#' @param dist.row distance between rows in the plot
#' @param label.font.size Size of the fonts
#' @param row.order Order the genes (rows) Default:TRUE. Genes with more mutations will be in the first rows
#' @param col.order Order columns. Default:TRUE.
#' @param annotation Matrix or data frame with the annotation.
#' Should have a column bcr_patient_barcode with the same ID of the mutation object
#' @param annotation.position Position of the annotation "bottom" or "top"
#' @param label.title Title of the label
#' @param annotation.legend.side Position of the annotation legend
#' @param heatmap.legend.side Position of the heatmap legend
#' @param information Which column to use as informastion from MAF.
#' Options: 1) "Variant_Classification" (The information will be "Frame_Shift_Del", "Frame_Shift_Ins",
#'         "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",  "Nonsense_Mutation",
#'              "Nonstop_Mutation",  "RNA",  "Silent" ,  "Splice_Site",  "Targeted_Region",  "Translation_Start_Site")
#' 2) "Variant_Type" (The information will be INS,DEL,SNP)
#' @importFrom ComplexHeatmap oncoPrint draw HeatmapAnnotation
#' @importFrom grid gpar grid.rect
#' @importFrom data.table dcast setDT setDF :=
#' @examples
#' mut <- GDCquery_Maf(tumor = "ACC", pipelines = "mutect")
#' TCGAvisualize_oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:10], rm.empty.columns = TRUE)
#' TCGAvisualize_oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:10],
#'                  filename = "onco.pdf",
#'                  color=c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"))
#' clin <- GDCquery_clinic("TCGA-ACC","clinical")
#' clin <- clin[,c("bcr_patient_barcode","disease","gender","tumor_stage","race","vital_status")]
#' TCGAvisualize_oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:20],
#'                 filename = "onco.pdf",
#'                 annotation = clin,
#'                 color=c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"),
#'                 rows.font.size=10,
#'                 heatmap.legend.side = "right",
#'                 dist.col = 0,
#'                 label.font.size = 10)
#'
#' @export
#' @return A oncoprint plot
TCGAvisualize_oncoprint <- function (mut,
                                     genes,
                                     filename,
                                     color,
                                     annotation.position = "bottom",
                                     annotation,
                                     height,
                                     width = 10,
                                     rm.empty.columns = FALSE,
                                     show.column.names = FALSE,
                                     show.row.barplot = TRUE,
                                     label.title = "Mutation",
                                     column.names.size = 8,
                                     label.font.size = 16,
                                     rows.font.size = 16,
                                     dist.col = 0.5,
                                     dist.row = 0.5,
                                     information = "Variant_Type",
                                     row.order = TRUE,
                                     col.order = TRUE,
                                     heatmap.legend.side = "bottom",
                                     annotation.legend.side = "bottom"){


    if(missing(mut))   stop("Missing mut argument")
    mut <- setDT(mut)
    mut$value <- 1
    if(rm.empty.columns == FALSE) all.samples <- unique(mut$Tumor_Sample_Barcode)

    mut$Hugo_Symbol <- as.character(mut$Hugo_Symbol)
    if(!missing(genes) & !is.null(genes)) mut <- subset(mut, mut$Hugo_Symbol %in% genes)

    if(!rm.empty.columns){
        formula <- paste0("Tumor_Sample_Barcode + Hugo_Symbol ~ ", information)
        suppressMessages({mat <- dcast(mut, as.formula(formula),value.var = "value",fill = 0,drop = FALSE)})
    } else {
        formula <- paste0("Tumor_Sample_Barcode + Hugo_Symbol ~ ", information)
        suppressMessages({mat <- dcast(mut, as.formula(formula),value.var = "value",fill = 0,drop = TRUE)})
    }

    # mutation in the file
    columns <- colnames(mat)[-c(1:2)]

    # value will be a collum with all the mutations
    mat$value <- ""

    for ( i in columns){
        mat[,i] <-  replace(mat[,i,with = FALSE],mat[,i,with = FALSE]>0,paste0(i,";"))
        mat[,i] <-  replace(mat[,i,with = FALSE],mat[,i,with = FALSE]==0,"")
        mat[,value:=paste0(value,get(i))]
    }

    # After the gene selection, some of the mutation might not exist
    # we will remove them to make the oncoprint work
    mutation.type <- c()
    for (i in columns){
        if(length(grep(i,mat$value)) > 0) mutation.type <- c(mutation.type,i)
    }

    # now we have a matrix with pairs samples/genes mutations
    # we want a matrix with samples vs genes mutations with the content being the value
    mat <- setDF(dcast(mat, Tumor_Sample_Barcode~Hugo_Symbol, value.var="value",fill=""))
    rownames(mat) <- mat[,1]
    mat <- mat[,-1]

    if(rm.empty.columns == FALSE) {
        aux <- data.frame(row.names = all.samples[!all.samples %in% rownames(mat)])
        if(nrow(aux) > 0) {
            aux[,colnames(mat)] <- ""
            mat <- rbind(mat,aux)
        }
    }


    alter_fun = function(x, y, w, h, v) {
        n = sum(v)
        h = h*0.9
        # use `names(which(v))` to correctly map between `v` and `col`
        if(n) {
            grid.rect(x, y - h*0.5 + 1:n/n*h,  w-unit(dist.col, "mm"), 1/n*h,
                      gp = gpar(fill = color[names(which(v))], col = NA), just = "top")
        } else {
            grid.rect(x, y, w-unit(dist.col, "mm"), h-unit(dist.row, "mm"), gp = gpar(fill = color["background"], col = NA))
        }
    }

    # get only the colors to the mutations
    # otherwise it gives errors

    if(missing(color)){
        color <- c(rainbow(length(mutation.type)), "#CCCCCC")
        names(color) <- c(mutation.type,"background")
    } else{
        if("background" %in% names(color)) {
            color <- color[c(mutation.type,"background")]
        } else {
            color <- c(color[mutation.type],"background"= "#CCCCCC")
        }
    }
    # header are samples, rows genes
    mat <- t(mat)

    if(!missing(height)) height <- length(genes)/2
    if(!missing(filename)) pdf(filename,width = width,height = height)

    if(missing(annotation)) annotation <- NULL
    if(!is.null(annotation)){
        if(!"bcr_patient_barcode" %in% colnames(annotation))
            stop("bcr_patient_barcode column should be in the annotation")
        idx <- match(substr(colnames(mat),1,12),annotation$bcr_patient_barcode)
        if(all(is.na(idx)))
            stop(" We couldn't match the columns names with the bcr_patient_barcode column in the annotation object")
        annotation <- annotation[idx,]

        annotation$bcr_patient_barcode <- NULL

        n.col <- sum(sapply(colnames(annotation), function(x) {
            length(unique(annotation[,x]))
        }))

        # add automatic colors: not working

        get.color <- function(df,col){
            idx <- which(colnames(df) == col)
            start <- 1
            if(idx != 1) start <- length(na.omit(unique(unlist(c(df[,1:(idx-1)]))))) + 1
            end <- start + length(na.omit(unique(df[,col]))) -1
            diff.colors <- c("purple","thistle","deeppink3","magenta4","lightsteelblue1","black",
                             "chartreuse","lightgreen","maroon4","darkslategray",
                             "lightyellow3","darkslateblue","firebrick1","aquamarine",
                             "dodgerblue4","bisque4","moccasin","indianred1",
                             "yellow","gray93","cyan","darkseagreen4",
                             "lightgoldenrodyellow","lightpink","sienna1",
                             "darkred","palevioletred","tomato4","blue",
                             "mediumorchid4","royalblue1","magenta2","darkgoldenrod1")
            return(diff.colors[start:end])
        }
        col.annot <- lapply(colnames(annotation), function(x) {
            #idx <- which(colnames(annotation) == x) - 1
            #print(idx/n.col)
            ret <- get.color(annotation,x)
            #ret <- rainbow(length(unique(annotation[,x])),start = idx/n.col,alpha=0.5)
            names(ret) <- as.character(na.omit(unique(annotation[,x])))
            return(ret)
        })
        names(col.annot) <-  colnames(annotation)

        annotHeatmap <- HeatmapAnnotation(df=annotation,
                                          col=col.annot,
                                          annotation_legend_param=list(title_gp=gpar(fontsize=label.font.size,
                                                                                     fontface="bold"),
                                                                       labels_gp=gpar(fontsize=label.font.size),#sizelabels
                                                                       grid_height=unit(8,"mm"))
        )
    }
    if(heatmap.legend.side == "bottom") {
        nrow <- 1
        title_position <- "leftcenter"
    } else {
        nrow <- 10
        title_position <- "topcenter"
    }

    if(is.null(annotation) & !row.order & !col.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       column_names_gp = gpar(fontsize = column.names.size),
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    } else if(!is.null(annotation) & annotation.position == "bottom" & !row.order & !col.order){

        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       column_names_gp = gpar(fontsize = column.names.size),
                       show_row_barplot = show.row.barplot,
                       show_column_names = show.column.names,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       bottom_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )

    } else if(!is.null(annotation) & annotation.position == "top" & !row.order  & !col.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       column_names_gp = gpar(fontsize = column.names.size),
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       top_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"),  # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    }  else if(is.null(annotation) & row.order  & !col.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       column_names_gp = gpar(fontsize = column.names.size),
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    } else if(!is.null(annotation) & annotation.position == "bottom" & row.order & !col.order){

        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       remove_empty_columns = FALSE,
                       column_names_gp = gpar(fontsize = column.names.size),
                       show_row_barplot = show.row.barplot,
                       show_column_names = show.column.names,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       bottom_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )

    } else if(!is.null(annotation) & annotation.position == "top" & row.order & !col.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       column_names_gp = gpar(fontsize = column.names.size),
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       top_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"),  # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    } else if(is.null(annotation) & !row.order & col.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       alter_fun = alter_fun, col = color,
                       column_names_gp = gpar(fontsize = column.names.size),
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    } else if(!is.null(annotation) & annotation.position == "bottom" & !row.order & col.order){

        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       column_names_gp = gpar(fontsize = column.names.size),
                       show_row_barplot = show.row.barplot,
                       show_column_names = show.column.names,
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       bottom_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )

    } else if(!is.null(annotation) & annotation.position == "top" & !row.order & col.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       column_names_gp = gpar(fontsize = column.names.size),
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       top_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"),  # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    }  else if(is.null(annotation) & row.order & col.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       alter_fun = alter_fun, col = color,
                       column_names_gp = gpar(fontsize = column.names.size),
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    } else if(!is.null(annotation) & annotation.position == "bottom" & row.order & col.order){

        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       remove_empty_columns = FALSE,
                       column_names_gp = gpar(fontsize = column.names.size),
                       show_row_barplot = show.row.barplot,
                       show_column_names = show.column.names,
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       bottom_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )

    } else if(!is.null(annotation) & annotation.position == "top" & row.order & col.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       alter_fun = alter_fun, col = color,
                       column_names_gp = gpar(fontsize = column.names.size),
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       top_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"),  # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    }

    draw(p, heatmap_legend_side = heatmap.legend.side, annotation_legend_side = annotation.legend.side)
    if(!missing(filename)) {
        dev.off()
        message(paste0("File saved as: ", filename ))
    }

}
