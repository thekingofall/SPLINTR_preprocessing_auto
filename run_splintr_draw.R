#!/usr/bin/env Rscript

# --- Clean Environment ---
rm(list = ls())

# --- Check and Install Packages ---
required_packages <- c("ggplot2", "data.table", "dplyr", "assertthat",
                       "seqinr", "fastmatch", #"venn", # Not strictly needed
                       "VennDiagram", "patchwork", "RColorBrewer",
                       "futile.logger", "grid", "scales", "UpSetR",
                       "argparse") # Added argparse

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, repos = "http://cran.us.r-project.org") # Specify repo for non-interactive install
  }
}
cat("Package check complete.\n\n")

# --- Load Libraries ---
cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(assertthat)
  library(seqinr)
  library(fastmatch)
  library(VennDiagram)
  library(patchwork)
  library(RColorBrewer)
  library(futile.logger)
  library(grid)
  library(scales)
  library(UpSetR)
  library(argparse)
})
cat("Libraries loaded.\n\n")

# --- Suppress VennDiagram Log ---
flog.threshold(ERROR, name = "VennDiagramLogger")

# --- Define Command Line Arguments ---
parser <- ArgumentParser(description = "SPLINTR Library Barcode Preprocessing and Analysis (Multiple Samples)")

parser$add_argument("--input_dir", type = "character",
                    default = "/home/YangZongmian/SPLINTR/X101SC24127971-Z01-J027/Processed_SPLINTR_Output_Parallel/Step2_Starcode",
                    help = "Directory containing the '*_starcode.txt' files from Step 2.")

parser$add_argument("--output_dir", type = "character",
                    default = "/home/YangZongmian/SPLINTR/X101SC24127971-Z01-J027/Processed_SPLINTR_Output_Parallel/Step3_SPLINTR_Analysis_Output",
                    help = "Main output directory for this analysis.")

parser$add_argument("--threshold", type = "integer",
                    default = 5,
                    help = "Minimum raw count threshold for filtering barcodes.")

# --- Parse Arguments ---
args <- parser$parse_args()

# --- Assign Parsed Arguments to Variables ---
starcode_input_dir <- args$input_dir
analysis_output_dir <- args$output_dir
raw_count_threshold <- args$threshold # This will be available globally for functions that check its existence

# --- Setup Environment (Using Parsed Arguments) ---

# Create output sub-directories
individual_plots_dir <- file.path(analysis_output_dir, "individual_plots")
results_dir <- file.path(analysis_output_dir, "results")
combined_plots_dir <- file.path(analysis_output_dir, "combined_plots")

cat("Creating output directories if they don't exist...\n")
dir.create(individual_plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(combined_plots_dir, showWarnings = FALSE, recursive = TRUE)

# Print configuration
cat("------------------------------------\nConfiguration:\n------------------------------------\n")
cat("Starcode Input directory:", starcode_input_dir, "\n")
cat("Analysis Output directory:", analysis_output_dir, "\n")
cat("  -> Individual Plots directory:", individual_plots_dir, "\n")
cat("  -> Results directory:", results_dir, "\n")
cat("  -> Combined Plots directory:", combined_plots_dir, "\n")
cat("Raw Count Threshold:", raw_count_threshold, "\n")
cat("------------------------------------\n\n")

# --- Load Input Data ---
cat("--- Loading Input Data ---\n")
starcode_files <- list.files(path = starcode_input_dir, pattern = "Sample.*_starcode\\.txt$", full.names = TRUE, recursive = FALSE)
if (length(starcode_files) == 0) {
  stop("FATAL ERROR: No 'Sample*_starcode.txt' files found in: ", starcode_input_dir)
}
cat("Found", length(starcode_files), "starcode input files:\n")
print(basename(starcode_files))

starcode_data_list <- list()
read_errors <- list()
for (f in starcode_files) {
  s_name <- gsub("_starcode\\.txt$", "", basename(f))
  tryCatch({
    file_info <- file.info(f)
    if (is.na(file_info$size) || file_info$size < 2) {
      warning("Skipping potentially empty/invalid file: ", basename(f))
      read_errors[[s_name]] <- "File empty or too small"
    } else {
      dt <- fread(f, col.names = c("Barcode", "Raw_count"), header = FALSE, sep = "\t")
      if (nrow(dt) > 0 && ncol(dt) == 2) {
        dt$Raw_count <- as.numeric(dt$Raw_count)
        dt <- dt[!is.na(dt$Raw_count), ]
        if(nrow(dt) > 0) starcode_data_list[[s_name]] <- dt
        else { read_errors[[s_name]] <- "No valid numeric data rows"; warning(paste(s_name,":", read_errors[[s_name]])) }
      } else if (nrow(dt) == 0) {
        read_errors[[s_name]] <- "No data rows in file"; warning(paste(s_name,":", read_errors[[s_name]]))
      } else {
        read_errors[[s_name]] <- "Incorrect column count"; warning(paste(s_name,":", read_errors[[s_name]]))
      }
    }
  }, error = function(e) {
    warning("Error reading file ", basename(f), ": ", e$message); read_errors[[s_name]] <- e$message
  })
}

sample_names <- names(starcode_data_list)
if (length(starcode_data_list) == 0) {
  stop("FATAL ERROR: No data could be loaded. Check input files and read errors.")
}
cat("\nSuccessfully loaded data for", length(sample_names), "samples:", paste(sample_names, collapse=", "), "\n")
if (length(read_errors) > 0) { cat("\nErrors/Warnings during loading:\n"); print(read_errors) }
cat("--- Input Data Loaded ---\n\n")


#----- Define Common Functions -----#

# Barcode frequency distribution plot function (with axis limits)
plot.barcode.dist <- function(barcodes_df, plot_title, output_filename, xlims = NULL, ylims = NULL){
  barcodes_df$Raw_count <- as.numeric(barcodes_df$Raw_count)
  barcodes_df_plot <- barcodes_df[barcodes_df$Raw_count > 0 & !is.na(barcodes_df$Raw_count), ]
  if (nrow(barcodes_df_plot) == 0) { warning("No valid barcodes for frequency plot: ", plot_title); return(NULL) }
  barcodes_df_plot <- barcodes_df_plot[order(barcodes_df_plot$Raw_count, decreasing = TRUE), ]
  barcodes_df_plot$Rank <- seq_len(nrow(barcodes_df_plot))
  
  p <- ggplot(barcodes_df_plot, aes(y = Raw_count, x = Rank)) +
    geom_point(stat = "identity", show.legend = FALSE, size = 0.5, alpha = 0.7) +
    scale_y_continuous(trans = 'log10', labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_cartesian(xlim = xlims, ylim = ylims, expand = FALSE) + # Use coord_cartesian
    # Use the globally defined raw_count_threshold from parsed args
    {if(exists("raw_count_threshold")) geom_hline(yintercept = raw_count_threshold, color = "blue", linetype="dashed", alpha=0.6) else NULL} +
    theme_bw(base_size = 10) + xlab("Barcode Rank") + ylab("Raw Count (log10 scale)") + ggtitle(plot_title) +
    theme(plot.title = element_text(size = 11, face = "bold"), axis.text = element_text(size=9), axis.title = element_text(size=10))
  
  suppressMessages(suppressWarnings( ggsave(output_filename, plot = p, dpi = 300, width = 6, height = 4, units = "in") ))
  return(p)
}

# Barcode cumulative sum plot function (fixed 0-1 axes)
plot.barcode.cumsum <- function(barcodes_df, plot_title, output_filename){
  barcodes_df$Raw_count <- as.numeric(barcodes_df$Raw_count)
  sorted_counts <- sort(barcodes_df$Raw_count[barcodes_df$Raw_count > 0 & !is.na(barcodes_df$Raw_count)], decreasing = TRUE)
  if (length(sorted_counts) == 0) { warning("No valid barcodes for cumulative sum plot: ", plot_title); return(NULL) }
  cumsum_val <- cumsum(as.numeric(sorted_counts)); max_cumsum <- max(cumsum_val)
  if (max_cumsum == 0) { warning("Total count is zero for cumulative sum: ", plot_title); return(NULL) }
  cumsum_prop <- cumsum_val / max_cumsum; barcode_rank_prop <- seq_along(cumsum_prop) / length(cumsum_prop)
  plot_data <- data.frame(barcode_rank_prop = barcode_rank_prop, proportion = cumsum_prop)
  
  p <- ggplot(plot_data, aes(y = proportion, x = barcode_rank_prop)) +
    geom_line(show.legend = FALSE, size=0.8) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + # Explicit 0-1 limits
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_bw(base_size = 10) + xlab("Barcode Rank (Proportion)") + ylab("Cumulative Read Sum (Proportion)") + ggtitle(plot_title) +
    theme(plot.title = element_text(size = 11, face = "bold"), axis.text = element_text(size=9), axis.title = element_text(size=10))
  
  suppressMessages(suppressWarnings( ggsave(output_filename, plot = p, dpi = 300, width = 6, height = 4, units = "in") ))
  return(p)
}

# Function to parse starcode output, filter, plot (with global limits), and save results
# Takes count_threshold explicitly
parse.starcode <- function(starcode_data, sample_id, count_threshold, individual_plot_dir, results_out_dir,
                           global_raw_xlims = NULL, global_raw_ylims = NULL, global_filt_xlims = NULL, global_filt_ylims = NULL){
  barcodes <- as.data.frame(starcode_data); colnames(barcodes) <- c("Barcode", "Raw_count"); barcodes$Raw_count <- as.numeric(barcodes$Raw_count)
  cat("Processing sample:", sample_id, "- Initial barcode count:", nrow(barcodes), "\n")
  p_raw_dist <- NULL; p_raw_cumsum <- NULL; p_filt_dist <- NULL; p_filt_cumsum <- NULL
  
  plot_title_raw <- paste(sample_id, ": Raw Frequency"); plot_filename_raw_dist <- file.path(individual_plot_dir, paste0(sample_id, "_01_raw_barcode_frequency.png"))
  p_raw_dist <- plot.barcode.dist(barcodes, plot_title_raw, plot_filename_raw_dist, xlims = global_raw_xlims, ylims = global_raw_ylims)
  plot_title_raw_cumsum <- paste(sample_id, ": Raw Cumulative Sum"); plot_filename_raw_cumsum <- file.path(individual_plot_dir, paste0(sample_id, "_02_raw_cumulative_sum.png"))
  p_raw_cumsum <- plot.barcode.cumsum(barcodes, plot_title_raw_cumsum, plot_filename_raw_cumsum)
  
  barcodes_filtered <- barcodes[which(barcodes$Raw_count >= count_threshold),]; n_filtered <- nrow(barcodes_filtered)
  cat("Sample:", sample_id, "- Filtered count (>=", count_threshold, "):", n_filtered, "\n")
  if (n_filtered == 0) {
    warning("No barcodes passed filter for sample: ", sample_id)
    return(list(filtered_df = NULL, plots = list(raw_dist=p_raw_dist, raw_cumsum=p_raw_cumsum, filt_dist=NULL, filt_cumsum=NULL)))
  }
  
  barcodes_filtered <- barcodes_filtered[order(barcodes_filtered$Raw_count, decreasing = TRUE), ]
  rank_prefix <- gsub("^Sample", "", sample_id); if (nchar(rank_prefix) == 0 || rank_prefix == sample_id) rank_prefix <- sample_id
  ranks <- paste0(rank_prefix, "_Barcode_", seq_len(n_filtered)); barcodes_filtered$Rank <- ranks
  
  plot_title_filt <- paste(sample_id, ": Filtered Freq (>=", count_threshold, ")"); plot_filename_filt_dist <- file.path(individual_plot_dir, paste0(sample_id, "_03_filtered_barcode_frequency.png"))
  p_filt_dist <- plot.barcode.dist(barcodes_filtered, plot_title_filt, plot_filename_filt_dist, xlims = global_filt_xlims, ylims = global_filt_ylims)
  plot_title_filt_cumsum <- paste(sample_id, ": Filtered CumSum (>=", count_threshold, ")"); plot_filename_filt_cumsum <- file.path(individual_plot_dir, paste0(sample_id, "_04_filtered_cumulative_sum.png"))
  p_filt_cumsum <- plot.barcode.cumsum(barcodes_filtered, plot_title_filt_cumsum, plot_filename_filt_cumsum)
  
  fasta_filename <- file.path(results_out_dir, paste0(sample_id, "_filtered_barcodes.fasta"))
  tryCatch({ seqinr::write.fasta(sequences = as.list(barcodes_filtered$Barcode), names = barcodes_filtered$Rank, file.out = fasta_filename, open = "w", as.string = TRUE); cat("Saved filtered FASTA to:", fasta_filename, "\n") }, error = function(e) warning("Error writing FASTA for ", sample_id, ": ", e$message))
  table_filename <- file.path(results_out_dir, paste0(sample_id, "_filtered_barcodes.tsv"))
  output_table <- barcodes_filtered[, c("Rank", "Barcode", "Raw_count")]
  tryCatch({ write.table(output_table, table_filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE); cat("Saved filtered table to:", table_filename, "\n") }, error = function(e) warning("Error writing TSV for ", sample_id, ": ", e$message))
  
  return(list(filtered_df = barcodes_filtered, plots = list(raw_dist = p_raw_dist, raw_cumsum = p_raw_cumsum, filt_dist = p_filt_dist, filt_cumsum = p_filt_cumsum)))
}

# Function to combine plots using patchwork (converts to grobs)
save_combined_plot <- function(plot_list, filename_base, p_dir, title, ncols = 2) {
  plot_list <- plot_list[!sapply(plot_list, is.null)]; n_plots <- length(plot_list)
  if (n_plots == 0) { cat("No valid plots to combine for:", title, "\n"); return() }
  cat("Combining", n_plots, "plots for:", title, "\n")
  plot_height_per_row <- 4; plot_width_per_col <- 6; nrows <- ceiling(n_plots / ncols)
  plot_height_combined <- plot_height_per_row * nrows; plot_width_combined <- plot_width_per_col * ncols
  grob_list <- list(); conversion_success <- TRUE
  for(i in seq_along(plot_list)) {
    plot_name <- names(plot_list)[i]; if(is.null(plot_name)) plot_name <- paste("Plot", i)
    grob_obj <- tryCatch(ggplotGrob(plot_list[[i]]), error = function(e) { cat("ERROR converting plot '", plot_name, "' to grob. Error:", e$message, "\n"); NULL })
    if (is.null(grob_obj)) conversion_success <- FALSE else grob_list[[i]] <- grob_obj
  }
  if (!conversion_success || length(grob_list) == 0) { cat("Skipping combination for '", title, "' due to grob conversion errors.\n"); return() }
  if(conversion_success) names(grob_list) <- names(plot_list)
  
  combined_plot_obj <- tryCatch({ patchwork::wrap_plots(grob_list, ncol = ncols, axes = 'keep', guides = 'collect') + plot_annotation(title = title, theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))) }, error = function(e) {
    cat("ERROR during patchwork combination for '", title, "'. Error:", e$message, "\n"); tryCatch(patchwork::wrap_plots(grob_list, ncol = ncols, axes = 'keep', guides = 'collect'), error = function(e2) { cat("Fallback wrap_plots failed. Error:", e2$message, "\n"); NULL }) })
  if (!is.null(combined_plot_obj)) {
    output_filename <- file.path(p_dir, paste0(filename_base, ".png"))
    suppressMessages(suppressWarnings( ggsave(filename = output_filename, plot = combined_plot_obj, width = plot_width_combined, height = plot_height_combined, dpi = 300, limitsize = FALSE, device = 'png') ))
    cat("Saved combined plot to:", output_filename, "\n")
  } else { cat("Failed to generate/save combined plot for:", title, "\n") }
}

# Function for Overlaid log10 count density plot
plot.overlaid.density <- function(filt_data_list, plot_title, output_filename) {
  if (length(filt_data_list) == 0) { cat("No filtered data for overlaid density.\n"); return(NULL) }
  combined_df <- bind_rows(filt_data_list, .id = "Sample") %>% filter(Raw_count > 0 & !is.na(Raw_count))
  if (nrow(combined_df) == 0) { cat("No positive counts for density plot.\n"); return(NULL) }
  p <- ggplot(combined_df, aes(x = log10(Raw_count), color = Sample)) + geom_density(show.legend = TRUE, size=0.8) +
    theme_bw(base_size = 10) + ggtitle(plot_title) + xlab("Raw Count (log10 scale)") + ylab("Density") +
    theme(plot.title = element_text(size = 12, face = "bold"), legend.position = "right") + scale_color_brewer(palette = "Paired")
  ggsave(output_filename, plot = p, dpi = 300, width = 7, height = 5, units = "in"); cat("Saved overlaid density plot to:", output_filename, "\n"); return(p)
}

# Function for Overlaid cumulative sum plot
plot.overlaid.cumsum <- function(filt_data_list, plot_title, output_filename) {
  if (length(filt_data_list) == 0) { cat("No filtered data for overlaid cumsum.\n"); return(NULL) }
  all_cumsum_data <- list()
  for (s_name in names(filt_data_list)) {
    df <- filt_data_list[[s_name]]; df$Raw_count <- as.numeric(df$Raw_count)
    sorted_counts <- sort(df$Raw_count[df$Raw_count > 0 & !is.na(df$Raw_count)], decreasing = TRUE)
    if (length(sorted_counts) > 0) {
      cumsum_val <- cumsum(as.numeric(sorted_counts)); max_cumsum <- max(cumsum_val)
      if (max_cumsum > 0) {
        cumsum_prop <- cumsum_val / max_cumsum; barcode_rank_prop <- seq_along(cumsum_prop) / length(cumsum_prop)
        all_cumsum_data[[s_name]] <- data.frame(Sample = s_name, barcode_rank_prop = barcode_rank_prop, proportion = cumsum_prop) } } }
  if (length(all_cumsum_data) == 0) { cat("Could not calculate cumulative data.\n"); return(NULL) }
  combined_cumsum_df <- bind_rows(all_cumsum_data)
  p <- ggplot(combined_cumsum_df, aes(x = barcode_rank_prop, y = proportion, color = Sample)) + geom_line(show.legend = TRUE, size=0.8) +
    theme_bw(base_size = 10) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    xlab("Barcode Rank (Proportion)") + ylab("Cumulative Read Sum (Proportion)") + ggtitle(plot_title) +
    theme(plot.title = element_text(size = 12, face = "bold"), legend.position = "right") + scale_color_brewer(palette = "Paired")
  ggsave(output_filename, plot = p, dpi = 300, width = 7, height = 5, units = "in"); cat("Saved overlaid cumulative plot to:", output_filename, "\n"); return(p)
}

#------- Analysis Pipeline -------#

# --- Calculate Global Axis Limits --- #
cat("--- Calculating Global Axis Limits ---\n")
global_max_raw_rank <- 0; global_max_raw_count <- 0
global_max_filt_rank <- 0; global_max_filt_count <- 0
temp_filtered_counts <- list()

cat("Pass 1: Calculating raw limits and filtering...\n")
for (s_name in names(starcode_data_list)) {
  df <- starcode_data_list[[s_name]]; df_valid <- df[df$Raw_count > 0 & !is.na(df$Raw_count), ]
  if(nrow(df_valid) > 0) { global_max_raw_rank <- max(global_max_raw_rank, nrow(df_valid)); global_max_raw_count <- max(global_max_raw_count, df_valid$Raw_count) }
  # Use the globally defined threshold from args
  df_filt <- df[which(df$Raw_count >= raw_count_threshold),]; if(nrow(df_filt) > 0) { temp_filtered_counts[[s_name]] <- df_filt$Raw_count } }

cat("Pass 2: Calculating filtered limits...\n")
if (length(temp_filtered_counts) > 0) { global_max_filt_rank <- max(sapply(temp_filtered_counts, length)); global_max_filt_count <- max(unlist(temp_filtered_counts))
} else { cat("Warning: No samples had data passing the filter threshold.\n") }

axis_buffer <- 1.05 # Add 5% buffer
global_raw_xlims <- c(0, global_max_raw_rank * axis_buffer); global_raw_ylims <- c(1, global_max_raw_count * axis_buffer)
global_filt_xlims <- c(0, global_max_filt_rank * axis_buffer); global_filt_ylims <- c(1, global_max_filt_count * axis_buffer)
cat("Global Raw Limits: Rank ~", global_max_raw_rank, ", Count ~", global_max_raw_count, "\n")
cat("Global Filt Limits: Rank ~", global_max_filt_rank, ", Count ~", global_max_filt_count, "\n")
cat(" -> Raw Axis Ranges: X=", paste(round(global_raw_xlims), collapse="-"), "; Y=", paste(round(global_raw_ylims), collapse="-"), "\n")
cat(" -> Filt Axis Ranges: X=", paste(round(global_filt_xlims), collapse="-"), "; Y=", paste(round(global_filt_ylims), collapse="-"), "\n")
cat("--- Global Limits Calculated ---\n\n"); rm(temp_filtered_counts)


# --- Process Each Sample Individually --- #
cat("--- Processing Individual Samples ---\n")
filtered_data_list <- list()
all_plots <- list(raw_dist = list(), raw_cumsum = list(), filt_dist = list(), filt_cumsum = list())
processed_samples <- names(starcode_data_list)

for (s_name in processed_samples) {
  cat("\nProcessing Sample:", s_name, "-------------------\n")
  # Pass the globally defined threshold from args
  processing_output <- parse.starcode(starcode_data = starcode_data_list[[s_name]], sample_id = s_name, count_threshold = raw_count_threshold,
                                      individual_plot_dir = individual_plots_dir, results_out_dir = results_dir,
                                      global_raw_xlims = global_raw_xlims, global_raw_ylims = global_raw_ylims,
                                      global_filt_xlims = global_filt_xlims, global_filt_ylims = global_filt_ylims)
  if (!is.null(processing_output$filtered_df)) filtered_data_list[[s_name]] <- processing_output$filtered_df
  if (!is.null(processing_output$plots$raw_dist)) all_plots$raw_dist[[s_name]] <- processing_output$plots$raw_dist
  if (!is.null(processing_output$plots$raw_cumsum)) all_plots$raw_cumsum[[s_name]] <- processing_output$plots$raw_cumsum
  if (!is.null(processing_output$plots$filt_dist)) all_plots$filt_dist[[s_name]] <- processing_output$plots$filt_dist
  if (!is.null(processing_output$plots$filt_cumsum)) all_plots$filt_cumsum[[s_name]] <- processing_output$plots$filt_cumsum
  cat("Finished processing:", s_name, "\n")
}

cat("\n--- Individual Sample Processing Summary ---\n"); cat("Attempted:", length(processed_samples), "; Successfully filtered:", length(filtered_data_list), "\n")
failed_filter_samples <- setdiff(processed_samples, names(filtered_data_list))
if (length(failed_filter_samples) > 0) cat("Samples with no filtered data:", paste(failed_filter_samples, collapse=", "), "\n")
cat("Individual plots/results saved to:", individual_plots_dir, "and", results_dir, "\n")
cat("--- End of Individual Sample Processing ---\n\n")
# rm(starcode_data_list); gc() # Optional cleanup

# --- Create Overlaid Plots --- #
cat("--- Generating Overlaid Plots ---\n")
# Pass the globally defined threshold from args for use in title
plot.overlaid.density(filtered_data_list, plot_title = paste("Overlaid Log10 Count Density (Filtered, Count >=", raw_count_threshold, ")"),
                      output_filename = file.path(combined_plots_dir, "Combined_06_Overlaid_Filtered_Density.png"))
plot.overlaid.cumsum(filtered_data_list, plot_title = paste("Overlaid Cumulative Sum (Filtered, Count >=", raw_count_threshold, ")"),
                     output_filename = file.path(combined_plots_dir, "Combined_07_Overlaid_Filtered_CumSum.png"))
cat("--- Finished Overlaid Plots ---\n\n")


# --- Combine Individual Plots (Patchwork) --- #
cat("--- Combining Individual Plots (Patchwork) ---\n")
num_samples_to_plot <- length(processed_samples)
ncols_combined <- if(num_samples_to_plot <= 4) 2 else if (num_samples_to_plot <= 9) 3 else 4

save_combined_plot(all_plots$raw_dist, "Combined_01_Raw_Barcode_Frequency", combined_plots_dir, "Raw Barcode Frequency (All Samples, Consistent Axes)", ncols = ncols_combined)
save_combined_plot(all_plots$raw_cumsum, "Combined_02_Raw_Cumulative_Sum", combined_plots_dir, "Raw Cumulative Sum (All Samples)", ncols = ncols_combined)
# Pass the globally defined threshold from args for use in title
save_combined_plot(all_plots$filt_dist, "Combined_03_Filtered_Barcode_Frequency", combined_plots_dir, paste("Filtered Barcode Frequency (Count >=", raw_count_threshold,", Consistent Axes)"), ncols = ncols_combined)
save_combined_plot(all_plots$filt_cumsum, "Combined_04_Filtered_Cumulative_Sum", combined_plots_dir, paste("Filtered Cumulative Sum (Count >=", raw_count_threshold, ")"), ncols = ncols_combined)
cat("--- Finished Combining Plots ---\n\n")


# --- Barcode Overlap Analysis and Visualization --- #
cat("--- Barcode Overlap Analysis ---\n")
if (length(filtered_data_list) < 2) {
  cat("Skipping overlap analysis: Less than 2 samples have filtered data.\n")
} else {
  cat("Analyzing overlap for filtered barcodes across", length(filtered_data_list), "samples...\n")
  cat("Samples included:", paste(names(filtered_data_list), collapse=", "), "\n")
  
  # Prepare barcode list, ensuring valid character vectors
  barcode_list_for_overlap <- lapply(filtered_data_list, function(df) {
    if("Barcode" %in% colnames(df)) {
      bc_vec <- as.character(df$Barcode); bc_vec <- bc_vec[!is.na(bc_vec) & nzchar(bc_vec)]; return(bc_vec)
    } else { return(character(0)) } })
  barcode_list_for_overlap <- barcode_list_for_overlap[sapply(barcode_list_for_overlap, length) > 0]
  
  if (length(barcode_list_for_overlap) < 2) {
    cat("Skipping overlap analysis: Less than 2 valid barcode lists remain after final checks.\n")
  } else {
    
    # Calculate Common/Unique (Independent of visualization)
    common_all_samples <- tryCatch(Reduce(intersect, barcode_list_for_overlap), error = function(e) { cat("Error calculating common barcodes:", e$message, "\n"); character(0) })
    num_common_all <- length(common_all_samples)
    cat("\nBarcodes common to ALL", length(barcode_list_for_overlap), "samples:", num_common_all, "\n")
    if (num_common_all > 0 && num_common_all < 50) { cat("Common barcodes (first few):", paste(head(common_all_samples, 10), collapse=", "), if(num_common_all > 10) "...", "\n") }
    unique_all_samples <- tryCatch(unique(unlist(barcode_list_for_overlap)), error = function(e) { cat("Error calculating unique barcodes:", e$message, "\n"); character(0) })
    num_unique_all <- length(unique_all_samples)
    cat("Total unique filtered barcodes across these samples:", num_unique_all, "\n")
    
    # Visualization: UpSet Plot or Venn Diagram
    num_sets <- length(barcode_list_for_overlap)
    
    # Use UpSet Plot for > 5 sets
    if (num_sets > 5) {
      cat("\nNumber of sets (", num_sets, ") > 5. Generating UpSet plot instead of Venn Diagram.\n")
      upset_plot_filename_pdf <- file.path(combined_plots_dir, "Combined_05_Filtered_Barcode_Overlap_UpSet.pdf")
      upset_plot_filename_png <- file.path(combined_plots_dir, "Combined_05_Filtered_Barcode_Overlap_UpSet.png")
      upset_data <- tryCatch(UpSetR::fromList(barcode_list_for_overlap), error = function(e){ cat("Error preparing data for UpSetR:", e$message, "\n"); NULL })
      
      if (!is.null(upset_data)){
        cat("Generating UpSet plot (saving as PDF and PNG)...\n")
        # Use the globally defined threshold from args for plot title
        main_title_overlap <- paste("Overlap (Filtered, Count >=", raw_count_threshold, ")")
        upset_params <- list(nsets = num_sets, nintersects = 40, order.by = "freq", decreasing = TRUE,
                             mb.ratio = c(0.6, 0.4), point.size = 2.5, line.size = 1,
                             mainbar.y.label = "Intersection Size", sets.x.label = "Set Size",
                             text.scale = c(1.3, 1.3, 1, 1, 1.2, 1.1), # Adjusted text scales
                             main.bar.color = "steelblue", sets.bar.color = "maroon", # Example colors
                             query.legend = "bottom") # Example query legend position
        
        # Save as PDF
        tryCatch({ pdf(upset_plot_filename_pdf, onefile = FALSE, width = 10, height = 6); print(do.call(upset, c(list(data=upset_data, queries=list(list(query=intersects, params=list(names(barcode_list_for_overlap)), color="red", active=T, query.name=paste0("Common to all ",num_sets)))), upset_params))); dev.off() # Added query for all common
          cat("Saved UpSet plot (PDF) to:", upset_plot_filename_pdf, "\n") }, error = function(e){ cat("Error saving UpSet plot (PDF):", e$message, "\n"); if(names(dev.cur()) != "null device") dev.off() })
        # Save as PNG
        tryCatch({ png(upset_plot_filename_png, width = 10*100, height = 6*100, res=100, type="cairo"); print(do.call(upset, c(list(data=upset_data, queries=list(list(query=intersects, params=list(names(barcode_list_for_overlap)), color="red", active=T, query.name=paste0("Common to all ",num_sets)))), upset_params))); dev.off() # Added query for all common
          cat("Saved UpSet plot (PNG) to:", upset_plot_filename_png, "\n") }, error = function(e){ cat("Error saving UpSet plot (PNG):", e$message, "\n"); if(names(dev.cur()) != "null device") dev.off() })
      } else { cat("Skipping UpSet plot generation due to data preparation error.\n") }
      
      # Use Venn Diagram for 2-5 sets
    } else if (num_sets >= 2 && num_sets <= 5) {
      venn_plot_filename <- file.path(combined_plots_dir, "Combined_05_Filtered_Barcode_Overlap_Venn.png")
      cat("\nGenerating Venn Diagram (", num_sets, " sets) to:", venn_plot_filename, "\n")
      venn_colors <- RColorBrewer::brewer.pal(max(3, num_sets), "Pastel1")[1:num_sets]
      # Use the globally defined threshold from args for plot title
      main_title_overlap <- paste("Overlap (Filtered, Count >=", raw_count_threshold, ")")
      venn_plot_object <- tryCatch(venn.diagram( x = barcode_list_for_overlap, filename = NULL, col = venn_colors, fill = venn_colors, alpha = 0.5, cex = 0.9, fontface = "plain", cat.cex = 1.0, cat.fontface = "bold", cat.col = venn_colors, cat.pos = 0, cat.dist = 0.06, margin = 0.1, disable.logging = TRUE, main = main_title_overlap, main.cex = 1.2, main.fontface = "bold"), error = function(e) { cat("ERROR generating Venn Diagram object:", e$message, "\n"); NULL })
      if (!is.null(venn_plot_object)) { png(filename = venn_plot_filename, width = 8, height = 8, units = "in", res = 300, type = "cairo"); grid.draw(venn_plot_object); dev.off(); cat("Venn Diagram saved successfully.\n") } else { cat("Venn Diagram could not be saved.\n") }
    } else { cat("\nUnexpected number of sets (", num_sets, "), skipping overlap visualization.\n") }
  }
}
cat("--- Finished Barcode Overlap Analysis ---\n\n")

#----- Analysis Complete -----#
cat("===============================================\n")
cat("SPLINTR Analysis Script Finished.\n")
cat("Outputs:\n")
cat(" -> Individual sample results (tables, fasta):", results_dir, "\n")
cat(" -> Individual sample plots:", individual_plots_dir, "\n")
cat(" -> Combined & Overlaid plots:", combined_plots_dir, "\n")
cat("===============================================\n")