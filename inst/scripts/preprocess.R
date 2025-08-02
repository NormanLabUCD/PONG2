#!/usr/bin/env Rscript

#Ensure the files are the right format
data_format <- function(geno, out=""){
    message("Reformasting the data to the right format")
    # Define source and destination directories
    #plink --bfile chr19 --chr 19 --from-bp 55235681 --to-bp 55248171 --make-bed --out kir3dl3
    plink_path <- Sys.which("plink")
    if (plink_path == "") {
    stop("PLINK not found in system's PATH. Please download and install PLINK.
         For more information, visit: https://www.cog-genomics.org/plink2")
    }
    
    src_dir <- geno
    dst_dir <- out
    

    # Read the .bim file
    bim_file <- paste0(src_dir, ".bim")
    if (!file.exists(bim_file)) {
    stop(paste("File not found:", bim_file))
    }
    bim_data <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)

  # Define the expected file name pattern
    chr19_pattern <- "^(chr)?19"
    #print(head(bim_data[,1]))
    # Check if the file names match the pattern
    if (!all(grepl(chr19_pattern, bim_data[, 1], ignore.case = TRUE))) {
    ("Error: Only files with chromosome 19 are allowed.")
    }

    # Create the new SNP ID column (CHR:POS)
    bim_data$V2 <- paste0(bim_data$V1, ":", bim_data$V4)

    # Select and reorder the columns
    #formatted_bim <- bim_data[, c(1, 6, 3, 4, 5, 2)]

    # Write the formatted data to a new .bim file
    write.table(bim_data, paste0(dst_dir, "/", basename(src_dir), ".bim"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    # Identify and write duplicated SNPs to a file
    message("removing duplicated SNPs")
    duplicated_snps <- bim_data[duplicated(bim_data[, c(1, 4)]), ]
    write.table(duplicated_snps, paste0(dst_dir, "/snp_list.txt"),
                row.names = FALSE, col.names = TRUE, quote = FALSE)

  # Loop through each file extension and copy the file
   exts <- c(".bed", ".fam")
  for (ext in exts) {
    file_name <- paste0(src_dir, ext)
    system(paste("cp", file_name, dst_dir))
  }

    # Display the number of duplicated SNPs
    message(paste(dim(duplicated_snps)[1], "duplicate were found."))
  # Specify the path to the PLINK executable
    #plink <- system.file(package = "PONG", "bin", "plink")

  # Run PLINK to exclude duplicated SNPs and create a new set of files
    system(paste("plink --bfile", paste0(dst_dir, "/", basename(src_dir)),
               "--exclude", paste0(dst_dir, "/", "snp_list.txt"),
               "--make-bed --out", paste0(dst_dir, "/", "chr19")))
  message("Completed data reformatting.")
  message("Output files are saved in:", dst_dir)

  return(paste0(dst_dir, "/", "chr19"))
  
}




#Ensure users provide the require parameters for the model to work
parameters <- function(){


}