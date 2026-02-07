# Given raw segment files,

# Directory name used as the label
target_label <- "Bald Eagle"

library(dplyr)

# Paths
root_dir <- "/Volumes/gioj_b2/OESF_processed/2023/segments"
dir_cache_out = "data/cache/helpers/dir_index.rds"

if (!file.exists(dir_cache_out)) {
  message("Finding all directories under the root. This may take some time...")
  dir.create(dirname(dir_cache_out), recursive = TRUE, showWarnings = FALSE)
  
  all_dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
  
  dir_index <- data.frame(
    path = all_dirs,
    name = basename(all_dirs),
    stringsAsFactors = FALSE
  )
  
  message("Caching result to ", dir_cache_out)
  saveRDS(dir_index, dir_cache_out)
  
} else {
  message("Loading cached directory index")
  dir_index = readRDS(dir_cache_out)
}
###################################

dest_dir <- "/Volumes/gioj/calibration/audio/OESF_samples/2023"
dest_dir = paste0(dest_dir, "/", target_label)

# 1. Find all directories with the target name
label_dirs = dir_index %>% filter(name == target_label) %>% pull(path)

# 2. List .wav files only inside those directories
message("Finding all .wav files for class '", target_label, "'")
wav_files <- unlist(
  lapply(label_dirs, function(d) {
    list.files(d, pattern = "\\.wav$", full.names = TRUE)
  }),
  use.names = FALSE
)

# 3. Copy files
message("Copying ", length(wav_files)," files to ", dest_dir)
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
result = file.copy(
  from = wav_files,
  to = dest_dir,
  overwrite = TRUE
)

if (all(result)) {
  message(crayon::green("Finished"))
} else {
  message(crayon::red("Error copying files"))
}
