library(iRfcb)

ifcb_path <- Sys.getenv("ifcb_path")

class_counts <- ifcb_count_mat_annotations(file.path(ifcb_path, "manual/Baltic"),
                                           file.path(ifcb_path, "config/class2use_Baltic.mat"),
                                           skip_class = 1)

write.table(class_counts,
            "output/class_counts.txt",
            na = "",
            row.names = FALSE,
            sep = "\t")
